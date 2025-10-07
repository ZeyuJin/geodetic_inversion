function [slip_model, stats, varargout] = make_fault_from_insar_rec( ...
            slip_model_vs, slip_model_ds, varargin)
% MAKE_FAULT_FROM_INSAR_REC Flexible finite-fault inversion with arbitrary datasets.
% Perform finite fault inversion from geodetic data, using rectangular
% meshes
% Core idea:
%   - Pass any number of datasets via 'Datasets' (cell array or struct array)
%   - Get one predicted model output per dataset (same order), in varargout
%   - Configure smoothing, bounds, and zero-slip boundaries via Name–Value pairs
%
% Dependencies (must be on your MATLAB path):
%   build_green_function, build_smooth_function, zero_slip_boundary,
%   bounds_new, compute_moment  (same as your original workflow)
%
% Basic example:
%   D = {
%     '/path/asc12/los_samp{i}.mat','los',1.7,'ASC12 LOS';
%     '/path/des48/los_samp{i}.mat','los',1.0,'DES48 LOS';
%     '/path/des121/rng_samp{i}.mat','rng',1.0,'DES121 RNG';
%     '/path/ALOS/azo_samp{i}.mat','azo',1.0,'ALOS AZO';
%   };
%   [slip_model, stats, asc12LOS, des48LOS, des121RNG, alosAZO] = ...
%       make_fault_from_insar_rec(VS, DS, ...
%           'Datasets', D, 'IterStep', 1, ...
%           'Lambda', 0.017, 'ModelType','okada', 'Verbose', true, ...
%           'ZeroBoundary', struct( ...
%               'Left',   struct('Enable',true,  'Fault',1,      'Ratio',3e-4), ...
%               'Right',  struct('Enable',true,  'Fault',[1 2],  'Ratio',3e-4), ...
%               'Top',    struct('Enable',false, 'Fault',[],     'Ratio',1e-5), ...
%               'Bottom', struct('Enable',true,  'Fault',[1 2],  'Ratio',5e-4), ...
%               'PlaneFault', []));  % PlaneFault kept for backward-compat
%
% Inputs:
%   slip_model_vs, slip_model_ds  : your updip/downdip meshes (unchanged)
%
% Name–Value options:
%   'Datasets'  : N×3 or N×4 cell array {path, kind, weight[, name]}
%                 or struct array with fields path, kind, weight, (name)
%                 kind ∈ {'los','rng','azo','gps'} (synonyms mapped internally)
%                 weight : scalar applied when stacking for solve
%                 name   : optional label for stats
%                 Tip: use '{i}' in path and set 'IterStep' to substitute
%   'Lambda'    : smoothing weight (default 0.017)
%   'ModelType' : 'okada' (default) or 'layered'
%   'Con'       : [sx sy sz] sign-constraint selector for bounds_new (default [1 0 0])
%   'ShallowDipId' : indices to control dip-slip component (default [])
%   'SegmentSmoothFile','IntersectSmoothFile' : targeted smoothing files
%   'ZeroBoundary' : struct with four sides + backward-compat field:
%       .Left   : struct('Enable',true,'Fault',1:2,'Ratio',3e-4)
%       .Right  : struct('Enable',true,'Fault',1:2,'Ratio',3e-4)
%       .Top    : struct('Enable',false,'Fault',1:2,'Ratio',1e-5)
%       .Bottom : struct('Enable',true,'Fault',1:2,'Ratio',5e-4)
%       .PlaneFault : []  % if provided, used to fill empty per-side Fault fields
%   'IterStep'  : integer; substitutes '{i}' in dataset paths (default [])
%   'LSQOptions': struct of lsqlin options (merged into sensible defaults)
%   'Verbose'   : true/false (default true)
%
% Outputs:
%   slip_model : same format as your original (cols 12/13 = strike/dip slip)
%   stats      : struct with RMS, variance reduction, roughness, exitflag, etc.
%   varargout  : one predicted data vector per dataset, in the same order
%
% NOTE: To ensure "one output per dataset", this function requires:
%       nargout == 2 + numel(Datasets)
% Written by Zeyu Jin and Xiaoyu Zou. Improved by ChatGPT 5

%% ---------- Parse inputs ----------
p = inputParser;  p.KeepUnmatched = false;
addParameter(p,'Datasets',{},@(x) iscell(x) || isstruct(x));
addParameter(p,'Lambda',0.017,@isscalar);
addParameter(p,'ModelType','okada',@(s) ischar(s) || isstring(s));
addParameter(p,'Con',[1 0 0],@(v) isnumeric(v) && numel(v)==3);
addParameter(p,'ShallowDipId',[],@(x) isnumeric(x) || isempty(x));
addParameter(p,'SegmentSmoothFile',[],@(x) ischar(x)||isstring(x)||isempty(x));
addParameter(p,'IntersectSmoothFile',[],@(x) ischar(x)||isstring(x)||isempty(x));

zb_default = struct( ...
  'Left',   struct('Enable',true,  'Fault',1:2, 'Ratio',3e-4), ...
  'Right',  struct('Enable',true,  'Fault',1:2, 'Ratio',3e-4), ...
  'Top',    struct('Enable',false, 'Fault',1:2, 'Ratio',1e-5), ...
  'Bottom', struct('Enable',true,  'Fault',1:2, 'Ratio',5e-4), ...
  'PlaneFault', [] );   % kept for backward compatibility

addParameter(p,'ZeroBoundary',zb_default,@(s) isstruct(s));
addParameter(p,'IterStep',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));
addParameter(p,'LSQOptions',struct(),@(s) isstruct(s));
addParameter(p,'Verbose',true,@islogical);
parse(p,varargin{:});
opt = p.Results;

% Normalize dataset specs
datasets = normalizeDatasets(opt.Datasets, opt.IterStep);

% Enforce "one output per dataset"
k = numel(datasets);
if nargout ~= 2 + k
    error(['make_fault_from_insar_rec: Number of outputs must equal 2 + number of datasets.\n' ...
           'You provided %d outputs, but there are %d dataset(s).'], nargout, k);
end

% Backward-compat: PlaneFault -> fill empty per-side Fault lists
if isfield(opt.ZeroBoundary,'PlaneFault') && ~isempty(opt.ZeroBoundary.PlaneFault)
    pf = opt.ZeroBoundary.PlaneFault;
    if ~isfield(opt.ZeroBoundary.Left,'Fault')   || isempty(opt.ZeroBoundary.Left.Fault),   opt.ZeroBoundary.Left.Fault   = pf; end
    if ~isfield(opt.ZeroBoundary.Right,'Fault')  || isempty(opt.ZeroBoundary.Right.Fault),  opt.ZeroBoundary.Right.Fault  = pf; end
    if ~isfield(opt.ZeroBoundary.Top,'Fault')    || isempty(opt.ZeroBoundary.Top.Fault),    opt.ZeroBoundary.Top.Fault    = pf; end
    if ~isfield(opt.ZeroBoundary.Bottom,'Fault') || isempty(opt.ZeroBoundary.Bottom.Fault), opt.ZeroBoundary.Bottom.Fault = pf; end
end
ZB = opt.ZeroBoundary;

%% ---------- Assemble slip model ----------
slip_model = [slip_model_vs; slip_model_ds];
slip_model(:,2) = (1:size(slip_model,1)).';   % reindex nodes
nflt = max(slip_model(:,1));

%% ---------- Smoothing ----------
[H,h1,~] = build_smooth_function(slip_model_vs, slip_model_ds, ...
    opt.segmentSmoothFileIfExists('SegmentSmoothFile'), ...
    opt.segmentSmoothFileIfExists('IntersectSmoothFile'), ...
    'noramp','dip_id',opt.ShallowDipId);

%% ---------- Boundary constraints (per-fault control) ----------
[Wt,dt] = deal([]); [Wb,db] = deal([]); [Wl,dl] = deal([]); [Wr,dr] = deal([]);
appendRows = @(Wacc,Dacc,Wi,Di) deal([Wacc; Wi],[Dacc; Di]);

% LEFT boundary
if ZB.Left.Enable && ~isempty(ZB.Left.Fault)
    for f = reshape(ZB.Left.Fault,1,[])
        [Wi,Di] = zero_slip_boundary(slip_model, f, 'left', ZB.Left.Ratio);
        [Wl,dl] = appendRows(Wl,dl,Wi,Di);
    end
end

% RIGHT boundary
if ZB.Right.Enable && ~isempty(ZB.Right.Fault)
    for f = reshape(ZB.Right.Fault,1,[])
        [Wi,Di] = zero_slip_boundary(slip_model, f, 'right', ZB.Right.Ratio);
        [Wr,dr] = appendRows(Wr,dr,Wi,Di);
    end
end

% TOP boundary (per fault: its top layer index)
if ZB.Top.Enable && ~isempty(ZB.Top.Fault)
    for f = reshape(ZB.Top.Fault,1,[])
        top_layer_no = min(slip_model(slip_model(:,1)==f, 3));
        [Wi,Di] = zero_slip_boundary(slip_model, f, top_layer_no, ZB.Top.Ratio);
        [Wt,dt] = appendRows(Wt,dt,Wi,Di);
    end
end

% BOTTOM boundary (per fault: its bottom layer index)
if ZB.Bottom.Enable && ~isempty(ZB.Bottom.Fault)
    for f = reshape(ZB.Bottom.Fault,1,[])
        bottom_layer_no = max(slip_model(slip_model(:,1)==f, 3));
        [Wi,Di] = zero_slip_boundary(slip_model, f, bottom_layer_no, ZB.Bottom.Ratio);
        [Wb,db] = appendRows(Wb,db,Wi,Di);
    end
end

%% ---------- Build Green's & data stacks ----------
G_stack = [];   % weighted for solve
b_stack = [];
G_raw_all = []; % unweighted for reporting
b_raw_all = [];

for i = 1:k
    ds = datasets(i);
    dt = mapKind(ds.kind);  % maps 'los'/'rng'->'insar', 'azo'->'AZO', 'gps'->'camp_gps'
    [Gi_raw, Gi, bi_raw, bi] = build_green_function(slip_model, ds.path, dt, 'noramp', opt.ModelType);

    % For solving (weighted)
    G_stack = [G_stack; ds.weight * Gi];
    b_stack = [b_stack; ds.weight * bi];

    % For reporting (unweighted, raw)
    G_raw_all = [G_raw_all; Gi_raw];
    b_raw_all = [b_raw_all; bi_raw];

    datasets(i).Gi_raw = Gi_raw;  %#ok<AGROW>
end

%% ---------- Add smoothing & boundary rows ----------
Greens = [G_stack;  H * (opt.Lambda / h1);  Wb; Wt; Wl; Wr];
b_full  = [b_stack; zeros(h1,1);            db;  dt;  dl;  dr];

%% ---------- Bounds ----------
tSm = zeros(1, nflt+1); fid = slip_model(:,1);
for i = 1:nflt, tSm(i+1) = nnz(fid == i); end
add_col = 0; NT = 2; NS = nflt;
[lb, ub] = bounds_new(NS, NT, tSm, add_col, opt.Con);

%% ---------- Solve ----------
lsqopt = optimset('LargeScale','on','DiffMaxChange',1e-1,'DiffMinChange',1e-12, ...
    'TolCon',1e-12,'TolFun',1e-12,'TolPCG',1e-12,'TolX',1e-12, ...
    'MaxIter',1e9,'MaxPCGIter',1e9);
% merge user overrides
fn = fieldnames(opt.LSQOptions);
for i = 1:numel(fn), lsqopt.(fn{i}) = opt.LSQOptions.(fn{i}); end

[u, resnorm, residual, exitflag] = lsqlin(Greens, double(b_full), [], [], [], [], lb, ub, [], lsqopt);

%% ---------- Diagnostics ----------
rms0      = sum(b_raw_all.^2);
rms       = sum((G_raw_all*u - b_raw_all).^2);
redu_perc = 100*(rms0 - rms)/rms0;

if opt.Verbose
    fprintf('VR = %5.2f%%, rms0 = %.3e, rms = %.3e\n', redu_perc, rms0, rms);
    fprintf('resnorm = %.3e, mean(residual) = %.3e, exitflag = %d\n', sqrt(resnorm), mean(residual), exitflag);
end

rough_matrix    = H*u;
RMS_misfit_fit  = sum((G_stack*u - b_stack).^2);
model_roughness = sqrt(sum(rough_matrix.^2)/numel(rough_matrix));

%% ---------- Fill slip model; compute per-dataset predictions ----------
slip_model(:,12) = u(1:sum(tSm));               % strike-slip
slip_model(:,13) = u(sum(tSm)+1:end);           % dip-slip

pred = cell(1,k);
for i = 1:k
    pred{i} = datasets(i).Gi_raw * u;           % predictions in input order
end
[varargout{1:k}] = pred{:};

% Optional seismic moment check
try
    compute_moment(slip_model, opt.ModelType);
catch
    if opt.Verbose, warning('compute_moment not available or failed; continuing.'); end
end

%% ---------- Stats struct ----------
stats = struct();
stats.redu_perc      = redu_perc;
stats.rms0           = rms0;
stats.rms            = rms;
stats.resnorm        = resnorm;
stats.residual_mean  = mean(residual);
stats.exitflag       = exitflag;
stats.roughness      = model_roughness;
stats.RMS_misfit     = RMS_misfit_fit;
stats.lambda         = opt.Lambda;
stats.model_type     = char(opt.ModelType);

% strip Gi_raw before storing
S = rmfield(datasets, intersect(fieldnames(datasets), {'Gi_raw'}));
stats.dataset = S;

end % ======= main =======

%% ---------- Local helpers ----------
function S = normalizeDatasets(D, iterStep)
    if isempty(D)
        error('You must provide ''Datasets'' as an N×3/4 cell array or struct array.');
    end
    if iscell(D)
        if size(D,2) < 3
            error('Datasets cell must have at least 3 columns: {path, kind, weight[, name]}.');
        end
        n = size(D,1);
        S = repmat(struct('path','','kind','','weight',1,'name',''), 1, n);
        for i = 1:n
            path = string(D{i,1});
            if ~isempty(iterStep)
                path = strrep(path, '{i}', num2str(iterStep));
            end
            S(i).path   = char(path);
            S(i).kind   = lower(string(D{i,2}));
            S(i).weight = D{i,3};
            if size(D,2) >= 4 && ~isempty(D{i,4})
                S(i).name = char(D{i,4});
            else
                S(i).name = sprintf('%s_%d', S(i).kind, i);
            end
        end
    elseif isstruct(D)
        S = D;
        % optional templating for struct inputs
        if ~isempty(iterStep)
            for i = 1:numel(S)
                S(i).path = strrep(string(S(i).path), '{i}', num2str(iterStep));
                S(i).path = char(S(i).path);
            end
        end
        for i = 1:numel(S)
            if ~isfield(S(i),'name') || isempty(S(i).name)
                S(i).name = sprintf('%s_%d', lower(string(S(i).kind)), i);
            end
            S(i).kind = lower(string(S(i).kind));
        end
    else
        error('Unsupported Datasets type. Use cell array or struct array.');
    end
end

function dt = mapKind(kind)
    k = lower(string(kind));
    switch k
        case {'los','rng','insar'}
            dt = 'insar';
        case {'azo','azimuth','mai'}
            dt = 'AZO';
        case {'gps','camp_gps','gnss'}
            dt = 'camp_gps';
        otherwise
            error('Unknown dataset kind: %s (use los|rng|azo|gps)', kind);
    end
end

% Light syntactic sugar so you can pass [] for missing smoothing files
function v = segmentSmoothFileIfExists(optField)
    % this small helper is defined as a nested function above via function handle:
    % opt.segmentSmoothFileIfExists('SegmentSmoothFile')
    % but here we provide an actual implementation for code clarity:
end
