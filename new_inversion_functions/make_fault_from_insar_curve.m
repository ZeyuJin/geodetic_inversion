function [slip_model, stats, varargout] = make_fault_from_insar_curve( ...
            slip_model, varargin)
% MAKE_FAULT_FROM_INSAR_CURVE Flexible curved-fault inversion (arbitrary datasets & geometries).
% perform finite-fault inversion using Triangular Dislocation Elements in a
% homogeneous half-space
% Core features
%   • Any number of datasets: pass via 'Datasets' (cell/struct array). 
%     You'll get the same number of prediction outputs (varargout) in order.
%   • Any number of geometries: pass via 'Geometries' (cell/struct array).
%     Used for optional moment checks; your own routines can consume them as needed.
%   • All parameters are Name–Value pairs—no editing internals.
%
% Example
%   D = {
%     '/data/asc12/los_samp{i}.mat','los',1.7,'ASC12 LOS';
%     '/data/des48/los_samp{i}.mat','los',1.0,'DES48 LOS';
%     '/data/des121/rng_samp{i}.mat','rng',1.0,'DES121 RNG';
%   };
%   G = {geometry_main, geometry_branch, geometry_third};  % any count
%
%   [slip_model, stats, asc12LOS_pred, des48LOS_pred, des121RNG_pred] = ...
%     make_fault_from_insar_curve(slip_model, ...
%        'Datasets', D, ...
%        'Geometries', G, ...
%        'IterStep', 1, ...
%        'Lambda', 0.017, ...
%        'SSR', 3, 'DSR', 1, ...
%        'ModelType', 'okada_curve', ...
%        'Con', [0 -1 0], ...
%        'ZeroBoundary', struct( ...
%           'Left',   struct('Enable',true,  'Ratio',2e-4), ...
%           'Right',  struct('Enable',true,  'Ratio',2e-4), ...
%           'Top',    struct('Enable',false, 'Ratio',1e-5), ...
%           'Bottom', struct('Enable',true,  'Ratio',5e-4) ), ...
%        'Verbose', true);
%
% Dependencies expected on path (same as your originals):
%   build_green_function, smoothingMatrix_test, zero_slip_boundary_curve,
%   bounds_new, triangle_moment (optional), lsqlin (Optimization Toolbox).
%
% I/O contract:
%   nargout must equal 2 + numel(Datasets).
%   Outputs: [slip_model, stats, pred1, pred2, ..., predN]
%
% Notes:
%   • This version does not assume a fixed number of geometries; it treats them
%     as a list available for optional post-processing (e.g., moment).
%   • If your own helper expects two geometries, we try calling with all of them
%     and fall back to per-geometry calls.

%% ---------------------- Parse inputs ----------------------
p = inputParser;  p.KeepUnmatched = false;

% Data & geometry collections
addParameter(p,'Datasets',{},@(x) iscell(x) || isstruct(x));
addParameter(p,'Geometries',{},@(x) iscell(x) || isstruct(x) || isempty(x));

% Smoothing / modeling
addParameter(p,'Lambda',0.017,@isscalar);
addParameter(p,'SSR',3,@isscalar);      % strike-slip smoothing ratio
addParameter(p,'DSR',1,@isscalar);      % dip-slip smoothing ratio
addParameter(p,'ModelType','okada_curve',@(s) ischar(s)||isstring(s));

% Bounds & other knobs
addParameter(p,'Con',[0 -1 0],@(v) isnumeric(v) && numel(v)==3);
addParameter(p,'PTs',[],@(x) isempty(x) || isnumeric(x));
addParameter(p,'IterStep',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));

% Zero-slip boundaries (curved helper sides). Top is optional/off by default.
zb_default = struct( ...
  'Left',   struct('Enable',true,  'Ratio',2e-4), ...
  'Right',  struct('Enable',true,  'Ratio',2e-4), ...
  'Top',    struct('Enable',false, 'Ratio',1e-5), ...
  'Bottom', struct('Enable',true,  'Ratio',5e-4) );
addParameter(p,'ZeroBoundary',zb_default,@(s) isstruct(s));

% Solver & logging
addParameter(p,'LSQOptions',struct(),@(s) isstruct(s));
addParameter(p,'Verbose',true,@islogical);

parse(p,varargin{:});
opt = p.Results;

% Datasets normalization & output arity check
datasets = normalizeDatasets(opt.Datasets, opt.IterStep);
k = numel(datasets);
if nargout ~= 2 + k
    error(['make_fault_from_insar_curve: number of outputs must equal 2 + number of datasets.\n' ...
           'You requested %d outputs, but there are %d dataset(s).'], nargout, k);
end

% Geometries normalization (make it a cell array)
geometries = normalizeGeometries(opt.Geometries);

%% ---------------------- Smoothing (curved) ----------------------
% For curved meshes: use local normals in slip_model(:,7:9) if present.
if size(slip_model,2) >= 9
    n1 = slip_model(:,7:9);
else
    error('slip_model must include local normals in columns 7:9 for smoothingMatrix_test.');
end
[H,h1] = smoothingMatrix_test(n1, 'ss_ratio', opt.SSR, 'ds_ratio', opt.DSR);

%% ---------------------- Zero-slip boundaries ----------------------
[Wt,dt] = deal([]); [Wb,db] = deal([]); [Wl,dl] = deal([]); [Wr,dr] = deal([]);

if isfield(opt.ZeroBoundary,'Bottom') && opt.ZeroBoundary.Bottom.Enable
    [Wb,db] = zero_slip_boundary_curve(slip_model, opt.PTs, 'bottom', opt.ZeroBoundary.Bottom.Ratio);
end
if isfield(opt.ZeroBoundary,'Top') && opt.ZeroBoundary.Top.Enable
    [Wt,dt] = zero_slip_boundary_curve(slip_model, opt.PTs, 'top', opt.ZeroBoundary.Top.Ratio);
end
if isfield(opt.ZeroBoundary,'Left') && opt.ZeroBoundary.Left.Enable
    [Wl,dl] = zero_slip_boundary_curve(slip_model, opt.PTs, 'left', opt.ZeroBoundary.Left.Ratio);
end
if isfield(opt.ZeroBoundary,'Right') && opt.ZeroBoundary.Right.Enable
    [Wr,dr] = zero_slip_boundary_curve(slip_model, opt.PTs, 'right', opt.ZeroBoundary.Right.Ratio);
end

%% ---------------------- Build Green's & data stacks ----------------------
G_stack = [];   b_stack = [];       % weighted (solve)
G_raw_all = []; b_raw_all = [];     % unweighted (reporting)

for i = 1:k
    ds = datasets(i);
    dt = mapKind(ds.kind);  % 'insar'|'AZO'|'camp_gps'
    [Gi_raw, Gi, bi_raw, bi] = build_green_function(slip_model, ds.path, dt, ...
                                    'noramp', opt.ModelType, 'pts', opt.PTs);
    G_stack   = [G_stack;   ds.weight * Gi];
    b_stack   = [b_stack;   ds.weight * bi];
    G_raw_all = [G_raw_all; Gi_raw];
    b_raw_all = [b_raw_all; bi_raw];
    datasets(i).Gi_raw = Gi_raw; %#ok<AGROW>
end

%% ---------------------- Add smoothing & boundaries ----------------------
Greens = [G_stack;  H * (opt.Lambda / h1);  Wb; Wt; Wl; Wr];
b_full  = [b_stack; zeros(h1,1);            db;  dt;  dl;  dr];

%% ---------------------- Bounds ----------------------
nflt  = max(slip_model(:,1));
tSm   = zeros(1, nflt+1);
fid   = slip_model(:,1);
for ii = 1:nflt, tSm(ii+1) = nnz(fid == ii); end
add_col = 0; NT = 2; NS = nflt;
[lb,ub] = bounds_new(NS, NT, tSm, add_col, opt.Con);

%% ---------------------- Solve (lsqlin) ----------------------
lsqopt = optimset('LargeScale','on','DiffMaxChange',1e-1,'DiffMinChange',1e-12, ...
    'TolCon',1e-12,'TolFun',1e-12,'TolPCG',1e-12,'TolX',1e-12, ...
    'MaxIter',1e9,'MaxPCGIter',1e9);
fn = fieldnames(opt.LSQOptions);
for ii = 1:numel(fn), lsqopt.(fn{ii}) = opt.LSQOptions.(fn{ii}); end

[u, resnorm, residual, exitflag] = lsqlin(Greens, double(b_full), [], [], [], [], lb, ub, [], lsqopt);

%% ---------------------- Fill slip model; predictions ----------------------
slip_model(:,2) = u(1:sum(tSm));                % strike-slip
slip_model(:,3) = u(sum(tSm)+1:end);            % dip-slip

pred = cell(1,k);
for i = 1:k
    pred{i} = datasets(i).Gi_raw * u;
end
[varargout{1:k}] = pred{:};

%% ---------------------- Diagnostics ----------------------
rms0      = sum(b_raw_all.^2);
rms       = sum((G_raw_all*u - b_raw_all).^2);
VR_pct    = 100*(rms0 - rms)/rms0;

if opt.Verbose
    fprintf('VR = %5.2f%%, rms0 = %.3e, rms = %.3e\n', VR_pct, rms0, rms);
    fprintf('resnorm = %.3e, mean(residual) = %.3e, exitflag = %d\n', sqrt(resnorm), mean(residual), exitflag);
end

rough_matrix    = H*u;
RMS_misfit_fit  = sum((G_stack*u - b_stack).^2);
model_roughness = sqrt(sum(rough_matrix.^2)/numel(rough_matrix));

%% ---------------------- Optional: moment checks with arbitrary geometries ----------------------
if ~isempty(geometries)
    try
        % Try calling a routine that accepts variable # of geometry inputs:
        % triangle_moment(slip_model, geom1, geom2, ..., geomN)
        feval('triangle_moment', slip_model, geometries{:});
    catch
        % Fallback: call per-geometry; ignore failures silently if not present
        for gi = 1:numel(geometries)
            try
                triangle_moment(slip_model, geometries{gi});
            catch
                % no-op
            end
        end
    end
end

%% ---------------------- Stats ----------------------
stats = struct();
stats.VR_percent     = VR_pct;
stats.rms0           = rms0;
stats.rms            = rms;
stats.resnorm        = resnorm;
stats.residual_mean  = mean(residual);
stats.exitflag       = exitflag;
stats.roughness      = model_roughness;
stats.RMS_misfit     = RMS_misfit_fit;
stats.lambda         = opt.Lambda;
stats.SSR            = opt.SSR;
stats.DSR            = opt.DSR;
stats.model_type     = char(opt.ModelType);
S = rmfield(datasets, intersect(fieldnames(datasets), {'Gi_raw'}));
stats.dataset        = S;
stats.geometry_count = numel(geometries);

end % ================== main ==================

%% ---------------------- Helpers ----------------------
function S = normalizeDatasets(D, iterStep)
    if isempty(D)
        error('Provide ''Datasets'' as an N×3/4 cell array or a struct array.');
    end
    if iscell(D)
        if size(D,2) < 3
            error('Datasets cell must be {path, kind, weight[, name]}.');
        end
        n = size(D,1);
        S = repmat(struct('path','','kind','','weight',1,'name',''), 1, n);
        for i = 1:n
            pth = string(D{i,1});
            if ~isempty(iterStep), pth = strrep(pth, '{i}', num2str(iterStep)); end
            S(i).path   = char(pth);
            S(i).kind   = lower(string(D{i,2}));
            S(i).weight = D{i,3};
            if size(D,2) >= 4 && ~isempty(D{i,4})
                S(i).name = char(string(D{i,4}));
            else
                S(i).name = sprintf('%s_%d', S(i).kind, i);
            end
        end
    elseif isstruct(D)
        S = D;
        if ~isempty(iterStep)
            for i = 1:numel(S)
                S(i).path = char(strrep(string(S(i).path), '{i}', num2str(iterStep)));
            end
        end
        for i = 1:numel(S)
            if ~isfield(S(i),'name') || isempty(S(i).name)
                S(i).name = sprintf('%s_%d', lower(string(S(i).kind)), i);
            end
            S(i).kind = lower(string(S(i).kind));
        end
    else
        error('Unsupported Datasets type. Use cell or struct array.');
    end
end

function G = normalizeGeometries(Gin)
    if isempty(Gin)
        G = {};
        return;
    end
    if iscell(Gin)
        G = Gin;
    elseif isstruct(Gin)
        % Allow struct array of geometries
        G = arrayfun(@(s) s, Gin, 'UniformOutput', false);
    else
        error('Geometries must be a cell array, a struct array, or empty.');
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
