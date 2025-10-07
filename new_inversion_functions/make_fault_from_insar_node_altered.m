function [slip_model, varargout] = make_fault_from_insar_node_altered(slip_model, varargin)
%MAKE_FAULT_FROM_INSAR_NODE_ALTERED  Flexible node-based finite-fault inversion.
% Perform finite fault inversion using the hybrid TDE-layered curved nodal
% model
%
% NEW MODE (recommended):
%   D = {
%     '/path/asc12/los_samp{i}.mat','los',1.7,'ASC12 LOS';
%     '/path/des48/los_samp{i}.mat','los',1.0,'DES48 LOS';
%     '/path/alos/rng_samp{i}.mat','rng',1.0,'ALOS RNG';
%   };
%   [slip_model, stats, pred1, pred2, pred3] = make_fault_from_insar_node_altered( ...
%       slip_model, 'Datasets', D, 'IterStep', 1, 'Lambda', 0.017, ...
%       'Smoothing', struct('Form','curvature','Mass','voronoi', ...
%                           'ss_ratio',3,'ds_ratio',1,'PerFaultNorm',true,'Verbose',true), ...
%       'ZeroBoundary', struct( ...
%           'Left',   struct('Enable',true,  'Fault',1:2,'Ratio',3e-4), ...
%           'Right',  struct('Enable',true,  'Fault',1:2,'Ratio',3e-4), ...
%           'Top',    struct('Enable',false, 'Fault',1:2,'Ratio',1e-5), ...
%           'Bottom', struct('Enable',true,  'Fault',1:2,'Ratio',5e-4), ...
%           'PlaneFault', []), ...
%       'Con',[0 -1 0], 'Builder','altered', 'Verbose',true);
%
% LEGACY MODE (backward compatibility):
%   If you request exactly 12 outputs, results are packed as:
%     [slip_model, RMS_misfit, model_roughness, insar_model1..8, rms]
%   (Predictions beyond 8 are ignored in legacy packing but still influence the inversion.)
%
% REQUIRED HELPERS ON PATH (as in your original workflow):
%   build_green_function_altered (or build_green_function), smoothingMatrix_laplace4,
%   zero_slip_boundary, bounds_new, lsqlin (Optimization Toolbox).

%% ---------------- Parse inputs ----------------
p = inputParser;  p.KeepUnmatched = false;
addParameter(p,'Datasets',{},@(x) iscell(x) || isstruct(x));
addParameter(p,'IterStep',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));
addParameter(p,'Lambda',0.017,@isscalar);                % original default 0.17e-1
addParameter(p,'Smoothing',struct(),@(s) isstruct(s));
addParameter(p,'ZeroBoundary',struct(),@(s) isstruct(s));
addParameter(p,'Con',[0 -1 0],@(v) isnumeric(v)&&numel(v)==3);
addParameter(p,'Builder','altered',@(s) any(strcmpi(s,{'altered','standard'})));
addParameter(p,'ModelType','okada_curve',@(s) ischar(s)||isstring(s));
addParameter(p,'PTs',[],@(x) isempty(x) || isnumeric(x));
addParameter(p,'LSQOptions',struct(),@(s) isstruct(s));
addParameter(p,'Verbose',true,@islogical);
% Back-compat alias
addParameter(p,'smoothness',[],@(x) isempty(x) || isscalar(x));
parse(p,varargin{:});
opt = p.Results;
if ~isempty(opt.smoothness), opt.Lambda = opt.smoothness; end

% Normalize datasets
datasets = normalizeDatasets(opt.Datasets, opt.IterStep);
k = numel(datasets);

% Default smoothing & boundaries
Sm = struct('Form','curvature','Mass','voronoi','ss_ratio',3,'ds_ratio',1, ...
            'PerFaultNorm',true,'Verbose',true);
Sm = mergeStructs(Sm, opt.Smoothing);

ZB = struct( ...
  'Left',   struct('Enable',true,  'Fault',1:2,'Ratio',3e-4), ...
  'Right',  struct('Enable',true,  'Fault',1:2,'Ratio',3e-4), ...
  'Top',    struct('Enable',false, 'Fault',1:2,'Ratio',1e-5), ...
  'Bottom', struct('Enable',true,  'Fault',1:2,'Ratio',5e-4), ...
  'PlaneFault', [] );
ZB = mergeStructs(ZB, opt.ZeroBoundary);
% Back-compat: if PlaneFault is supplied, fill missing per-side Fault lists
if isfield(ZB,'PlaneFault') && ~isempty(ZB.PlaneFault)
    if ~isfield(ZB.Left,'Fault')   || isempty(ZB.Left.Fault),   ZB.Left.Fault   = ZB.PlaneFault;   end
    if ~isfield(ZB.Right,'Fault')  || isempty(ZB.Right.Fault),  ZB.Right.Fault  = ZB.PlaneFault;   end
    if ~isfield(ZB.Top,'Fault')    || isempty(ZB.Top.Fault),    ZB.Top.Fault    = ZB.PlaneFault;   end
    if ~isfield(ZB.Bottom,'Fault') || isempty(ZB.Bottom.Fault), ZB.Bottom.Fault = ZB.PlaneFault;   end
end

%% ---------------- Smoothing (laplace4 with contact + normals) ----------------
% Original convention: normals in cols 6:8, contact matrix in cols 9:14.
lambda    = opt.Lambda * 1e-4;       % same scale as original
normals   = slip_model(:,6:8);
contactMat= slip_model(:,9:14);

[H,h1,~,~,~,~] = smoothingMatrix_laplace4( ...
    contactMat, normals, ...
    'Form',Sm.Form, 'Mass',Sm.Mass, ...
    'ss_ratio',Sm.ss_ratio, 'ds_ratio',Sm.ds_ratio, ...
    'PerFaultNorm', Sm.PerFaultNorm, ...
    'FaultID', slip_model(:,1), ...
    'Verbose', Sm.Verbose);

%% ---------------- Zero-slip boundaries (per side & per fault) ----------------
[Wt,dt] = deal([]); [Wb,db] = deal([]); [Wl,dl] = deal([]); [Wr,dr] = deal([]);
appendRows = @(Wacc,Dacc,Wi,Di) deal([Wacc; Wi],[Dacc; Di]);

% TOP per fault (layer index = min layer on that fault)
if ZB.Top.Enable && ~isempty(ZB.Top.Fault)
    for f = reshape(ZB.Top.Fault,1,[])
        top_layer_no = min(slip_model(slip_model(:,1)==f,3));
        [Wi,Di] = zero_slip_boundary(slip_model, f, top_layer_no, ZB.Top.Ratio);
        [Wt,dt] = appendRows(Wt,dt,Wi,Di);
    end
end
% BOTTOM per fault (layer index = max layer on that fault)
if ZB.Bottom.Enable && ~isempty(ZB.Bottom.Fault)
    for f = reshape(ZB.Bottom.Fault,1,[])
        bottom_layer_no = max(slip_model(slip_model(:,1)==f,3));
        [Wi,Di] = zero_slip_boundary(slip_model, f, bottom_layer_no, ZB.Bottom.Ratio);
        [Wb,db] = appendRows(Wb,db,Wi,Di);
    end
end
% LEFT / RIGHT per fault (string side)
if ZB.Left.Enable && ~isempty(ZB.Left.Fault)
    for f = reshape(ZB.Left.Fault,1,[])
        [Wi,Di] = zero_slip_boundary(slip_model, f, 'left', ZB.Left.Ratio);
        [Wl,dl] = appendRows(Wl,dl,Wi,Di);
    end
end
if ZB.Right.Enable && ~isempty(ZB.Right.Fault)
    for f = reshape(ZB.Right.Fault,1,[])
        [Wi,Di] = zero_slip_boundary(slip_model, f, 'right', ZB.Right.Ratio);
        [Wr,dr] = appendRows(Wr,dr,Wi,Di);
    end
end

%% ---------------- Build Green's & data stacks ----------------
G_stack = []; b_stack = [];       % weighted (solve)
G_raw_all = []; b_raw_all = [];   % unweighted (VR/reporting)
Gi_raw_holder = cell(1,k);

for i = 1:k
    ds = datasets(i);
    switch lower(opt.Builder)
        case 'altered'
            [Gi_raw, Gi, bi_raw, bi] = build_green_function_altered(ds.path);
        otherwise
            dt = mapKind(ds.kind);
            [Gi_raw, Gi, bi_raw, bi] = build_green_function( ...
                slip_model, ds.path, dt, 'noramp', opt.ModelType, 'pts', opt.PTs);
    end
    G_stack   = [G_stack;   ds.weight * Gi];
    b_stack   = [b_stack;   ds.weight * bi];
    G_raw_all = [G_raw_all; Gi_raw];
    b_raw_all = [b_raw_all; bi_raw];
    Gi_raw_holder{i} = Gi_raw; %#ok<AGROW>
end

%% ---------------- Add smoothing & boundaries ----------------
Greens = [G_stack;  H * (lambda / h1);  Wb; Wt; Wl; Wr];
b_full  = [b_stack; zeros(h1,1);        db;  dt;  dl;  dr];

%% ---------------- Bounds ----------------
nflt = max(slip_model(:,1));
tSm  = zeros(1, nflt+1);
fid  = slip_model(:,1);
for ii = 1:nflt, tSm(ii+1) = nnz(fid==ii); end
[lb,ub] = bounds_new(nflt, 2, tSm, 0, opt.Con);  % NS=nflt, NT=2, add_col=0

%% ---------------- Solve (lsqlin) ----------------
lsqopt = optimset('LargeScale','on','DiffMaxChange',1e-1,'DiffMinChange',1e-12, ...
    'TolCon',1e-12,'TolFun',1e-12,'TolPCG',1e-12,'TolX',1e-12, ...
    'MaxIter',1e9,'MaxPCGIter',1e9);
fn = fieldnames(opt.LSQOptions);
for ii = 1:numel(fn), lsqopt.(fn{ii}) = opt.LSQOptions.(fn{ii}); end

[u, resnorm, residual, exitflag] = lsqlin(Greens, double(b_full), [], [], [], [], lb, ub, [], lsqopt);

%% ---------------- Fill slip model; predictions ----------------
% Node-based convention: strike/dip in cols 4/5.
slip_model(:,4) = u(1:sum(tSm));
slip_model(:,5) = u(sum(tSm)+1:end);

pred = cell(1,k);
for i = 1:k
    pred{i} = Gi_raw_holder{i} * u;
end

%% ---------------- Diagnostics / stats ----------------
rms0      = sum(b_raw_all.^2);
rms       = sum((G_raw_all*u - b_raw_all).^2);
VR_pct    = 100 * (rms0 - rms) / rms0;

RMS_misfit     = sum((G_stack*u - b_stack).^2);   % data-only misfit on solve stack
rough_matrix   = H*u;
model_roughness= sqrt(sum(rough_matrix.^2) / max(1,numel(rough_matrix)));

if opt.Verbose
    fprintf('VR = %5.2f%%, rms0 = %.3e, rms = %.3e\n', VR_pct, rms0, rms);
    fprintf('resnorm = %.3e, mean(residual) = %.3e, exitflag = %d\n', sqrt(resnorm), mean(residual), exitflag);
end

stats = struct();
stats.VR_percent     = VR_pct;
stats.rms0           = rms0;
stats.rms            = rms;
stats.resnorm        = resnorm;
stats.residual_mean  = mean(residual);
stats.exitflag       = exitflag;
stats.RMS_misfit     = RMS_misfit;
stats.roughness      = model_roughness;
stats.lambda         = opt.Lambda;
stats.smoothing      = Sm;
stats.con            = opt.Con;
stats.builder        = char(opt.Builder);
stats.model_type     = char(opt.ModelType);
stats.dataset        = datasets;  % already normalized

%% ---------------- Package outputs ----------------
if nargout == 12
    % LEGACY MODE: [slip_model, RMS_misfit, model_roughness, insar_model1..8, rms]
    maxLegacy = 8;
    outPred = cell(1, maxLegacy);
    for i = 1:maxLegacy
        if i <= k
            outPred{i} = pred{i};
        else
            outPred{i} = [];
        end
    end
    varargout = [{RMS_misfit, model_roughness}, outPred, {rms}];
else
    % NEW MODE: [slip_model, stats, pred1, pred2, ...]
    varargout = [{stats}, pred];
end
end % ===== main =====

%% ================= Helpers =================
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
            S(i).kind   = char(lower(string(D{i,2})));
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
                S(i).name = sprintf('%s_%d', char(lower(string(S(i).kind))), i);
            end
            S(i).kind = char(lower(string(S(i).kind)));
        end
    else
        error('Unsupported Datasets type. Use cell or struct array.');
    end
end

function out = mergeStructs(a,b)
    out = a;
    if isempty(b), return; end
    f = fieldnames(b);
    for i = 1:numel(f)
        k = f{i};
        if isfield(out,k) && isstruct(out.(k)) && isstruct(b.(k))
            out.(k) = mergeStructs(out.(k), b.(k));
        else
            out.(k) = b.(k);
        end
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
