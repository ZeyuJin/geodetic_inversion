function pfit = fit_sinusoid_tmp(p,x,y)
% fit the LOS as the function:
% p = a*sin[(bx+cy)+phi]+d
% without topo as the input

    indx_good=~isnan(p);
    phi=p(indx_good);
    xin=x(indx_good);
    yin=y(indx_good);

    indx2=~isinf(phi);
    phi=phi(indx2);
    xin=xin(indx2);
    yin=yin(indx2);
    
    Nin=length(phi);
    G = zeros(Nin,2);
    G(:,1) = xin(:);
    G(:,2) = yin(:);
    
    % non-linear fitting
    C0 = pi / 25;
    C1 = 1;
    model_sine = @(b,G)b(1) + b(2).*sin((b(3).*G(:,1) + b(4).*G(:,2)) + b(5));
    H = mean(phi);
    amp = (max(phi) - min(phi)) / 2;
    a0 = C1 / sqrt(C1^2+1) * C0;
    b0 = -1 / sqrt(C1^2+1) * C0;
    beta0 = [H,amp,a0,b0,0];
    opts = optimset('LargeScale','on','DiffMaxChange',1e-1,'DiffMinChange',1e-12, ...
        'TolCon',1e-12,'TolFun',1e-12,'TolPCG',1e-12,'TolX',1e-12,'MaxIter',1e9,'MaxPCGIter',1e9);
    mdl = fitnlm(G,phi,model_sine,beta0,'Options',opts);
    
    % predicted values
    G_raw = [x(:),y(:)];
    pfit = predict(mdl,G_raw);
    pfit = reshape(pfit,size(x));
end