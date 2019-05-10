function [ Xout,r,infos ] = topological_beamforming_altmin( Ksize,X0,params, A, b )
%This function is created by Kai Yang at 6:50pm on Aug. 12, 2016
tstart=tic;
if ~isfield(params,'tol')
    tol=1e-4;
else
    tol = params.tol;
end
if ~isfield(params,'costtol')
    costtol=1e-3;
else
    costtol = params.costtol;
end
if ~isfield(params,'verbosity')
    verb=1;
else
    verb = params.verbosity;
end
m=Ksize(1); n = Ksize(2);

cost_func = @(X)norm(A*X(:)-b)/sqrt(m);
fx = @(X)norm(A*X(:)-b);
r0=1;

if isempty(X0)
    X0.U = (randn(m,r0)+1i*randn(m,r0))/sqrt(2*m); X0.V = (randn(n,r0)+1i*randn(n,r0))/(sqrt(2*n));
else
    tmp.L=X0.U; tmp.R = X0.V; X0 = tmp; 
end
mycost = cost_func(X0.U*X0.V');infos.cost=[];r=r0-1;
xsol=X0;
while(mycost>costtol && r<=min(m,n))
%     myerr = mycost; iter = 0; 
    r=r+1;
    if r~=r0
        X0.U = (randn(m,r)+1i*randn(m,r))/sqrt(2*m); X0.V = (randn(n,r)+1i*randn(n,r))/(sqrt(2*n));
    end
    
    [xsol,infos_in,infos_all] = altmin_fixedrank_algorithm(r, Ksize, X0, params, A, b);
    mycost = cost_func(xsol.U*xsol.V');
    infos.cost = [infos.cost;mycost];
    if verb
        fprintf('Iteration %2.3d, mycost %7.3e.\n', r, mycost);
    end
end
Xout.L=xsol.U; Xout.R=xsol.V;
infos.time = toc(tstart);
infos.lastcost = mycost;
end

