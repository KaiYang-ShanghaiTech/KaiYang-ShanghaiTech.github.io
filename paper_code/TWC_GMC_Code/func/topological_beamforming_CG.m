function [ Xout,r,infos ] = topological_beamforming_CG( Ksize,X0, params,A,b )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
tstart = tic;
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
if ~isfield(params,'tolgradnorm')
    params.tolgradnorm=1e-6;
end
if ~isfield(params,'maxiter')
    params.maxiter=500;
end
m=Ksize(1);  n=Ksize(2);
r0 = 1;
if isempty(X0)
    X0.L = (randn(m,r0)+1i*randn(m,r0))/sqrt(2*m); X0.R = (randn(n,r0)+1i*randn(n,r0))/sqrt(2*n);
end
cost_func = @(X)norm(A*reshape(X.L*X.R',[m*n,1])-b)/sqrt(m);

infos.cost = [];

for r = 1:min(m,n)
        X0.L = (randn(m,r)+1i*randn(m,r))/sqrt(2*m); X0.R = (randn(n,r)+1i*randn(n,r))/sqrt(2*n);
%         X0.L = randn(m,r)/sqrt(m); X0.R = randn(n,r)/sqrt(n);
        [Xcg, infos_iter]=Riemannian_fixedrank_CG(r, Ksize, X0, params, A, b);
        mycost=cost_func(Xcg);
        infos.cost = [infos.cost;mycost];
        if mycost<costtol
            break;
        end
        if verb
            fprintf('r: %3d, mycost: %7.3e \n', r, mycost);
        end
%         X0.L = (randn(K,r+1)+1i*randn(K,r+1))/sqrt(2*K); X0.R = (randn(K^2,r+1)+1i*randn(K^2,r+1))/(K*sqrt(2));
%% Rank-r update: Riemannian
end
Xout = Xcg;
infos.time = toc(tstart);
infos.lastcost = cost_func(Xout);
end