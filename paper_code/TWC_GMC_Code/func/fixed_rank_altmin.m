function [ Xout,r,infos ] = fixed_rank_altmin( r,K,X0,params, A, b )
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
cost_func = @(X)norm(A*X(:)-b)/sqrt(K);
fx = @(X)norm(A*X(:)-b);

m=K;
n = size(A,2)/K;
if isempty(X0)
    U0 = (randn(m,r)+1i*randn(m,r))/sqrt(2*m); V0 = (randn(n,r)+1i*randn(n,r))/(sqrt(2*n));
else
    U0=X0.U; V0=X0.V;
end
mycost = cost_func(U0*V0');infos.cost=[];r=r0-1;

    myerr = mycost; iter = 0; 
    if r~=r0
        U0 = (randn(m,r)+1i*randn(m,r))/sqrt(2*m); V0 = (randn(n,r)+1i*randn(n,r))/(sqrt(2*n));
%         scales = norm(U0*V0','fro');
%         U0=U0/sqrt(scales);
%         V0=V0/sqrt(scales);
    end
    while(~(myerr<tol && iter>0)&&iter<20)
        %% initialize
        iter = iter+1;
        %% V update
        cvx_begin quiet
        variable V(n,r) complex
        minimize fx(U0*V')
        cvx_end
        
        %% U update
        cvx_begin quiet
        variable U(m,r) complex
        minimize fx(U*V')
        cvx_end
        
        myerr = abs(cost_func(U0*V0')-cost_func(U*V'));
        if verb==2
            fprintf('     r: %2.3d, Iteration: %3d, myerror %7.3e.\n', r, iter, myerr);
        end
        U0=U;V0=V;
    end
    mycost = cost_func(U*V');
    infos.cost = [infos.cost;mycost];
    if verb
        fprintf('Iteration %2.3d, mycost %7.3e.\n', r, mycost);
    end

Xout.L=U; Xout.R=V;
infos.time = toc(tstart);
end

