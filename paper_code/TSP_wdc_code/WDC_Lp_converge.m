function [X,mycost,infos]=WDC_Lp_converge( m,n,Ax,y,r, X0,options)
if ~isfield(options,'maxiter'); options.maxiter = 20; end
if ~isfield(options,'inner_maxiter'); options.inner_maxiter = 500; end
if ~isfield(options,'verbosity'); options.verbosity = 2; end
if ~isfield(options,'ranktol'); options.ranktol = 1e-6; end
maxiter = options.maxiter;
verb = options.verbosity;
inner_maxiter = options.inner_maxiter;
ranktol = options.ranktol;
p = 0.5;

infos(1).time = 0;
infos(maxiter).time = [];
infos(maxiter).cost = [];
t=tic;

[u,s,v] = svd(X0);
gamma = 1e-1;%2*norm(X0,'fro')^2;
[Q,R]=qr(Ax,0);
Aps = R'*((R*R')\Q');

infos(1).cost = s(r+1,r+1);
s0 = gamma^(1-p/2);
for iter = 1:maxiter
    W = (X0'*X0+gamma*eye(n))^(p/2-1);
    Xtemp = X0;
    for inner_iter = 1:inner_maxiter
        X_minus_grad = Xtemp-2*s0*Xtemp*W;
        X = X_minus_grad- reshape(Aps*(Ax*X_minus_grad(:)-y),[m,n]);
        inner_cost = sum_square_abs(vec(X-Xtemp));%/st_cost;
        if verb >=3
            fprintf(' iter:%d, inner iter:%d, inner cost : %.3e\n', iter, inner_iter,inner_cost);
        end
%         if inner_cost<1e-6
%             break;
%         end
        Xtemp = X;
    end
    
    [u,s,v] = svd(full(X));
    mycost = s(r+1,r+1);
    infos(iter+1).cost = mycost;
    infos(iter+1).time = toc(t);
%     r = sum(s(:) > ranktol);
    if verb>=2
        fprintf(' iter:%d, r:%d, cost:%.3e,gamma:%.3e\n', iter, r,s(r+1,r+1),gamma);
    end
    if mycost<ranktol %|| abs(mycost-infos(iter).cost)/infos(iter).cost<1e-6
        break;
    end
    X0=X;
    if gamma>0%1e-12
        gamma = gamma/1.01;
    end
    s0 = gamma^(1-p/2);
end
infos = infos(1:iter+1);
end