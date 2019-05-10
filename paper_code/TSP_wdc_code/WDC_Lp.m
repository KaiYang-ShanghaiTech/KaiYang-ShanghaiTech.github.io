function [X,r,mycost,mytime]=WDC_Lp( m,n,Ax,y, options)
if ~isfield(options,'maxiter'); options.maxiter = 20; end
if ~isfield(options,'inner_maxiter'); options.inner_maxiter = 500; end
if ~isfield(options,'verbosity'); options.verbosity = 2; end
if ~isfield(options,'ranktol'); options.ranktol = 1e-6; end
maxiter = options.maxiter;
inner_maxiter = options.inner_maxiter;
verb = options.verbosity;
ranktol = options.ranktol;
p = 0.5;

t0 = tic;

X0 = (randn(m,n)+1i*randn(m,n))/(sqrt(2*m));
% s = svd(X0);
gamma = 1e-2*norm(X0,'fro')^2;%1e-2;%*s(end);%1e-6*norm(X0)^2;
[Q,R]=qr(Ax,0);
Aps = R'*((R*R')\Q');
s0 = gamma^(1-p/2);
for iter = 1:maxiter
    W = (X0'*X0+gamma*eye(n))^(p/2-1);
%     if sum(isnan(W(:)))
%         fprintf('W is nan\n');
%     end
    Xtemp = X0;
%     st_cost = norm(Xtemp,'fro')^2;
    for inner_iter = 1:inner_maxiter
        X_minus_grad = Xtemp-2*s0*Xtemp*W;
%         XA = Xtemp- reshape(Aps*(Ax*Xtemp(:)-y),[m,n]);
%         XAc = reshape(Aps*(Ax*X_minus_grad(:)),[m,n]);
%         X = XA+XAc;
        X = X_minus_grad- reshape(Aps*(Ax*X_minus_grad(:)-y),[m,n]);
        inner_cost = sum_square_abs(vec(X-Xtemp));%/st_cost;
        if verb >=3
            fprintf(' iter:%d, inner iter:%d, inner cost : %.3e\n', iter, inner_iter,inner_cost);
        end
        if inner_cost<1e-6
            break;
        end
        Xtemp = X;
    end
    if verb>=2
        [u,s,v] = svd(full(X));
        r = sum(s(:) > ranktol);
        s(r+1:end,r+1:end)=0;
        Xo = u*s*v';
        mycost  = 1/2*sum_square_abs(Ax*Xo(:)-y);
        fprintf(' iter:%d, r:%d, cost:%.3e,gamma:%.3e, min singular value:%.3e\n', iter, r,mycost,gamma,s(r,r));
    end
    if norm(full(X-X0))<1e-8
        break;
    end
%     if sum(isnan(X(:)))
%         fprintf('X is nan\n');
%     end
    X0=X;
    if gamma>1e-10
        gamma = gamma/1.05;
    end
    s0 = gamma^(1-p/2);
end
% if sum(isnan(full(X)))
%     fprintf('isnan\n');
% end
[u,s,v] = svd(full(X));
r = sum(s(:) > ranktol);
if r<min(m,n)
    s(r+1:end,r+1:end)=0;
end
Xo = u*s*v';
mycost  = norm(Ax*Xo(:)-y)/sqrt(m);
mytime = toc(t0);
end