function [X,mycost,infos] = WDC_DCF_converge(m,n,Ax,y, r, X0, options)
if ~isfield(options,'maxiter'); options.maxiter = 20; end
if ~isfield(options,'verbosity'); options.verbosity = 20; end
if ~isfield(options,'ranktol'); options.ranktol = 1e-6; end
maxiter = options.maxiter;
verb = options.verbosity;
ranktol = options.ranktol;

infos(1).time = 0;
infos(maxiter).time = [];
infos(maxiter).cost = [];
t=tic;

[Q,R]=qr(Ax,0);
[u0,s0,v0] = svd(X0);
infos(1).cost = s0(r+1,r+1);
X1=X0;
for iter = 1:maxiter
    s_g = s0;
    s_g(r+1:end,r+1:end)=0;
    partial_X = u0*2*s_g*v0';

    Aps = R'*((R*R')\Q');  % right inverse AA^+ = I
    X = reshape(partial_X(:)/2-Aps*(Ax*partial_X(:))/2+Aps*y,[m,n]);
    [u0,s0,v0] = svd(full(X));
    mycost = s0(r+1,r+1);
    infos(iter+1).cost = mycost;
    infos(iter+1).time = toc(t);
%     err = sum_square(X(:)-X1(:))/sum_square(X1(:));
    if verb>=2
        fprintf('r:%d,iter:%.3d,mycost:%.3e\n',r,iter,mycost);
    end
    if mycost<ranktol || abs(mycost-infos(iter).cost)/infos(iter).cost<1e-6
        break;
    end
    X1 = X;
end
[u,s,v] = svd(full(X));
    r = rank(X,ranktol);
    s(r+1:end,r+1:end)=0;
    X=u*s*v';
    mycost = norm(Ax*X(:) - y)^2;
    infos = infos(1:iter+1);
end