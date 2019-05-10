function [X,mycost,infos] = WDC_DC_converge(m,n,Ax,y, r,X0,options)
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
        [u0,s0,v0] = svd(full(X0));
        infos(1).cost = s0(r+1,r+1);
        
        obj = 0;
        s0 = zeros(m,n);
        s0(1:m+1:(r-1)*m+r) = 1;
        partial_X = u0*s0*v0';
        for iter = 1:maxiter
            cvx_begin quiet
            variable X(m,n)
            minimize(norm_nuc(X)-real(partial_X(:)'*X(:)))
            subject to
                Ax*X(:) == y;
            cvx_end
            [u0,s0,v0] = svd(full(X));
            mycost = s0(r+1,r+1);
            infos(iter+1).cost = mycost;
            infos(iter+1).time = toc(t);
            err = abs(obj - cvx_optval);
            if verb>=2
                fprintf('r:%d,iter:%.3d,err:%.3e,obj:%.3e,mycost:%.3e\n',r,iter,err,cvx_optval,mycost);
            end
            if mycost<ranktol
                break;
            end
            obj = cvx_optval;
            s0 = zeros(m,n);
            s0(1:m+1:(r-1)*m+r) = 1;
            partial_X = u0*s0*v0';
        end
    [u,s,v] = svd(full(X));
    r = rank(X,ranktol);
    s(r+1:end,r+1:end)=0;
    X=u*s*v';
    mycost = norm(Ax*X(:) - y)^2;
    infos = infos(1:iter+1);
end