function [X,r,mycost,mytime] = WDC_DC(m,n,Ax,y, options)
if ~isfield(options,'maxiter'); options.maxiter = 20; end
if ~isfield(options,'verbosity'); options.verbosity = 20; end
if ~isfield(options,'ranktol'); options.ranktol = 1e-6; end
maxiter = options.maxiter;
verb = options.verbosity;
ranktol = options.ranktol;
    max_r = min(m,n);
t0 = tic;
lastcost=1;
    for r = 1:max_r
        X0 = (randn(m,n)+1i*randn(m,n))/(sqrt(2*m));
        [u0,s0,v0] = svd(X0);
        lastcostin = 1;
        for iter = 1:maxiter
            s0 = zeros(size(s0));
            s0(1:m+1:(r-1)*m+r) = 1;
            partial_X = u0*s0*v0';
            cvx_begin quiet
            variable X(m,n)
            minimize(norm_nuc(X)-real(partial_X(:)'*X(:)))
            subject to
                Ax*X(:) == y;
            cvx_end
            
            try
                [u0,s0,v0] = svd(full(X),'econ');
                if size(s0,2)<r+1
                    mycostin = 0;
                else
                    mycostin = s0(r+1,r+1);%sum(s0(r*m+r+1:m+1:end));
                end
            catch
                mycostin=1;
                fprintf('err for svd');
            end
            
            mycostin = s0(r*m+r+1);%sum(s0(r*m+r+1:m+1:end));
            err = abs(mycostin-lastcostin)/lastcostin;
            if verb>=2
                fprintf('r:%d,iter:%.3d,err:%.3e,obj:%.3e,mycost:%.3e\n',r,iter,err,cvx_optval,mycostin);
            end
            if mycostin<ranktol || err<1e-8
                break;
            end
            lastcostin = mycostin;
            X0=X;
        end
        rx = rank(X,ranktol);
        mycost = mycostin;
        rel_cost = abs(mycost-lastcost)/lastcost;
        if verb
            fprintf('rx:%d,mycost:%.3e,rel_cost:%.3e\n',r,mycost,rel_cost);
        end
        if strcmp(cvx_status,'Solved') && rx<=r%mycost<ranktol
            break;
        end
        lastcost = mycost;
    end
    [u,s,v] = svd(X);
    r = rank(X,ranktol);
    s(r+1:end,r+1:end)=0;
    X=u*s*v';
    mycost = norm(Ax*X(:) - y)/sqrt(m);
    mytime = toc(t0);
end