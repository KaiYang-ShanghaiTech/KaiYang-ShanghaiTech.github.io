function [X,r,mycost,mytime] = WDC_DCF(m,n,Ax,y, options)
if ~isfield(options,'maxiter'); options.maxiter = 20; end
if ~isfield(options,'verbosity'); options.verbosity = 20; end
if ~isfield(options,'ranktol'); options.ranktol = 1e-6; end
maxiter = options.maxiter;
verb = options.verbosity;
ranktol = options.ranktol;
t0 = tic;
[Q,R]=qr(Ax,0);
Aps = R'*((R*R')\Q');
    max_r = min(m,n);
    X0 = randn(m,n);
    lastcost = 1;
    for r = 1:max_r
        X0 = (randn(m,n)+1i*randn(m,n))/(sqrt(2*m));%(randn(m,r)+1i*randn(m,r))/(sqrt(2*m))* (randn(r,n)+1i*randn(r,n))/(sqrt(2*n));%(randn(m,n)+1i*randn(m,n))/sqrt(2*m);
        [u0,s0,v0] = svd(X0);
        X1=X0;
        lastcostin = 1;
        for iter = 1:maxiter
            s_g = s0;
            s_g(r+1:end,r+1:end)=0;
            partial_X = u0*2*s_g*v0';
            X = reshape(partial_X(:)/2-Aps*(Ax*partial_X(:))/2+Aps*y,[m,n]);
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
            err = abs(mycostin-lastcostin)/lastcostin;%sum_square(X(:)-X1(:))/sum_square(X1(:));
            if verb>=2
                fprintf('r:%d,iter:%.3d,err:%.3e,mycost:%.3e\n',r,iter,err,mycostin);
            end
            if mycostin<ranktol || err<1e-8
                break;
            end
            X1 = X;
            lastcostin = mycostin;
        end
        rx = rank(X,ranktol);
        mycost = mycostin;
        rel_cost = abs(mycost-lastcost)/lastcost;
        if verb
            fprintf('rx:%d,mycost:%.3e,rel_cost:%.3e\n',r,mycost,rel_cost);
        end
        if rx<=r%mycost<ranktol
            break;
        end
%         if rel_cost<1e-6
%             X0 = X+1e-4*randn(m,n);
%         else
%             X0 = X;
%         end
        lastcost = mycost;
        %X0 = X;
    end
    [u,s,v] = svd(X);
    r = rank(X,ranktol);
    s(r+1:end,r+1:end)=0;
    X=u*s*v';
    mycost = norm(Ax*X(:) - y)/sqrt(m);
    mytime = toc(t0);
end