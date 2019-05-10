function [X,r,mycost,mytime]=WDC_nuc( m,n,Ax,y, options)
if ~isfield(options,'verbosity'); options.verbosity = 2; end
if ~isfield(options,'ranktol'); options.ranktol = 1e-6; end
verb = options.verbosity;
ranktol = options.ranktol;

t0=tic;

cvx_begin quiet
cvx_solver sdpt3
variable X(m,n)
minimize(norm_nuc(X))
subject to
    Ax*X(:) == y;
cvx_end

[u,s,v] = svd(full(X));
r = sum(s(:) > ranktol);
s(r+1:end,r+1:end)=0;
Xo = u*s*v';
mycost  = 1/2*sum_square_abs(Ax*Xo(:)-y);
if verb>=2
    fprintf(' r:%d, cost:%.3e\n', r,mycost);
end
mytime = toc(t0);
end