function [Xout, infos]=Riemannian_fixedrank_TR(r, Ksize, X0, options, A, b)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if ~isfield(options,'tolgradnorm'); options.gradtol=1e-6; end
if ~isfield(options,'verbosity'); options.verbosity=1; end
m=Ksize(1);  n=Ksize(2);
% Pick the manifold of fixed-rank matrices
    problem.M = symfixedrankYYcomplexfactory(m+n, r);

    %% Define the problem cost function f(X) = 1/2 * || P.*(X-A) ||^2
    problem.cost = @cost;  % The input X is a structure with fields U, S, V representing a rank k matrix as U*S*V'.
    function [f, store] = cost(X, store)
        if ~isfield(store,'Xmat')
            store.Xmat = X(1:m,:)*X(m+1:end,:)';
        end
        f = .5*norm(A*store.Xmat(:) - b)^2;
    end

    %% Define the Euclidean gradient of the cost function nabla f(X) = P.*(X-A)
    problem.egrad = @egrad;
    function [g,store] = egrad(X,store)
        % Same comment here about Xmat.
        if ~isfield(store,'Xmat')
            store.Xmat = X(1:m,:)*X(m+1:end,:)';
        end
        XL = X(1:m,:); XR = X(m+1:end,:);
        if ~isfield(store,'G')
            store.G = reshape(A'*(A*store.Xmat(:) - b),[m,n]);
        end
        g.L= store.G*XR;
        g.R = store.G'*XL;
        g = [g.L;g.R];
    end


    %% Define the Euclidean Hess  
    problem.ehess = @ehess;
    function [Hess,store] = ehess(X, eta, store)
        if ~isfield(store,'Xmat')
            store.Xmat = X(1:m,:)*X(m+1:end,:)';
        end
        XL = X(1:m,:); XR = X(m+1:end,:);
        if ~isfield(store,'G')
            store.G = reshape(A'*(A*store.Xmat(:) - b),[m,n]);
        end
        etaL = eta(1:m,:); etaR = eta(m+1:end,:);
        
        tmp = etaL*XR' + XL*etaR';
        Pdot  = reshape(A'*A*tmp(:),[m,n]);
        
        Hess.L = store.G * etaR + Pdot*XR;
        Hess.R = store.G' * etaL + Pdot'*XL;
        
        Hess = [Hess.L;Hess.R];
    end

%      figure,checkgradient(problem);
%      figure,checkhessian(problem);
    %% If mycost exists, we use our stopping criteron, else default
options.stopfun = @mystop;
    function y = mystop(problem, x, info, last)
       y = sqrt(2*info(last).cost/m)<options.costtol | info(end).time >800;
    end

leakage_func = options.leakage_func;
options.statsfun = @mystatsfun;
function stats = mystatsfun(problem, X, stats)
       thisleakage = leakage_func(X(1:m,:)*X(m+1:end,:)');
    stats.leakage = thisleakage;
end
    %% 
       Xs = [X0.L;X0.R];
       [Xcg, xcost, infos, options] = trustregions(problem, Xs, options);
       Xout.L = Xcg(1:m,:);
       Xout.R = Xcg(m+1:end,:);
end

