function [xsol,infos,infos_all] = altmin_fixedrank_algorithm(r, Ksize, X0, options, A, b)
    % function exhaustive_admission_algorithm(K, r, params)
    %
    % Input:
    % K is the number of nodes/receivers/transmitters.
    % r is the rank.
    %
    % params.verbosity:  verbosity. Default is 2.
    % params.maxiter:    maximum number of T-R iterations. Default is 200.
    %
    %
    % Output:
    % xsol:             a structure with low-rank factors L and R.
    % infos:            structure containing the stats.
    %
    
    % Input
%     if r>K;error('r must be less than K!\n');end
    m=Ksize(1); n = Ksize(2);
    problem.cost = @cost;  % The input X is a structure with fields U, S, V representing a rank k matrix as U*S*V'.
    function f = cost(X)
        Xmat = X.U*X.V';
        f = .5*norm(A*Xmat(:) - b)^2;
    end

    %% Define the Euclidean gradient of the cost function nabla f(X) = P.*(X-A)
    problem.gradu = @egradu;
    function gu = egradu(X)
        % Same comment here about Xmat.
        Xmat = X.U*X.V';
        G = reshape(A'*(A*Xmat(:) - b),[m,n]);
        gu= G*X.V;
    end

    problem.gradv = @egradv;
    function gv = egradv(X)
        % Same comment here about Xmat.
        Xmat = X.U*X.V';
        G = reshape(A'*(A*Xmat(:) - b),[m,n]);
        gv = G'*X.U;
    end

    options.stopfun = @mystop;
    function y = mystop(problem, x, info, last)
        if ~mod(last,10)
            y = sqrt(2*info(last).cost/m)<options.costtol | info(last).time >800;
        else
            y =info(last).time >800;
        end
%        cost_func = @(X)norm(A*reshape(X.L*X.R',[m*n,1])-b)/sqrt(m); 
    end

    [ xsol,~,infos,infos_all ] = altmin_algo( X0,problem,options);
    %[xsol, ~, infos, ~] = trustregions(problem, X0, options);
    Xsol = xsol.U*xsol.V';
end
