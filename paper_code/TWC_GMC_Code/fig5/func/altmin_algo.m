function [ Xout,thiscost,infos,infos_all ] = altmin_algo( X0,problem,params)
%This function is created by Kai Yang at 6:50pm on Aug. 12, 2016
tstart=tic;
if ~isfield(params,'tol')
    intol=1e-6;
else
    intol = params.tol;
end
if ~isfield(params,'tol')
    tol=1e-6;
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
if ~isfield(params,'inmaxiter')
    inmaxiter = 1000;
else
    inmaxiter = params.inmaxiter;
end
if ~isfield(params,'maxiter')
    maxiter = 30;
else
    maxiter = params.maxiter;
end
if ~isfield(params,'stepsize')
    stepsize = 10;
else
    stepsize = params.stepsize;
end
leakage_func = params.leakage_func;
c_constant = 0.5;
tau = 0.8;

% rate_func = params.rate_func;
cost_func = problem.cost;%@(X)0.5*norm(Af(X)-b)^2;
grad_u = problem.gradu;
grad_v = problem.gradv;

% r0=1; 
% m = size(X0.U,1);
% n = size(X0.V,1);
r = size(X0.U,2);
infos(maxiter).iter=maxiter-1;
infos(maxiter).cost=[];%r=r0-1;
infos(maxiter).time=[];

    riter = 0; cost0 = cost_func(X0); X = X0; infos(1).iter = riter; infos(1).cost = cost0; infos(1).time = toc(tstart); 
    thisleakage = leakage_func(X0.U*X0.V');
    infos(1).leakage =thisleakage;
    infos_all = infos; overall_iter = 0; stopflag=0; 
    while(riter<maxiter)
        riter = riter+1;
        %% V update
%         g1 = @(X)grad_v(X);
        inneriter=0; objcost = cost_func(X0); lastcost = objcost;
        while(1)
            overall_iter = overall_iter + 1; % for the overall infos
            inneriter = inneriter+1;
            alpha = stepsize;
            gv = grad_v(X0);
            t = c_constant*norm(gv,'fro')^2;
            j = 0;
            while(1)
               j = j + 1;
               X.V = X0.V-alpha*gv;
               objcost = cost_func(X);
               if objcost <= lastcost-alpha*c_constant*norm(gv,'fro')^2
                   break;
               else
                   alpha = tau*alpha;
               end
            end
            thiscost = objcost;
            myerr = abs(thiscost-lastcost)/abs(lastcost);
            infos_all(overall_iter+1).iter = overall_iter; 
            infos_all(overall_iter+1).cost = thiscost;
            infos_all(overall_iter+1).time = toc(tstart);
            thisleakage = leakage_func(X0.U*X.V');
            infos_all(overall_iter+1).leakage = thisleakage;
            if verb>=3
                fprintf('r: %d, Iteration: %.2d,%.3d, mycost %.5e, myerror %7.3e, leakage %.3e.\n', r, riter, inneriter, thiscost,myerr, thisleakage);
            end
%             if thiscost<costtol; break; end
            if (inneriter>inmaxiter || myerr<intol)
                break;
            end
            
            if isfield(params,'stopfun') && params.stopfun(problem,[],infos_all,overall_iter+1)
                if verb>=2
                    fprintf('user defined stop function satisfied.\n');
                end
                stopflag = 1;
                break;
            end
                
            X0.V = X.V;
            lastcost = thiscost;
        end
       % U update
%         g2 = @(U)grad_u(U,V);
        if stopflag
            infos(riter+1).iter = riter; 
            infos(riter+1).cost = thiscost;
            infos(riter+1).time = toc(tstart);
            thisleakage = leakage_func(X.U*X.V');
            infos(riter+1).leakage = thisleakage;
            break;
        end
        X0 = X;
        inneriter=0; objcost = cost_func(X0); lastcost = objcost;
        while(1)
            overall_iter = overall_iter + 1; % for the overall infos
            inneriter = inneriter+1;
            alpha = stepsize;
            gu = grad_u(X0);
            while(1)
               X.U = X0.U-alpha*gu;
               objcost = cost_func(X);
               if objcost <= lastcost-alpha*c_constant*norm(gu,'fro')^2;
                   break;
               else
                   alpha = tau*alpha;
               end
%                fprintf('alpha=%.3e\n',alpha);
            end
            thiscost = objcost;
            myerr = abs(thiscost-lastcost)/abs(lastcost);
            infos_all(overall_iter+1).iter = overall_iter; 
            infos_all(overall_iter+1).cost = thiscost;
            infos_all(overall_iter+1).time = toc(tstart);
            thisleakage = leakage_func(X0.U*X.V');
            infos_all(overall_iter+1).leakage = thisleakage;
            if verb>=3
                fprintf('r: %d, Iteration: %.2d,%.3d, mycost %.5e, myerror %7.3e,  leakage %.3e.\n', r, riter, inneriter, thiscost,myerr,thisleakage);
            end
%             if thiscost<costtol; break; end
            if (inneriter>inmaxiter || myerr<intol)
                break;
            end
            
            if isfield(params,'stopfun') && params.stopfun(problem,[],infos_all,overall_iter+1)
                if verb>=2
                    fprintf('user defined stop function satisfied.\n');
                end
                stopflag = 1;
                break;
            end
            
            X0.U = X.U;
            lastcost = thiscost;
        end
        if stopflag
            infos(riter+1).iter = riter; 
            infos(riter+1).cost = thiscost;
            infos(riter+1).time = toc(tstart);
            thisleakage = leakage_func(X.U*X.V');
            infos(riter+1).leakage = thisleakage;
            break;
        end
        infos(riter+1).iter = riter; 
        infos(riter+1).cost = thiscost;
        infos(riter+1).time = toc(tstart);
%         thissumrate = rate_func(X.U*X.V');
        thisleakage = leakage_func(X.U*X.V');
        infos(riter+1).leakage = thisleakage;
        outerr = abs(thiscost-cost0)/abs(cost0);
        if verb>=2
            fprintf('Iteration: %.2d, mycost %.3e, outerror %7.3e, leakage %.3e.\n', riter, thiscost,outerr,thisleakage);
        end
        if outerr<tol ; break; end
        cost0 = thiscost;
    end
Xout = X;
infos = infos(1:riter+1);
infos_all = infos_all(1:overall_iter+1);
end

