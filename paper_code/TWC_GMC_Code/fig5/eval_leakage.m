function sum_leakage= eval_leakage(trans_cache_set,connected_links,H_channel,find_submat,X,r)
    K = size(trans_cache_set,1);
    P = 1;
    interference_leakage = zeros(K,1);
    [u,s,v] = svds(X,r);
    U = sqrt(s)*u';
    V = sqrt(s)*v';
    cache_set_size = cumsum([0,sum(trans_cache_set)]);

    tmp = 0;
    for i = 1:K
        stmp = norm(V(:,cache_set_size(i)+1:cache_set_size(i+1)),'fro');
        if tmp < stmp
            tmp = stmp;
        end
    end
    V = V/tmp*sqrt(P);

    X = U'*V;
    for k=1:K
        tmp2 = 0;
        for i = [1:k-1,k+1:K]
            for j=1:K
                if connected_links(k,j) &&  trans_cache_set(i,j)
                    [subrow,subcol] = find_submat(k,j,i);
                    tmp2 = tmp2+norm(H_channel(k,j)*X(subrow,subcol)*H_channel(k,j)'*X(subrow,subcol)','fro')^2;
                end
            end
            interference_leakage(k) = interference_leakage(k) + tmp2;
        end
    end
    sum_leakage = sum(interference_leakage);
end



