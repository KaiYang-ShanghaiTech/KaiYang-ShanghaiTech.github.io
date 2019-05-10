function [sum_rate,rates,sinr, signal_power, interference_power, noise_power]= eval_rates(trans_cache_set,connected_links,H_channel,find_index,mysnr,Xout)
    K = size(trans_cache_set,1);
    P = 10^(mysnr/10);
    signal_power = zeros(K,1);
    interference_power = zeros(K,1);
    noise_power = zeros(K,1);
    sinr = zeros(K,1);
    rates = zeros(K,1);
    if isfield(Xout,'L')
        X = Xout.L*Xout.R';
        r = size(Xout.L,2);
    else
        X=Xout;
        r = rank(X);
    end
    [u,s,v] = svds(X,r);
    U = sqrt(s)*u';
    V = sqrt(s)*v';
    cache_set_size = cumsum([0,sum(trans_cache_set)]);
%     for i = 1:K
%         tmp = norm(V(:,cache_set_size(i)+1:cache_set_size(i+1)),'fro');
%         V(:,cache_set_size(i)+1:cache_set_size(i+1)) = V(:,cache_set_size(i)+1:cache_set_size(i+1))/tmp*sqrt(P);
%     end
    tmp = 0;
    for i = 1:K
        stmp = norm(V(:,cache_set_size(i)+1:cache_set_size(i+1)),'fro');
        if tmp < stmp
            tmp = stmp;
        end
    end
    V = V/tmp*sqrt(P);
%     X = Xout.L*Xout.R';
%     r= size(Xout.L,2);
%     [u,s,v] = svds(X,r);
%     V = v'*sqrt(P);
%     U =u';
%     X = U'*V;

%     V= V*sqrt(P);
    X = U'*V;
    for k=1:K
        tmp1 = 0;
        for j=1:K
            if trans_cache_set(k,j) &&  connected_links(k,j)
                tmp1 = tmp1+H_channel(k,j)*X(find_index(k,j,k));
            end
        end
        signal_power(k) = norm(tmp1)^2;
        tmp2 = 0;
        for i = [1:k-1,k+1:K]
            for j=1:K
                if connected_links(k,j) &&  trans_cache_set(i,j)
                    tmp2 = tmp2+H_channel(k,j)*X(find_index(k,j,i));
                end
            end
            interference_power(k) = interference_power(k) + norm(tmp2)^2;
        end
        noise_power(k) = norm(U(:,k))^2;
        sinr(k) = signal_power(k)/(interference_power(k)+noise_power(k));
        rates(k) = log2(1+sinr(k))/r;
    end
    sum_rate = sum(rates);
end



