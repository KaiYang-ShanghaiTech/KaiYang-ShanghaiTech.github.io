function [trans_cache_set, receiver_cache_set] = make_transeiver_cache_size(K,trans_cache_size, receiver_cache_size)
% This function is created by Kai Yang at 21:22pm on Oct. 17, 2016
% It's used to create a full transmitter cooperation-full connected model with fixed cache size in users of K-user interference
% channel:
% omega = make_cache_size(K,cache_size)
% K: the user number
% cache_size: cache_size in each user (equal to each other)
% modified at 21:22pm on Oct. 17, 2016
    tcache = int16(round(K*trans_cache_size));
    if tcache<1
        tcache = 1;
    end
    trans_cache_set = false(K,K);
    for k = 1:K
        tmp = randperm(K-1,tcache-1);
        tmp2 = [1:k-1,k+1:K];
        trans_cache_set(tmp2(tmp),k)=true;
        trans_cache_set(k,k)=true;
    end
    rcache = int16(round(K*receiver_cache_size));
    if rcache==K
        rcache = K-1;
    end
    receiver_cache_set = false(K,K);
    for k = 1:K
        tmp = randperm(K-1,rcache);
        tmp2 = [1:k-1,k+1:K];
        receiver_cache_set(tmp2(tmp),k) = true;
    end
end