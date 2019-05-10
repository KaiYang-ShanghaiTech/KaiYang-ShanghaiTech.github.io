function cache_set = make_cache_size(K,cache_size)
% This function is created by Kai Yang at 21:22pm on Oct. 17, 2016
% It's used to create a full transmitter cooperation-full connected model with fixed cache size in users of K-user interference
% channel:
% omega = make_cache_size(K,cache_size)
% K: the user number
% cache_size: cache_size in each user (equal to each other)
% modified at 21:22pm on Oct. 17, 2016
    cache_set = rand(K,K)<cache_size;
    cache_set(1:K+1:end)=true;
%     a = false(1,K);
%     for j=1:K
%         tmp = [1:j-1,j+1:K];
%         xi = randperm(K-1,cache_size-1);
%         b = a;
%         b(tmp(xi))=true;
%         b(j)=true;
%         cache_set(:,j)=b;
%     end
end