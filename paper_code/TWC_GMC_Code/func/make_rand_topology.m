function connected_links = make_rand_topology(K,p)
% This function is created by Kai Yang at 6:19pm on Aug. 12, 2016
% It's used to create a random topology model of K-user interference
% channel:
% connected_links = make_rand_topology(K,links)
% K: the user number
% links: number of links
% modified at 10:30 am, Aug. 13,2016
    if p<0||p>1
        error('Links must be [0,1]!');
    end
    omega = rand(K,K)<p;
    omega(1:K+1:end)=1;
%     omega = eye(K);  % at least one different link is connected for each user
%     interference_links = links-K;
%     tmp=randperm(K^2-K,interference_links);
%     obs_set = setdiff(1:K^2,1:K+1:K^2);
%     omega(obs_set(tmp))=1;
    connected_links = logical(omega);
end