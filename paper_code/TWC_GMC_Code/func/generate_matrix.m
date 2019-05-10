function [A,b,m,n,find_index] = generate_matrix(connected_links,trans_cache_set,flag)
% This function is created by Kai Yang at 6:19pm on Aug. 12, 2016
% It's used to create a random topology model of K-user interference
% channel:
% omega = make_rand_topology(K,interference_links)
% K: the user number
% interference_links: number of links
% modified at 10:30 am, Aug. 13,2016
if nargin<3;flag=1;end

if flag %full cooperation
    K = size(trans_cache_set,1);
    cache_set_size = sum(trans_cache_set);
    m=K;n = sum(cache_set_size);
    cache_set = cell(K,1); link_set = cell(K,1);
    for i=1:K
        cache_set{i} = find(trans_cache_set(:,i));
        link_set{i} = find(connected_links(i,:));
    end
%     cache_set(1:K+1:K^2)=false;
%     cache_link_set = ~cache_link_set;
    zeros_num = 0;
    for k=1:K
        for i=link_set{k}
            tmp = trans_cache_set(:,i); 
%             tmp(k)=0;
            zeros_num = zeros_num+sum(tmp)-tmp(k);
        end
    end
%     zeros_num = connected_links.*sum(cache_set)
%     sum(sum(connected_links,2) .* sum(cache_set,2));
    A = zeros(K+zeros_num,K*sum(cache_set_size));
    b = [ones(K,1);zeros(zeros_num,1)];
    find_index = @(k,i,j)(sum(cache_set_size(1:i-1))+find(cache_set{i}==j)-1)*K+k;%   ((i-1)*K+j-1)*K+k;
    zeros_ind = K;
%     con_set = cell(K,1); tmp2=cell(K,1);
    for k=1:K
        tmp = find(connected_links(k,:));
        for i=tmp
            if sum(cache_set{i}==k)
                con_set = find_index(k,i,k);
                A(k,con_set)=1;
            end
        end
        
        for i=tmp
            tmp2 = find(trans_cache_set(:,i));
%             tmp3 = setdiff(tmp2,k);
            for j = tmp2'
                if j~=k
                    zeros_ind = zeros_ind+1;
                    A(zeros_ind,find_index(k,i,j))=1;
                end
            end
        end
    end
    A = sparse(A);
    b = sparse(b);
end