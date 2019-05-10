function [A,b,m,n] = generate_matrix_multidatastream(connected_links, trans_cache_set, receiver_cache_set, data_stream)
% This function is created by Kai Yang at 16:41 on July 20, 2017
% It's used to create a random topology model of K-user interference

    K = size(trans_cache_set,1);
    trans_cache_set_size = sum(trans_cache_set);
    m = K*data_stream;  n=sum(trans_cache_set_size)*data_stream;
    
    cache_set = cell(K,1); link_set = cell(K,1);
    for i=1:K
        cache_set{i} = find(trans_cache_set(:,i));
        link_set{i} = find(connected_links(i,:));
    end
    
    zeros_num = 0;
    for k=1:K
        for j= 1:K
            if connected_links(k,j)
                for i = 1:K
                    if trans_cache_set(i,j) && i~=k && ~receiver_cache_set(i,k)
                        zeros_num = zeros_num+data_stream^2;
%                         fprintf('zeros_num:%d\n',zeros_num);
                    end
                end
            end
        end
    end
    
    
    eqnum = K*data_stream^2+zeros_num;
    A = sparse(eqnum,m*n);
    b = [vec(repmat(eye(data_stream),[1,K]));zeros(zeros_num,1)];
    find_index = @(k,j,i,di,dj)(sum(trans_cache_set_size(1:j-1))+find(cache_set{j}==i)-1)*K*data_stream^2+(dj-1)*K*data_stream+(k-1)*data_stream+di;%   ((i-1)*K+j-1)*K+k;
    one_ind = 0;
    zeros_ind = K*data_stream^2;
    for k=1:K
        for j=1:K
            if connected_links(k,j) && trans_cache_set(k,j) 
                    con_set = zeros(data_stream,data_stream);
                    for ss = 1:data_stream
                        con_set(:,ss) = find_index(k,j,k,1,ss):find_index(k,j,k,data_stream,ss);
                    end
                    A = A | sparse(one_ind+1:one_ind+data_stream^2,con_set(:),1,eqnum,m*n);
            end
            if connected_links(k,j)
                for i = 1:K
                    if trans_cache_set(i,j) && i~=k && ~receiver_cache_set(i,k)
                        con_set = zeros(data_stream,data_stream);
                        for ss = 1:data_stream
                            con_set(:,ss) = find_index(k,j,i,1,ss):find_index(k,j,i,data_stream,ss);
                        end
                        A = A| sparse(zeros_ind+1:zeros_ind+data_stream^2,con_set(:),1,eqnum,m*n);
                        zeros_ind = zeros_ind+data_stream^2;
                    end
%                     fprintf('zeros_ind:%d\n',zeros_ind);
                end
            end
        end
        one_ind = one_ind+data_stream^2;
    end
    A = sparse(double(A));
    b = sparse(b);
end