function [Ax,y,row_num,col_num]=generate_parameters(Tset,Rset,L,d,K,T,H)
%generate_model 此处显示有关此函数的摘要
%   此处显示详细说明
row_num = sum(Rset(:))*L*d;
col_num = sum(Tset(:))*L*d;
cum_Rset = cumsum([0,sum(Rset)]);
cum_Tset = cumsum([0,sum(Tset)]);

find_index = @myfind_index;
    function ind = myfind_index(k,l,i,j,m,n)
        rows = ((cum_Rset(k)+sum(Rset(1:l-1,k)))*L+m-1)*d+1:((cum_Rset(k)+sum(Rset(1:l-1,k)))*L+m)*d; % (((k-1)*T+l-1)*L+m-1)*d+1:(((k-1)*T+l-1)*L+m)*d;
        cols =((cum_Tset(i)+sum(Tset(1:j-1,i)))*L+n-1)*d+1:((cum_Tset(i)+sum(Tset(1:j-1,i)))*L+n)*d; % (((i-1)*T+l-1)*L+n-1)*d+1:(((i-1)*T+l-1)*L+n)*d;
        ind = zeros(d,d);
        for ss = 1:d
            ind(:,ss) = (cols(ss)-1)*row_num+rows(1):(cols(ss)-1)*row_num+rows(d);
        end
    end

one_len = sum(Rset(:))*d^2;
zero_len =sum((T-sum(Tset)-1).*sum(Rset))*d^2;
Ax = sparse(one_len+zero_len,row_num*col_num);
y = sparse(one_len+zero_len,1);
% y = sparse(1:one_len,1,1,one_len+zero_len,1);
ind = 0;
for k = 1:K
    for l = 1:T
        if Rset(l,k)
            for i = 1:K
                if Tset(l,i)
                    for m = 1:L
                        for n = 1:L
                            tmp = find_index(k,l,i,l,m,n);
                            Ax = Ax+sparse(ind+1:ind+d^2,tmp(:),H(k,i,m,n),one_len+zero_len,row_num*col_num);
                        end
                    end
                end
            end
            y(ind+1:ind+d^2)=vec(eye(d));
            ind = ind+d^2;
        end
    end
end
 
for k = 1:K
    for l = 1:T
        if Rset(l,k)
            for j = 1:T
                if ~Tset(j,k)&& j~=l
                    for i = 1:K
                        if Tset(j,i)
                            for m = 1:L
                                for n = 1:L
                                    tmp = find_index(k,l,i,j,m,n);
                                    Ax = Ax+sparse(ind+1:ind+d^2,tmp(:),H(k,i,m,n),one_len+zero_len,row_num*col_num);
                                end
                            end
                        end
                    end
                    ind = ind+d^2;
                end
            end
        end
    end
end

end
