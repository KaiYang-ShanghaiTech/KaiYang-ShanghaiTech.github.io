function [Tset,Rset]=generate_model_file_even(K,F,cache_len)
%generate_model 此处显示有关此函数的摘要
%   此处显示详细说明
muk = cache_len*K/F; %每个file被多少个user存
fcounter = zeros(F,1);
FTset = zeros(F,K);
for i =1:K
    [val,ind] = sort(fcounter);
%     ind = find(fcounter<muk);
%     fk = randperm(length(ind),cache_len);
%     ft = ind(fk);
    ft = ind(1:cache_len);
    FTset(ft,i)=1;
    fcounter(ft) = fcounter(ft)+1;
end
    

% H1 = zeros(F,K);
% for j = 1:F/cache_len
%     H1((j-1)*cache_len+1:j*cache_len,(j-1)*muk+1:j*muk) = ones(cache_len,muk);
% end
%  n1= randperm(F);
% n2= randperm(K);
% FTset = H1(n1,n2);

% batch_num = nchoosek(K,muk);
% batch_size = F/batch_num;
% batch_ind = reshape(1:F,[batch_size,batch_num]);
% k = 1;
% FTset = false(F,K);
% for t=nchoosek(1:K,muk)'
%     msg=batch_ind(:,k);
%     FTset(msg,t)=true;
%     k = k+1;
% end

T = F*K;
Tset = false(T,K); Rset = false(T,K);
for k =1:K
    for j = 1:F
        if FTset(j,k)
            tmp1 = (j-1)*K+1:j*K;
            Tset(tmp1,k) = true;
        else
            Rset((j-1)*K+k,k) = true;
        end
    end
end
end
