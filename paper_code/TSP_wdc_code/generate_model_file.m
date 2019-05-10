function [Tset,Rset]=generate_model_file(K,F,cache_len)
%generate_model 此处显示有关此函数的摘要
%   此处显示详细说明

while(1)
    FTset = false(F,K);
    for k =1:K
        tmp = randperm(F,F);
        FTset(tmp(1:cache_len),k) = true;
    end
    if isempty(find(~sum(double(FTset),2),1))
        break;
    end
end
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
