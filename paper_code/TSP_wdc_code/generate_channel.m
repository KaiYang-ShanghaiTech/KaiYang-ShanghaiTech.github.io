function H = generate_channel(distance,M,N,L1,L2)
    H = zeros(M,N,L1,L2);
    for i = 1:M
        for j = 1:N
            H(i,j,:,:) = 10^((-128.1-37.6*log10(distance(i,j)))/20)*(randn(L1,L2)/sqrt(2)+1i*randn(L1,L2)/sqrt(2));
        end
    end
end