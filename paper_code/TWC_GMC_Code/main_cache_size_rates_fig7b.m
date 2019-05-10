% Code of Fig 7b in paper: K. Yang, Y. Shi, and Z. Ding, ¡°Generalized low-rank
% optimization for topological cooperation in ultra-dense networks,¡± IEEE
% Trans. Wireless Commun., 2019
% Written by Kai Yang
clear;addpath('./func');
addpath(genpath('./manopt'));
rng('default');rng(1);%('seed',2); randn('seed',2);
K = 20;
testnum =100;
trans_cache_size_set = [0,0.3,0.5,0.7,1];
trans_cache_size_num = length(trans_cache_size_set);
trans_power =  -20:10:60; noise_power = -120; snr_set = trans_power-noise_power;
snr_len = length(snr_set);
links = 0.2;
params.verbosity =0;
params.costtol = 1e-3;
params.maxiter = 500;
params.tolgradnorm = 1e-6;
optionsCG = params;
optionsCG.maxiter = 2e3;
options_altmin = params;
options_altmin.maxiter = 50;
options_altmin.inmaxiter = 500;
verbose = 1;
all_rank_altmin = zeros(trans_cache_size_num,snr_len,testnum);
all_rank_TR = zeros(trans_cache_size_num,snr_len,testnum);
all_rank_CG = zeros(trans_cache_size_num,snr_len,testnum);

all_time_altmin = zeros(trans_cache_size_num,snr_len,testnum);
all_time_TR = zeros(trans_cache_size_num,snr_len,testnum);
all_time_CG = zeros(trans_cache_size_num,snr_len,testnum);

all_rate_altmin = zeros(trans_cache_size_num,snr_len,testnum);
all_rate_TR = zeros(trans_cache_size_num,snr_len,testnum);
all_rate_CG = zeros(trans_cache_size_num,snr_len,testnum);


% parpool(testnum);
t=tic;
for i=1:trans_cache_size_num
    parfor testi = 1:testnum
        connected_links = make_rand_topology(K,links);
        distance = rand(K,K)*0.1+0.1;  %100m - 200m
        H_channel = (10.^((-128.1-37.6*log10(distance))/20)).*(randn(K,K)/sqrt(2)+1i*randn(K,K)/sqrt(2));
        trans_cache_size = trans_cache_size_set(i);
        trans_cache_set = make_cache_size(K,trans_cache_size);
        [A,b,m,n,find_index] = generate_matrix(connected_links,trans_cache_set);
        [ Xout_TR,rank_TR,infos_TR ] = topological_beamforming_TR( [m,n],[], params,A,b );
        for  j = 1:snr_len
            mysnr = snr_set(j);
            all_time_TR(i,j,testi) = infos_TR.time; all_rank_TR(i,j,testi)=rank_TR;
            all_rate_TR(i,j,testi) = eval_rates2(trans_cache_set,connected_links,H_channel,find_index,mysnr,Xout_TR);
            if verbose >1
                fprintf('testnum:%2.d  cache_size:%2.3d, trans power:%.3d,  rank: %2.3d,%2.3d,%2.3d,  time:%2.4d,%2.4d,%2.4d, rates:%.1f,%.1f,%.1f\n', testi, trans_cache_size, trans_power(j), all_rank_altmin(i,j,testi),all_rank_TR(i,j,testi),all_rank_CG(i,j,testi),all_time_altmin(i,j,testi),all_time_TR(i,j,testi),all_time_CG(i,j,testi),all_rate_altmin(i,j,testi),all_rate_TR(i,j,testi),all_rate_CG(i,j,testi));
            elseif verbose>0
                fprintf('testnum:%2.d  cache_size:%2.3d, trans power:%.3d,  rank: %2.3d,%2.3d,%2.3d,  rates:%.1f,%.1f,%.1f\n', testi, trans_cache_size, trans_power(j), all_rank_altmin(i,j,testi),all_rank_TR(i,j,testi),all_rank_CG(i,j,testi),all_rate_altmin(i,j,testi),all_rate_TR(i,j,testi),all_rate_CG(i,j,testi));
            end
        end
    end
end
mytime = toc(t);
% save('cache_size_all_rate.mat');

co=[0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.4940    0.1840    0.5560];
set(groot,'defaultAxesColorOrder',co);
figure,plot(trans_power(1:7),mean(all_rate_TR(1,1:7,:),3),'*-','LineWidth',1.3);hold on;
plot(trans_power(1:7),mean(all_rate_TR(4,1:7,:),3),'s-','LineWidth',1.3);hold on;
plot(trans_power(1:7),mean(all_rate_TR(5,1:7,:),3),'o-','LineWidth',1.3);hold on;
xlabel('SNR [dB]','FontSize',14);ylabel('Sum Rate [bps/Hz]','FontSize',14);axis tight;
legend('q=0 (without message sharing)','q=0.7','q=1 (full cooperation)');
%legend('TIM (q=1/K)','q=0.6','full cooperation (q=1)');%'q=0.3',
set(gca,'FontSize',14,'FontName','Arial');

