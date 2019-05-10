% Code of Fig 7a in paper: K. Yang, Y. Shi, and Z. Ding, ¡°Generalized low-rank
% optimization for topological cooperation in ultra-dense networks,¡± IEEE
% Trans. Wireless Commun., 2019
% Written by Kai Yang
clear;addpath('./func');
addpath(genpath('./manopt'));
rng('default');rng(1);%('seed',2); randn('seed',2);
K = 20;
testnum =500;
trans_cache_size_set = [0:0.2:1];
trans_cache_size_num = length(trans_cache_size_set);
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

all_rank_altmin = zeros(trans_cache_size_num,testnum);
all_rank_TR = zeros(trans_cache_size_num,testnum);
all_rank_CG = zeros(trans_cache_size_num,testnum);

all_time_altmin = zeros(trans_cache_size_num,testnum);
all_time_TR = zeros(trans_cache_size_num,testnum);
all_time_CG = zeros(trans_cache_size_num,testnum);

all_rate_altmin = zeros(trans_cache_size_num,testnum);
all_rate_TR = zeros(trans_cache_size_num,testnum);
all_rate_CG = zeros(trans_cache_size_num,testnum);

trans_power = -80; noise_power = -102; mysnr = trans_power-noise_power;
% parpool(testnum);
t=tic;
for  j = 1: trans_cache_size_num
    trans_cache_size = trans_cache_size_set(j);
    for i=1:testnum
        trans_cache_set = make_cache_size(K,trans_cache_size);
        connected_links = make_rand_topology(K,links);
        H_channel = randn(K,K);
        [A,b,m,n,find_index] = generate_matrix(connected_links,trans_cache_set);
            [Xout_altmin,rank_altmin,infos_altmin] = topological_beamforming_altmin([m,n],[],options_altmin,A,b);
            all_time_altmin(j,i) = infos_altmin.time; all_rank_altmin(j,i)=rank_altmin;
            all_rate_altmin(j,i) = eval_rates(trans_cache_set,connected_links,H_channel,find_index,mysnr,Xout_altmin);
            
            [ Xout_CG,rank_CG,infos_CG ] = topological_beamforming_CG( [m,n],[], optionsCG,A,b );
            all_time_CG(j,i) = infos_CG.time; all_rank_CG(j,i)=rank_CG;
            all_rate_CG(j,i) = eval_rates(trans_cache_set,connected_links,H_channel,find_index,mysnr,Xout_CG);

            [ Xout_TR,rank_TR,infos_TR ] = topological_beamforming_TR( [m,n],[], params,A,b );
            all_time_TR(j,i) = infos_TR.time; all_rank_TR(j,i)=rank_TR;
            all_rate_TR(j,i) = eval_rates(trans_cache_set,connected_links,H_channel,find_index,mysnr,Xout_TR);

            fprintf('testnum:%2.d  cache_size:%2.3d  rank: %2.3d,%2.3d,%2.3d,  time:%2.4d,%2.4d,%2.4d\n', i, trans_cache_size,all_rank_altmin(j,i),all_rank_TR(j,i),all_rank_CG(j,i),all_time_altmin(j,i),all_time_TR(j,i),all_time_CG(j,i));
    end
end
mytime = toc(t);
save('cache_size_all_500.mat');
% exit
altmin_result = mean(all_rank_altmin,2);
TR_result = mean(all_rank_TR,2);
CG_result = mean(all_rank_CG,2);

altmin_time = sum(mean(all_time_altmin,2));
TR_time = sum(mean(all_time_TR,2));
CG_time = sum(mean(all_time_CG,2));

co=[0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.4940    0.1840    0.5560];
set(groot,'defaultAxesColorOrder',co);
figure,plot(trans_cache_size_set,1./altmin_result,'*-','LineWidth',1.3);hold on
plot(trans_cache_size_set,1./CG_result,'s-','LineWidth',1.3);
plot(trans_cache_size_set,1./TR_result,'o-','LineWidth',1.3);
xlabel('Cooperation Level','FontSize',14);ylabel('DoF','FontSize',14);axis tight;
legend('AltMin','RCG','RTR');
% ylim([0.2,0.35]);
set(gca,'FontSize',14,'FontName','Arial');