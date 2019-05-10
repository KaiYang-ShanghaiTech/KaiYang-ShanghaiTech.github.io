% Code of Fig 6 in paper: K. Yang, Y. Shi, and Z. Ding, ¡°Generalized low-rank
% optimization for topological cooperation in ultra-dense networks,¡± IEEE
% Trans. Wireless Commun., 2019
% Written by Kai Yang
clear;addpath('./func');%addpath(genpath('./online_toolbox'));
rand('seed',1); randn('seed',1);
K = 20;
testnum = 100;
links_size_set = [0.1:0.1:0.6];
links_size_len = length(links_size_set);
trans_cache_size = 1;

params.verbosity =0;
params.costtol = 1e-3;
params.maxiter = 500;
params.tolgradnorm = 1e-6;

optionsCG = params;
optionsCG.maxiter = 2e3;
options_altmin = params;
options_altmin.maxiter = 50;
options_altmin.inmaxiter = 500;

all_rank_altmin = zeros(links_size_len,testnum);
all_rank_TR = zeros(links_size_len,testnum);
all_rank_CG = zeros(links_size_len,testnum);

all_time_altmin = zeros(links_size_len,testnum);
all_time_TR = zeros(links_size_len,testnum);
all_time_CG = zeros(links_size_len,testnum);

t=tic;
trans_cache_set = make_cache_size(K,trans_cache_size);
for  j = 1: links_size_len
    link_size = links_size_set(j);
    parfor i=1:testnum
        connected_links = make_rand_topology(K,link_size);
        [A,b] = generate_matrix(connected_links,trans_cache_set);
        
        [X_altmin,rank_altmin,infos_altmin] = topological_beamforming_altmin(K,[],options_altmin,A,b);
        all_time_altmin(j,i) = infos_altmin.time; all_rank_altmin(j,i)=rank_altmin;
 
        [ X_CG,rank_CG,infos_CG ] = topological_beamforming_CG( K,[], optionsCG,A,b );
        all_time_CG(j,i) = infos_CG.time; all_rank_CG(j,i)=rank_CG;
        
        [ X_TR,rank_TR,infos_TR ] = topological_beamforming_TR( K,[], params,A,b );
        all_time_TR(j,i) = infos_TR.time; all_rank_TR(j,i)=rank_TR;
        
        fprintf('testnum:%2.d  links_size:%2.3d  rank: %2.3d,%2.3d,%2.3d,  time:%2.4d,%2.4d,%2.4d\n', i, link_size,all_rank_altmin(j,i),all_rank_TR(j,i),all_rank_CG(j,i),all_time_altmin(j,i),all_time_TR(j,i),all_time_CG(j,i));
    end
end
mytime = toc(t)
altmin_result = mean(all_rank_altmin,2);
TR_result = mean(all_rank_TR,2);
CG_result = mean(all_rank_CG,2);

save('topology.mat');

co=[0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.4940    0.1840    0.5560];
set(groot,'defaultAxesColorOrder',co);
figure,plot(links_size_set(2:end),1./altmin_result(2:end),'*-','LineWidth',1.3);hold on
plot(links_size_set(2:end),1./CG_result(2:end),'s-','LineWidth',1.3);
plot(links_size_set(2:end),1./TR_result(2:end),'o-','LineWidth',1.3);
xlabel('Number of Links','FontSize',14);ylabel('Achievable DoF','FontSize',14);axis tight;
legend('AltMin','RCG','RTR');
set(gca,'FontSize',14,'FontName','Arial');
% savefig('topology.fig');