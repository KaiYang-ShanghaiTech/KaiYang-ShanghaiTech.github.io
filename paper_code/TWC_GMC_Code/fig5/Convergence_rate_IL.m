% Code of Fig 5 in paper: K. Yang, Y. Shi, and Z. Ding, ¡°Generalized low-rank
% optimization for topological cooperation in ultra-dense networks,¡± IEEE
% Trans. Wireless Commun., 2019
% Written by Kai Yang
clear;addpath('./func');%addpath(genpath('./online_toolbox'));
rng('default'); rng(1);
K = 20;

optionsTR.verbosity =0;
optionsTR.costtol = 1e-8;
optionsTR.maxiter = 200;
optionsTR.tolgradnorm = 1e-10;

optionsCG = optionsTR;
optionsCG.maxiter = 15000;%1e5;
options_altmin = optionsTR;
options_altmin.maxiter = 50;
options_altmin.inmaxiter = 500;

data_stream = 3;
link_size = 0.3;
trans_cache_size = 1;
receiver_cache_size = 0;
connected_links = make_rand_topology(K,link_size);

distance = rand(K,K)*0.1+0.1;  %100m - 200m
H_channel = (randn(K,K)/sqrt(2)+1i*randn(K,K)/sqrt(2));%10^(140/20)*(10.^((-128.1-37.6*log10(distance))/20)).*

[trans_cache_set, receiver_cache_set] = make_transeiver_cache_size(K,trans_cache_size, receiver_cache_size);
 [A,b,m,n,find_index] = generate_matrix_multidatastream(connected_links,trans_cache_set,receiver_cache_set,data_stream);
r=12;
leakage_func = @(X)eval_leakage(trans_cache_set, connected_links, H_channel, find_index, X, r);
options_altmin.leakage_func = leakage_func;
optionsCG.leakage_func = leakage_func;
optionsTR.leakage_func = leakage_func;

X0.L = (randn(m,r)+1i*randn(m,r))/sqrt(2*m); X0.R = (randn(n,r)+1i*randn(n,r))/sqrt(2*n);
X1.U = X0.L; X1.V = X0.R;
 
[ X_altmin,infos_altmin,infosall_altmin ] = altmin_fixedrank_algorithm(r, [m,n],X1,options_altmin,A,b);

[X_CG, infos_CG]=Riemannian_fixedrank_CG_iter(r,  [m,n], X0, optionsCG, A,b);

[X_TR,infos_TR]=Riemannian_fixedrank_TR_iter(r,  [m,n], X0, optionsTR, A,b);

% save('convergence.mat');
co=[0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.4940    0.1840    0.5560];
set(groot,'defaultAxesColorOrder',co);
time1 = [infosall_altmin.time]; leakage1 = [infosall_altmin.leakage];
time2 = [infos_CG.time]; leakage2 = [infos_CG.leakage];
time3 = [infos_TR.time]; leakage3 = [infos_TR.leakage];

figure,clf,semilogy(1:length(leakage1),leakage1,'LineWidth',1.3);hold on
semilogy(1:length(leakage2),leakage2,'LineWidth',1.3);hold on
semilogy(1:length(leakage3),leakage3,'LineWidth',1.3);hold on
xlabel('Iteration','FontSize',14);ylabel('Interference Leakage','FontSize',14);axis tight;
xlim([0,3000]);
ylim([1e-15,inf]);
legend('AltMin','RCG','RTR');
set(gca,'FontSize',14,'FontName','Arial');

figure,clf,semilogy(time1,leakage1,'LineWidth',1.3);hold on
semilogy(time2,leakage2,'LineWidth',1.3);hold on
semilogy(time3,leakage3,'LineWidth',1.3);hold on
xlabel('Time [s]','FontSize',14);ylabel('Interference Leakage','FontSize',14);axis tight;
xlim([0,100]);
ylim([1e-15,inf]);
legend('AltMin','RCG','RTR');
set(gca,'FontSize',14,'FontName','Arial');