% Code for Figure 3 of paper:  K. Yang, Y. Shi, and Z. Ding, ¡°Data shuffling in wireless distributed
% computing via low-rank optimization,¡± IEEE Trans. Signal Process.,2019.
% Written by Kai Yang
clear all;
rng('default');
rng(1);
K=5;
F = 10;
T=K*F;
d = 1;
cache_len= 5;  
L=1;
M=1;
verb =2;
options.maxiter =300;
options.verbosity = verb;
options.costtol = 1e-6;
options.ranktol = 1e-6;
options2 = options;
options2.maxiter =5000;
options2.inner_maxiter = 500;

H = randn(K,K,L,L);
[Tset,Rset]=generate_model_file(K,F,cache_len);
[Ax,y,m,n] = generate_parameters(Tset,Rset,L,d,K,T,H);
r =13;
X0 =  (randn(m,n)+1i*randn(m,n))/(sqrt(2*m));

[X2,mycost2,infos2] = WDC_DC_converge(m,n,Ax,y, r,X0,options);
[X3,mycost3,infos3]=WDC_DCF_converge( m,n,Ax,y, r,X0,options2);
[X1,mycost1,infos1] = WDC_Lp_converge(m,n,Ax,y, r,X0,options2);
cost1 = [infos1(1:end).cost];   
cost2 = [infos2(1:end).cost];
cost3 = [infos3(1:end).cost];
time1 = [infos1(1:end).time];
time2 = [infos2(1:end).time];
time3 = [infos3(1:end).time];

allcost1 = cost1;
allcost2 = cost2;
allcost3 = cost3;
iter1 = length(allcost1);
iter2 = length(allcost2);
iter3 = length(allcost3);
save('convergence_K=10.mat');
co=[0.9290    0.6940    0.1250;
    116/255   52/255    129/255;
     0.8500    0.3250    0.0980];
set(groot,'defaultAxesColorOrder',co);
figure, semilogy(1:iter1,cost1(1:iter1),'-','LineWidth',2); hold on
semilogy(1:iter2,cost2(1:iter2),'-','LineWidth',2); hold on
semilogy(1:iter3,cost3(1:iter3),'-','LineWidth',2); hold on
% xlim([1,100])
legend('IRLS','DC-Nuc', 'Proposed DC');
xlabel('Iteration');
ylabel('Cost');
set(gca,'FontSize',14,'FontName','Times New Roman');
ylim([1e-6,inf]);
% set(gca,'YTick',[1e-4,1e-3,1e-2,1e-1,1,10]);

figure, semilogy(time1,cost1(1:iter1),'-','LineWidth',2); hold on
semilogy(time2,cost2(1:iter2),'-','LineWidth',2); hold on
semilogy(time3,cost3(1:iter3),'-','LineWidth',2); hold on
% xlim([0,8])
ylim([1e-6,inf]);
legend('IRLS','DC-Nuc', 'Proposed DC');
xlabel('Time [s]');
ylabel('Cost');
set(gca,'FontSize',14,'FontName','Times New Roman');
%  set(gca,'XTick',[0:1:8]);
%  ylim([1e-4,inf]);