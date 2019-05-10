% Code for Figure 5 of paper:  K. Yang, Y. Shi, and Z. Ding, ¡°Data shuffling in wireless distributed
% computing via low-rank optimization,¡± IEEE Trans. Signal Process.,2019.
% Written by Kai Yang
clear all;
rng('default');
rng(1);
K=8;
F = 4;
T = F*K;
d = 1;
cache_len = 1;
antenna_size_set = [1:1:4];
antenna_size_len = length(antenna_size_set);
verb = 0;
options.maxiter = 100;
options.verbosity = verb;
options.ranktol = 1e-5;
options.costtol = 1e-5;
options_lp = options;
options_lp.maxiter = 10000;
options2 = options;
options2.maxiter =10000;
testnum = 100;
DoF_vec1 = zeros(antenna_size_len,testnum);
DoF_vec2 = zeros(antenna_size_len,testnum);
DoF_vec3 = zeros(antenna_size_len,testnum);

[Tset,Rset]=generate_model_file(K,F,cache_len);
tstart=tic;
for ci = 1:antenna_size_len
    L = antenna_size_set(ci);
    parfor ti = 1:testnum
        H = randn(K,K,L,L);
        [Ax,y,m,n] = generate_parameters(Tset,Rset,L,d,K,T,H);
        r1=20;r2=20;r3=20;
        [~,r1,mycost1]=WDC_nuc( m,n,Ax,y, options);
        DoF_vec1(ci,ti) = d/r1;
        try
            [~,r2,mycost2]=WDC_Lp( m,n,Ax,y, options_lp);
            DoF_vec2(ci,ti) = d/r2;
        catch
            DoF_vec2(ci,ti) = nan;
            mycost2 = nan;
        end
        [~,r3,mycost3]=WDC_DCF( m,n,Ax,y, options2);
        DoF_vec3(ci,ti) = d/r3;
        fprintf('Number of antennas:%d, testnum:%d, DoF:%d,%d,%d, cost:%.3e,%.3e,%.3e\n',antenna_size_set(ci),ti,DoF_vec1(ci,ti),DoF_vec2(ci,ti),DoF_vec3(ci,ti),mycost1,mycost2,mycost3);
        
    end
end
mytime = toc(tstart);
save('wdc_antenna.mat');
co=[ 0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250];
set(groot,'defaultAxesColorOrder',co);
figure, plot(antenna_size_set,mean(DoF_vec1,2),'*-','LineWidth',1.5); hold on
plot(antenna_size_set,mean(DoF_vec2,2,'omitnan'),'s-','LineWidth',1.5);hold on
plot(antenna_size_set,mean(DoF_vec3,2),'o-','LineWidth',1.5);hold on
legend('Nuclear norm','IRLS','Proposed DC');
xlabel('Number of Antennas');
ylabel('Achievable DoF');
set(gca,'FontSize',14,'FontName','Times New Roman');
set(gca,'XTick',1:4);
ylim([0,0.3]);