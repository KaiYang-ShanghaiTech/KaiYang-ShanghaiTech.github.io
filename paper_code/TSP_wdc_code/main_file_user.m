% Code for Figure 6 of paper:  K. Yang, Y. Shi, and Z. Ding, ¡°Data shuffling in wireless distributed
% computing via low-rank optimization,¡± IEEE Trans. Signal Process.,2019.
% Written by Kai Yang
clear all;
rng('default');
rng(1);
F = 5;
d = 1;
cache_len = 2;
L = 1;
user_size_set = [5:5:20];
user_size_len = length(user_size_set);
verb =0;
options.maxiter = 100;
options.verbosity = verb;
options.costtol = 1e-5;
options.ranktol = 1e-5;
options_lp = options;
options_lp.maxiter = 10000;
options2 = options;
options2.maxiter =10000;
testnum = 100;
DoF_vec1 = zeros(user_size_len,testnum);
DoF_vec2 = zeros(user_size_len,testnum);
DoF_vec3 = zeros(user_size_len,testnum);
for ci = 1:user_size_len
    K = user_size_set(ci);
    T = F*K;
    parfor ti = 1:testnum
        [Tset,Rset]=generate_model_file_even(K,F,cache_len);
        H = randn(K,K,L,L);
        [Ax,y,m,n] = generate_parameters(Tset,Rset,L,d,K,T,H);
%         r1=1; r2=1; r3=1;
%         mycost1=0; mycost2=0; mycost3=0; mycost4=0;
        [~,r1,mycost1]=WDC_nuc( m,n,Ax,y, options);
        DoF_vec1(ci,ti) = d/r1;
        try
            [~,r2,mycost2]=WDC_Lp( m,n,Ax,y, options_lp);
            DoF_vec2(ci,ti) = d/r2;
        catch err
            r2 = nan; mycost2 = nan;
            DoF_vec2(ci,ti) = nan;
        end
        [~,r3,mycost3]=WDC_DCF( m,n,Ax,y, options2);
        DoF_vec3(ci,ti) = d/r3;
        fprintf('Number of users:%d, testnum:%d, DoF:%d,%d,%d, cost:%.3e,%.3e,%.3e\n',K,ti,DoF_vec1(ci,ti),DoF_vec2(ci,ti),DoF_vec3(ci,ti),mycost1,mycost2,mycost3);
    end
end
save('wdc_users3.mat');
co=[ 0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250];
set(groot,'defaultAxesColorOrder',co);
figure, plot(user_size_set,mean(DoF_vec1,2),'*-','LineWidth',1.2); hold on
plot(user_size_set,mean(DoF_vec2,2,'omitnan'),'s-','LineWidth',1.2); hold on
plot(user_size_set,mean(DoF_vec3,2),'o-','LineWidth',1.2);hold on
legend('Nuclear norm','IRLS','Proposed DC');
xlabel('Number of Users');
ylabel('Achievable DoF');
set(gca,'FontSize',14);
ylim([0,0.3]);