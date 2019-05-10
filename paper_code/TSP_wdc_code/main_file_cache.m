% Code for Figure 4 of paper:  K. Yang, Y. Shi, and Z. Ding, ¡°Data shuffling in wireless distributed
% computing via low-rank optimization,¡± IEEE Trans. Signal Process.,2019.
% Written by Kai Yang
clear all;
rng('default');
rng(1);
K=5;
F = 10;
T=K*F;
d = 1;
cache_size_set =[5:1:F-1];
cache_size_len = length(cache_size_set);
L=1;
M=1;
verb  = 0;

options.verbosity = verb;
options.costtol = 1e-5;
options.ranktol = 1e-5;
options.maxiter = 100;
options_lp = options;
options_lp.maxiter = 10000;
options2.verbosity = verb;
options2.ranktol = 1e-5;
options2.maxiter =10000;
options3.verbosity = verb;
options3.ranktol = 1e-5;
options3.maxiter =30;

testnum = 100;
DoF_vec1 = zeros(cache_size_len,testnum);
DoF_vec2 = zeros(cache_size_len,testnum);
DoF_vec3 = zeros(cache_size_len,testnum);
DoF_vec4 = zeros(cache_size_len,testnum);

mytime_vec1 = zeros(cache_size_len,testnum);
mytime_vec2 = zeros(cache_size_len,testnum);
mytime_vec3 = zeros(cache_size_len,testnum);
mytime_vec4 = zeros(cache_size_len,testnum);
for ci = 1:cache_size_len
    parfor ti = 1:testnum
        H = randn(K,K,L,L);
        [Tset,Rset]=generate_model_file(K,F,cache_size_set(ci));
        [Ax,y,m,n] = generate_parameters(Tset,Rset,L,d,K,T,H);
%         r1=1; r2=1; r3=1;
%         mycost1=0; mycost2=0; mycost3=0;
        [~,r1,mycost1,mytime1]=WDC_nuc( m,n,Ax,y, options);
        DoF_vec1(ci,ti) = d/r1; mytime_vec1(ci,ti) = mytime1;
        try
            [~,r2,mycost2,mytime2 ] = WDC_Lp( m,n,Ax,y,options_lp);
            DoF_vec2(ci,ti) = d/r2; mytime_vec2(ci,ti) = mytime2;
        catch err
            r2 = nan; mycost2 = nan;mytime2=nan;
            DoF_vec2(ci,ti) = nan; mytime_vec2(ci,ti) = nan;
        end
        [~,r3,mycost3,mytime3]=WDC_DCF( m,n,Ax,y, options2);
        DoF_vec3(ci,ti) = d/r3; mytime_vec3(ci,ti) = mytime3;
        
        [~,r4,mycost4,mytime4]=WDC_DC( m,n,Ax,y, options3);
        DoF_vec4(ci,ti) = d/r4; mytime_vec4(ci,ti) = mytime4;
        fprintf('cache size:%d, testnum:%d, DoF:%d,%d,%d,%d cost:%.3e,%.3e,%.3e,%.3e time:%.1e,%.1e,%.1e,%.1e\n',cache_size_set(ci),ti,DoF_vec1(ci,ti),DoF_vec2(ci,ti),DoF_vec3(ci,ti),DoF_vec4(ci,ti),mycost1,mycost2,mycost3,mycost4,mytime1,mytime2,mytime3,mytime4);
    end
end
save('wdc_cache.mat');
co=[ 0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    116/255   52/255    129/255];
set(groot,'defaultAxesColorOrder',co);
figure, plot(cache_size_set,mean(DoF_vec1,2),'*-','LineWidth',2); hold on
plot(cache_size_set,mean(DoF_vec2,2,'omitnan'),'s-','LineWidth',2); 
plot(cache_size_set,mean(DoF_vec3,2),'o-','LineWidth',2); 
plot(cache_size_set,mean(DoF_vec4,2),'d-','LineWidth',2); 
legend('Nuclear norm','IRLS','Proposed DC','DC-Nuc');
xlabel('Local Storage Size');
ylabel('Achievable DoF');
set(gca,'XTick',5:9);
set(gca,'YTick',0:0.1:1);
set(gca,'FontSize',14,'FontName','Times New Roman');


figure, plot(cache_size_set,mean(mytime_vec1,2),'*-','LineWidth',2); hold on
plot(cache_size_set,mean(mytime_vec2,2,'omitnan'),'s-','LineWidth',2); 
plot(cache_size_set,mean(mytime_vec3,2),'o-','LineWidth',2); 
plot(cache_size_set,mean(mytime_vec4,2),'d-','LineWidth',2); 
legend('Nuclear norm','IRLS','Proposed DC','DC-Nuc [13]');
xlabel('Local Storage Size');
ylabel('Time [s]');
set(gca,'XTick',5:9);
% set(gca,'YTick',0:0.1:1);
set(gca,'FontSize',14,'FontName','Times New Roman');