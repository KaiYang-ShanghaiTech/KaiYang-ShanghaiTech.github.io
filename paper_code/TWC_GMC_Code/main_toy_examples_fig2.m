% Code of Fig 2 in paper: K. Yang, Y. Shi, and Z. Ding, ¡°Generalized low-rank
% optimization for topological cooperation in ultra-dense networks,¡± IEEE
% Trans. Wireless Commun., 2019
% Written by Kai Yang
clear;addpath('./func');%addpath(genpath('./online_toolbox'));
%% This function recovers the result in 
% Yi, Xinping, and David Gesbert. "Topological interference management with transmitter cooperation." IEEE Transactions on Information Theory 61.11 (2015): 6107-6130.
rng('default'); rng(1);
receiver_cache_size = 0;
params.verbosity =0;
params.costtol = 1e-5;
params.maxiter = 1000;
params.tolgradnorm = 1e-8;
options_altmin = params;
options_altmin.maxiter = 50;
options_altmin.inmaxiter = 1000;
optionsCG = params;
optionsCG.maxiter = 1e3;

flag_altmin = 0;  flag_tr = 1;  flag_cg = 1;
flag_exam1 = 1; flag_exam4 = 1; flag_exam6 = 1; flag_exam7 = 1;
%%--------------------------------------Example in Fig. 2--------------------------------------
    K = 5;
    connected_links = [1,0,1,1,0;
                                    0,1,1,1,0;
                                    1,0,1,0,1;
                                    1,0,0,1,1;
                                    0,1,0,0,1];
    
    % non cooperation optimal achievable DoF 2/5
    data_stream = 1;
    trans_cache_size = 0;
    [trans_cache_set, receiver_cache_set] = make_transeiver_cache_size(K,trans_cache_size, receiver_cache_size);
    [A,b,m,n] = generate_matrix_multidatastream(connected_links,trans_cache_set,receiver_cache_set,data_stream);
    if flag_altmin
        [~,rank_altmin,infos_altmin] = topological_beamforming_altmin( [m,n],[],options_altmin,A,b);
        all_time_altmin = infos_altmin.time; all_rank_altmin=rank_altmin;
        fprintf(' AltMin: DoF %d',data_stream/all_rank_altmin);
    end
    if flag_cg
        [X_cg,rank_CG,infos_CG ] = topological_beamforming_CG( [m,n],[], optionsCG,A,b );
        all_time_CG = infos_CG.time; all_rank_CG=rank_CG;
        fprintf(' CG: DoF %d Time %.1f',data_stream/all_rank_CG,all_time_CG);
    end
    if flag_tr
        [X_tr,rank_TR,infos_TR ] = topological_beamforming_TR( [m,n],[], params,A,b );
        all_time_TR = infos_TR.time; all_rank_TR=rank_TR;
        fprintf(' TR: DoF %d Time %.1f',data_stream/all_rank_TR,all_time_TR);
    end
    fprintf('\n');
% full cooperation achievable DoF 1/2
    data_stream = 1;
    receiver_cache_set = zeros(K,K);
    trans_cache_set=connected_links;
    [A,b,m,n] = generate_matrix_multidatastream(connected_links,trans_cache_set,receiver_cache_set,data_stream);
    if flag_altmin
        [~,rank_altmin,infos_altmin] = topological_beamforming_altmin( [m,n],[],options_altmin,A,b);
        all_time_altmin = infos_altmin.time; all_rank_altmin=rank_altmin;
        fprintf(' AltMin: DoF %d',data_stream/all_rank_altmin);
    end
    if flag_cg
        [ ~,rank_CG,infos_CG ] = topological_beamforming_CG( [m,n],[], optionsCG,A,b );
        all_time_CG = infos_CG.time; all_rank_CG=rank_CG;
        fprintf(' CG: DoF %d Time %.1f',data_stream/all_rank_CG,all_time_CG);
    end
    if flag_tr
        [ ~,rank_TR,infos_TR ] = topological_beamforming_TR( [m,n],[], params,A,b );
        all_time_TR = infos_TR.time; all_rank_TR=rank_TR;
        fprintf(' TR: DoF %d Time %.1f',data_stream/all_rank_TR,all_time_TR);
    end
    fprintf('\n');
% save('example.mat');