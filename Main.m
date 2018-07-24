%% Main: call em' all 
tic; 
%% Experimetnal Group 
[ES_meso, ES_therm, HW_meso, HW_therm, KD_meso, KD_therm] = readindata(); 

% Check hydrophobicity distirbutions of each pair for similarit using
% KStest2

%[ks_KD] = check_KS(KD_meso, KD_therm); 

%[ks_HW] = check_KS(HW_meso, HW_therm); 

%[ks_ES] = check_KS(ES_meso, ES_therm); 


[my_lengths] = read_lengths_proteins(KD_meso,KD_therm); 

% Determine the Power Series values of the wavelet packet decomposition

[KD_meso_C, KD_therm_C] = WP_extract(KD_meso, KD_therm,my_lengths); 

[HW_meso_C, HW_therm_C] = WP_extract(HW_meso, HW_therm, my_lengths); 

[ES_meso_C, ES_therm_C] = WP_extract(ES_meso, ES_therm,my_lengths); 

% Perform KStest2 on each pair's powerseries 
%[ks_KD_C] = check_KS(KD_meso_C, KD_therm_C); 

%[ks_HW_C] = check_KS(HW_meso_C, HW_therm_C); 

%[ks_ES_C] = check_KS(ES_meso_C, ES_therm_C); 



[KD_meso_leaf_skew, KD_therm_leaf_skew] = WP_leaves_band_skew(KD_meso, KD_therm,my_lengths); 

[HW_meso_leaf_skew, HW_therm_leaf_skew] = WP_leaves_band_skew(HW_meso, HW_therm,my_lengths); 

[ES_meso_leaf_skew, ES_therm_leaf_skew] = WP_leaves_band_skew(ES_meso, ES_therm,my_lengths); 


%% Control Group 

% Generate control data set- maintain frequency remove position dependence 
% 
% [KD_meso_shuffle, KD_therm_shuffle] = bootstrap(KD_meso, KD_therm); 
% 
% [HW_meso_shuffle, HW_therm_shuffle] = bootstrap(HW_meso, HW_therm); 
% 
% [ES_meso_shuffle, ES_therm_shuffle] = bootstrap(ES_meso, ES_therm);
% 
% % KS test 2 on each control pair 
% 
% %[ks_KD_shuffle] = check_KS(KD_meso_shuffle, KD_therm_shuffle); 
% 
% %[ks_HW_shuffle] = check_KS(HW_meso_shuffle, HW_therm_shuffle); 
% 
% %[ks_ES_shuffle] = check_KS(ES_meso_shuffle, ES_therm_shuffle); 
% 
% % Compute power series using wavelet packet decomposition 'db1'
% 
% [KD_meso_C_shuffle, KD_therm_C_shuffle] = WP_extract(KD_meso_shuffle, KD_therm_shuffle, my_lengths); 
% 
% [HW_meso_C_shuffle, HW_therm_C_shuffle] = WP_extract(HW_meso_shuffle, HW_therm_shuffle, my_lengths); 
% 
% [ES_meso_C_shuffle, ES_therm_C_shuffle] = WP_extract(ES_meso_shuffle, ES_therm_shuffle, my_lengths); 
% 
% % KS Test 2 on power series 
% %[ks_KD_C_shuffle] = check_KS(KD_meso_C_shuffle, KD_therm_C_shuffle); 
% 
% %[ks_HW_C_shuffle] = check_KS(HW_meso_C_shuffle, HW_therm_C_shuffle); 
% 
% %[ks_ES_C_shuffle] = check_KS(ES_meso_C_shuffle, ES_therm_C_shuffle); 
% 
% % Determine which pairs rejected the null, across all three hydrophobicity
% % scales: for reference there are 5 pairs that rejected the null: 
% % [107, 110, 379, 388, 531] 
% 
% %find_overlaps_in_rejections; 
% 
% 
% [KD_meso_contorl_leaf_skew, KD_therm_control_leaf_skew] = WP_leaves_band_skew(KD_meso_shuffle, KD_therm_shuffle,my_lengths); 
% 
% [HW_meso_control_leaf_skew, HW_therm_control_leaf_skew] = WP_leaves_band_skew(HW_meso_shuffle, HW_therm_shuffle,my_lengths); 
% 
% [ES_meso_control_leaf_skew, ES_therm_control_leaf_skew] = WP_leaves_band_skew(ES_meso_shuffle, ES_therm_shuffle,my_lengths); 
% 

clear vars j my_list my_name my_names my_lenghts check_1 check_2 count
clear vars HW_rejects KD_rejects index_lengths ES_rejects 
toc; 

