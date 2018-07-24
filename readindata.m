function [ES_meso, ES_therm, HW_meso, HW_therm, KD_meso, KD_therm] = readindata()
%% Read data into matlab and make histograms 
tic; 
file_meso_ID = fopen('540_Mesophiles.txt'); 

file_therm_ID = fopen('540_Thermophiles.txt'); 

mesophiles = textscan(file_meso_ID, '%q'); 
mesophiles = mesophiles{1,1}; 

thermophiles = textscan(file_therm_ID, '%q'); 
thermophiles = thermophiles{1,1}; 

% So, yes, there are only 20 amino acids, but PDB files aren't perfect and
% we have a handfull of instances out of hundreds of thousands that are
% incorrect entries. 
key_set = {'A', 'C', 'D', 'E', 'F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X','U'}; 
KD_num = [ 1.8,2.5 , -3.5, -3.5, 2.8, -0.40,  -3.20, 4.5, -3.9, 3.8, 1.9, -3.5,  -1.6,  -3.5,  -4.5,  -0.8,  -0.7,  4.20, -0.9, -1.3,0,0]; 
HW_num = [ -0.5,  -1.0,  3.0, 3.0, -2.5, 0.0, -0.5, -1.8, 3.00, -1.8, -1.3, 0.2, 0.0, 0.20,  3.00, 0.3, -0.4, -1.5, -3.40, -2.3,0,0]; 
ES_num = [ 1.6, 2.0, -9.2, -8.2, 3.7, 1, -3.0, 3.10,-8.80, 2.8, 3.4, -4.8, -0.2, -4.1, -12.3, 0.60, 1.20, 2.60,  1.9, -0.7,0,0]; 

KD_container = containers.Map(key_set,KD_num); 
HW_container = containers.Map(key_set,HW_num); 
ES_container = containers.Map(key_set,ES_num); 

KD_meso = {1,abs(length(mesophiles))}; 
HW_meso = {1,abs(length(mesophiles))}; 
ES_meso = {1,abs(length(mesophiles))}; 

for j = 1:abs(length(mesophiles))
    
    my_seq = mesophiles{j,1}; 
    
    kd_seq = zeros(1,abs(length(my_seq))); 
    hw_seq = zeros(1,abs(length(my_seq))); 
    es_seq = zeros(1,abs(length(my_seq))); 
    
    for k = 1:abs(length(my_seq))
        
       kd_seq(1,k) = KD_container(my_seq(1,k)); 
       hw_seq(1,k) = HW_container(my_seq(1,k)); 
       es_seq(1,k) = ES_container(my_seq(1,k)); 
        
    end 
    
    KD_meso{1,j} = kd_seq; 
    HW_meso{1,j} = hw_seq; 
    ES_meso{1,j} = es_seq; 
    
end 


KD_therm = {1,abs(length(thermophiles))}; 
HW_therm = {1,abs(length(thermophiles))}; 
ES_therm = {1,abs(length(thermophiles))}; 

for j = 1:abs(length(thermophiles))
    
    my_seq = thermophiles{j,1}; 
    
    kd_seq = zeros(1,abs(length(my_seq))); 
    hw_seq = zeros(1,abs(length(my_seq))); 
    es_seq = zeros(1,abs(length(my_seq))); 
    
    for k = 1:abs(length(my_seq))

       kd_seq(1,k) = KD_container(my_seq(1,k)); 
       hw_seq(1,k) = HW_container(my_seq(1,k)); 
       es_seq(1,k) = ES_container(my_seq(1,k)); 
        
    end 
    
    KD_therm{1,j} = kd_seq; 
    HW_therm{1,j} = hw_seq; 
    ES_therm{1,j} = es_seq; 
    
end 
clear vars ES_container ES_num es_seq file_meso_ID file_therm_ID HW_container
clear vars HW_num hw_seq j k KD_container KD_num kd_seq key_set mesophiles my_seq thermophiles
toc; 
end 