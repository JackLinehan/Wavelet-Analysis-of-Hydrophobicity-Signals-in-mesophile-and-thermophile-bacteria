function [meso_coeff, therm_coeff] = WP_extract(meso,therm)
%% Wavelet Packet analysis, power series approximation 
%% Get my lengths so we can compute levels 

protein_lengths; 
my_min_lengths = min(my_lengths,[],2);  
dwtmode('per'); 


meso_coeff = {1,540}; 

therm_coeff = {1,540}; 

for j = 1:abs(length(meso))
    
    level = wmaxlev(my_min_lengths(j,1),'db1'); 
    
    WPT_meso = wpdec(meso{1,j},level-4,'db1'); 
    
    WPT_therm = wpdec(therm{1,j},level-4,'db1'); 
    
    [~, T_meso_seq] = otnodes(WPT_meso); 
    
    [~, T_therm_seq] = otnodes(WPT_therm); 
    
    for k = 1:abs(length(T_meso_seq))
        
        coeff_meso = mean(wpcoef(WPT_meso,[level,T_meso_seq(k)]).^2); 
        
        coeff_therm = mean(wpcoef(WPT_therm,[level, T_therm_seq(k)]).^2); 
        
    end 
    
end 