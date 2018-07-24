function [meso_coeff, therm_coeff] = WP_extract(meso,therm, my_lengths)
%% Wavelet Packet analysis, power series approximation 
% This function looks at the scale determined by the shorter of the two
% proteins, and goes two to four levels above that. 
%% Get my lengths so we can compute levels 

%protein_lengths; 
my_min_lengths = min(my_lengths,[],2);  
dwtmode('zpd'); 


meso_coeff = {1,540}; 

therm_coeff = {1,540}; 

for j = 1:abs(length(meso))
    
    level = wmaxlev(my_min_lengths(j,1),'db1'); 
    
    WPT_meso = wpdec(meso{1,j},level-2,'db1'); 
    
    WPT_therm = wpdec(therm{1,j},level-2,'db1'); 
    
    [~, T_meso_seq, ~,J_meso] = otnodes(WPT_meso,'dp'); 
    
    [~, T_therm_seq, ~,J_therm] = otnodes(WPT_therm,'dp'); 
    
    coeff_meso = zeros(1,abs(length(T_meso_seq))); 
    coeff_therm = zeros(1,abs(length(T_therm_seq))); 
    
    J_meso = J_meso -1; 
    J_therm = J_therm - 1; 
    
    for k = 1:abs(length(J_meso))
        
        coeff_meso(1,k) = mean(wpcoef(WPT_meso,[level-2,J_meso(k)]).^2); 
        
        coeff_therm(1,k) = mean(wpcoef(WPT_therm,[level-2, J_therm(k)]).^2); 
        
    end 
    
    meso_coeff{1,j} = coeff_meso; 
    
    therm_coeff{1,j} = coeff_therm; 
    
    
end 