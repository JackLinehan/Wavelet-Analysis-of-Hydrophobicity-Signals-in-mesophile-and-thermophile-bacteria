global ES_meso ES_therm HW_meso HW_therm KD_meso KD_therm
my_lengths = zeros(abs(length(KD_meso)),2); 
for j = 1: abs(length(KD_meso))
    
    my_lengths(j,1:2) = [length(KD_meso{1,j}), length(KD_therm{1,j})]; 
    
end 
    
% The smallest mesophile length is 49 amino acids
% The smallest thermophile length is 50 amino acids 

% The longest mesophile amino acid chain is 591 proteins long 
% The longest thermophile amino acid chain is 594 



