function [my_results] = check_KS( meso, therm) 
%% I pass the signals through the KS test and check the results 

my_results = zeros(1,abs(length(meso))); 

for j = 1:abs(length(meso))
    
    my_results(1,j) = kstest2(meso{1,j},therm{1,j}); 

end 

