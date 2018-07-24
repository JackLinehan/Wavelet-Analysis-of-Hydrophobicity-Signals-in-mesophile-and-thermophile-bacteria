function [output_m, output_t] = bootstrap(input_m, input_t) 
%% This function uses the bootstrapping method to shuffle the hydrophobicity values along the signal 

cell_num = length(input_m); 

output_m = {1,cell_num}; 

output_t = {1,cell_num}; 

for j = 1:cell_num
    
    meso_seq = input_m{1,j}; 
    
    meso_seq_sort = sort(meso_seq); % Change the initial distribution
    
    therm_seq = input_t{1,j}; 
    
    therm_seq_sort = sort(therm_seq); % Change the initial distribution
    
    count = 0; 
    
    % I'd like to use a measure to determine this, matlab uses merinee
    % twister the generate PRNs. Not sure what the optimal number of
    % shuffles is though. 
    
    while (count < 1000) % shuffle it so many times that every permuatation is equally likely to occur
        
        
       meso_seq_sort = meso_seq_sort(randperm(length(meso_seq))); 
       
       therm_seq_sort = therm_seq_sort(randperm(length(therm_seq))); 
        
       count = count +1; 
        
    end 
    
    output_m{1,j} = meso_seq_sort; 
    
    output_t{1,j} = therm_seq_sort; 


    
end 