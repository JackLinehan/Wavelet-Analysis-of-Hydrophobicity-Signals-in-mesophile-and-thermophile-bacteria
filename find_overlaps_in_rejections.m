%% Finding the pairs that Reject the KS test 

% hyrophobic sequence results from kstest2 

ES_rejects = find(ks_ES == 1); 

HW_rejects = find(ks_HW == 1); 

KD_rejects = find(ks_KD == 1); 

%% Checking for overlap between hydrophobicity scales 

my_lengths = [abs(length(ES_rejects)), abs(length(HW_rejects)), abs(length(KD_rejects))]; 

my_names = {ES_rejects, HW_rejects, KD_rejects}; 

index_lengths = find(my_lengths == max(my_lengths)); % find the hydro-scale with the most rejects 

my_list = my_names{index_lengths}; % this is the group with the most rejections, and I'll check its values against the others 

my_name = cell(1,length(my_lengths)-1); 

count = 1;

for j = 1:abs(length(my_names))
    
    if (j ~= index_lengths)

        my_name{1,count} = my_names{1,j}; 
        
        count = count + 1; 
        
    end 
    
end 

count = 1; 

for j = 1:abs(length(my_list))
    
    check_1 = isempty(find(my_name{1} == my_list(1,j))); 
    
    check_2 = isempty(find(my_name{2} == my_list(1,j))); 
    
        
    if (check_1 == 0 && check_2 == 0)
        
        pairs_rejecting_ks_all_hydscales(count) = my_list(1,j); 
        disp(pairs_rejecting_ks_all_hydscales(count)); 
        count = count + 1; 
        
    end 
        
end 



