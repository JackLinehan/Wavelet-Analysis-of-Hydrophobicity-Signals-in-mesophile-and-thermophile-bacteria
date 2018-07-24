%% Make historgrams 

% r = rand(1,2); 
% 
% r = round(r.*540 + 1); 
% 
% 
% 
% for j = 1:abs(length(r))
%   
%     hold on 
%     
%     histogram(KD_meso_C{1,r(j)}); 
%     
%     histogram(KD_therm_C{1,r(j)}); 
%     
%     hold off 
%     
%     
% end 
    
function [] = make_histograms(j,KD_meso_C,KD_therm_C)

hold on

plot(KD_meso_C{1,j},'+'); 

plot(KD_therm_C{1,j}); 

hold off 



end 
    