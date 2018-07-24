%% Generate Figures of qqplots for power series 

fname = 'qqplots for the pairs Power Series Bootstrap; wmaxlev - 4'; 
mkdir(fname); 
mdir = pwd; 
fullfig; 
name = 'QQ Plot of Sample Data versus Standard Normal'; 
for j = 1:540
    
    figure(clf); 
    
    subplot(3,2,1); 
    qqplot(KD_meso_C_shuffle{1,j}); 
    title([name, ' - Mesophile, KD scale']); 
    subplot(3,2,3); 
    qqplot(HW_meso_C_shuffle{1,j}); 
    title([name, ' - Mesophile, HW scale']);
    subplot(3,2,5); 
    qqplot(ES_meso_C_shuffle{1,j}); 
    title([name, ' - Mesophile, ES scale']);
    
    
    subplot(3,2,2); 
    qqplot(KD_therm_C_shuffle{1,j}); 
    title([name, ' - Thermophile, KD scale']); 
    subplot(3,2,4); 
    qqplot(HW_therm_C_shuffle{1,j}); 
    title([name, ' - Thermophile, HW scale']);
    subplot(3,2,6); 
    qqplot(ES_therm_C_shuffle{1,j}); 
    title([name, ' - Thermophile, ES scale']);
    
    cd(fname); 
    
    page_name = ['QQ plot for bootstrapped pair' num2str(j)]; 
    
    print(page_name,'-dtiff'); 
    
    cd(mdir); 
    
end 