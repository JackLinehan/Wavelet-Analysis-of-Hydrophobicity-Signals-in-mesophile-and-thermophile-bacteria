%% Generate Files: I want to actually look at my data 

fname = 'db1 power series overlay plots'; 

mkdir(fname); 

mdir = pwd; 
fullfig; 
for j = 1:540 
    
    figure(clf); 
    
    % Kyte Doolittle Scale 
    
    meso = KD_meso_C{1,j}; 
    
    therm = KD_therm_C{1,j}; 
    
    
    subplot(3,1,1); 
    hold on; 
    
    plot(meso,'+-'); plot(therm,'*-'); 
    title('Kyte-Doolittle Scale: Mesophile and Thermophile Overlay'); 
    xlabel('Time'); 
    ylabel('Power'); 
    legend('mesophile +','thermophile *'); 
    hold off; 
    
   % Hopp-Woods Scale 
   
    meso = HW_meso_C{1,j}; 
    
    therm = HW_therm_C{1,j}; 
    
    
    subplot(3,1,2); 
    hold on; 
    
    plot(meso,'+-'); plot(therm,'*-'); 
    title('Hopp-Woods Scale: Mesophile and Thermophile Overlay'); 
    xlabel('Time'); 
    ylabel('Power');
    legend('mesophile +','thermophile *'); 
    hold off; 
    
    % Engelman-Steiz 
    
    meso = ES_meso_C{1,j}; 
    
    therm = ES_therm_C{1,j}; 
    
    
    subplot(3,1,3); 
    hold on; 
    
    plot(meso,'+-'); plot(therm,'*-'); 
    title('Engelman-Steitz Scale: Mesophile and Thermophile Overlay'); 
    xlabel('Time'); 
    ylabel('Power');
    
    hold off; 
    
    legend('mesophile +','thermophile *'); 
    
    cd(fname)
    
    page_name = [num2str(j) 'Pair Power series']; 
    
    print(page_name,'-dtiff'); 
    
    cd(mdir); 
    
    
end 