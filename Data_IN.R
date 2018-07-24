# Load in amino acid chain sequences for meso and thermo bacteria and generate hydrophobicity signals 
library(tictoc)
tic()
meso_seq <- strsplit(scan("540_Mesophiles.txt", what = "list", sep = " "), split = NULL) 

therm_seq <- strsplit(scan("540_Thermophiles.txt", what = "list", sep = " "), split = NULL) 

# and now I have two lists of 540 lists, with accessible elements WOOT! Scan reads in the data from the txt files as lists, then stringsplit comes in and seperates each entry into vectors, so that we can access the amino acids in each protein. structure looks like this: list(protein(list(aminoacid)))

######################################################################################################################################
# Here I define some protein names for an index that was going to be used for a much smaller data set. These are no longer relevant. 

#meso_protein_name <- c("ACP_2", "ACP", "Acyl", "ADK_2", "ADK", "Chemotaxis", "CheW", "CheY", "ColdShock_3","ColdShock", "ColdShock_2", "FKBP", "HGCS", "HPr", "IF1A", "IF5A", "NusB", "RNaseH", "RNaseP", "RRF", "Sigma", "Thioredoxin", "TIF3", "UNK")

#therm_protein_name <- c("ACP_2", "ACP", "Acyl", "ADK_2", "ADK", "Chemotaxis", "CheW", "CheY", "ColdShock_2", "ColdShock_3", "ColdShock", "FKBP", "HGCS", "Hpr", "IF1A", "IF5A", "NusB", "RNaseH", "RNaseP", "RRF", "Sigma", "Thioredoxin", "TIF3", "UNK")
######################################################################################################################################
# Now let's generate some hydrophobicity scales: this is done by generating a "pseudo"- dictionary in R that mimicks the dictionary in python, by making a list: as well I generate numerical scales and find the mean values

KD <- list( A = 1.8, C = 2.5 , D = -3.5, E = -3.5, F = 2.8, G = -0.40, H = -3.20, I = 4.5, K = -3.9, L = 3.8, M = 1.9, N = -3.5, P = -1.6, Q = -3.5, R = -4.5, S = -0.8, T = -0.7, V = 4.20, W = -0.9, Y = -1.3)

KD_num <- c( 1.8,2.5 , -3.5, -3.5, 2.8, -0.40,  -3.20, 4.5, -3.9, 3.8, 1.9, -3.5,  -1.6,  -3.5,  -4.5,  -0.8,  -0.7,  4.20, -0.9, -1.3)

mean_kd <- mean(KD_num)

HW <- list( A = -0.5, C = -1.0, D = 3.0, E = 3.0, F = -2.5, G = 0.0, H = -0.5, I = -1.8, K = 3.00, L = -1.8, M = -1.3, N = 0.2, P = 0.0, Q = 0.20, R = 3.00, S = 0.3, T = -0.4, V = -1.5, W = -3.40, Y = -2.3)
  
  
HW_num <- c( -0.5,  -1.0,  3.0, 3.0, -2.5, 0.0, -0.5, -1.8, 3.00,  -1.8, -1.3, 0.2, 0.0, 0.20,  3.00, 0.3, -0.4,  -1.5, -3.40,  -2.3)

mean_hw <- mean(HW_num)
  

ES <- list(A = 1.6, C = 2.0, D = -9.2, E = -8.2, F = 3.7, G = 1, H = -3.0, I = 3.10, K = -8.80, L = 2.8, M = 3.4, N = -4.8, P = -0.2, Q = -4.1, R = -12.3, S = 0.60, T = 1.20, V = 2.60, W = 1.9, Y = -0.7)


ES_num <- c( 1.6, 2.0, -9.2, -8.2, 3.7, 1, -3.0, 3.10,-8.80, 2.8, 3.4, -4.8, -0.2,-4.1, -12.3, 0.60, 1.20, 2.60,  1.9,-0.7)

mean_ES <- mean(ES_num)

# This first loop generates signals for mesophile bacteria 

KD_meso <- matrix(list(), 1, 540)
HW_meso <- matrix(list(), 1, 540)
ES_meso <- matrix(list(), 1, 540)


for (j in 1:length(meso_seq)) { 
	
	x <- meso_seq[j]
	my_seq <- unlist(x)
	num_residues <- length(my_seq)
	KD_seq <- numeric(num_residues) 
	HW_seq <- numeric(num_residues)
	ES_seq <- numeric(num_residues)
	
	for (k in 1:num_residues){
		
		KD_seq[k] <- KD[my_seq[k]]
		HW_seq[k] <- HW[my_seq[k]]
		ES_seq[k] <- ES[my_seq[k]]

	}
	
	KD_meso[[1,j]] = unlist(KD_seq)
	HW_meso[[1,j]] = unlist(HW_seq)
	ES_meso[[1,j]] = unlist(ES_seq)
	
	
	}
	
	
	# And now for the thermophiles 
	
	KD_therm <- matrix(list(), 1, 540)
    HW_therm <- matrix(list(), 1, 540)
    ES_therm <- matrix(list(), 1, 540)


for (j in 1:length(therm_seq)) { 
	
	x <- therm_seq[j]
	my_seq <- unlist(x)
	num_residues <- length(my_seq)
	KD_seq <- numeric(num_residues) 
	HW_seq <- numeric(num_residues)
	ES_seq <- numeric(num_residues)
	
	for (k in 1:num_residues){
		
		KD_seq[k] <- KD[my_seq[k]]
		HW_seq[k] <- HW[my_seq[k]]
		ES_seq[k] <- ES[my_seq[k]]

	}
	
	KD_therm[[1,j]] = unlist(KD_seq)
	HW_therm[[1,j]] = unlist(HW_seq)
	ES_therm[[1,j]] = unlist(ES_seq)
	
	
	}
	toc()
	
	# So now we have the following data sets: 
	# meso_seq, therm_seq - list of proteins with character sequences 
	# KD_therm/meso, HW_therm/meso, ES_therm/meso, matricies of proteins with vectors of primary chains  
	
unique_meso <- sapply(KD_meso, function(x) length(unique(x)))
unique_therm <- sapply(KD_therm, function(x) length(unique(x)))
unique_meso/unique_therm
	
# playing with overlapping histograms 

#meso <- unlist(KD_meso[1])
#therm <- unlist(KD_therm[1])

#hist(meso, breaks = 17, col = rgb(1,0,0,0.5), main = "overlapping Histogram")
#hist(therm, breaks = 17, col = rgb(0,0,1,0.5), add = T)
#box()

## Mean and Variance -----------------------------------------------------------------------------------------------------------------------------

# Kyte-Doolittle 
mean_kd_meso <- sapply(KD_meso,mean)

var_kd_meso <- sapply(KD_meso,var)

mean_kd_therm <- sapply(KD_therm,mean)

var_kd_meso <- sapply(KD_therm,var)

# Hopp Woods 
mean_HW_meso <- sapply(HW_meso,mean)

var_HW_meso <- sapply(HW_meso,var)

mean_HW_therm <- sapply(HW_therm,mean)

var_HW_meso <- sapply(HW_therm,var)

#Engleman-Steitz 
mean_ES_meso <- sapply(ES_meso,mean)

var_ES_meso <- sapply(ES_meso,var)

mean_ES_therm <- sapply(ES_therm,mean)

var_ES_meso <- sapply(ES_therm,var)

# ratio of the means 

ratio_kd_meso <- mean_kd_meso/mean_kd 
ratio_kd_therm <- mean_kd_therm / mean_kd 

ratio_hw_meso <- mean_HW_meso / mean_hw
ratio_hw_therm <- mean_HW_therm / mean_hw 

ratio_ES_meso <- mean_ES_meso / mean_ES
ratio_ES_therm <- mean_ES_therm / mean_ES

plot(ratio_kd_meso,xlab = "indexing",ylab = "ratio of mean value per protien with KD mean",main = "KD mesophile mean/ mean(KD)")
lines(rep(1,540))

dev.new()
plot(ratio_kd_therm,xlab = "indexing", ylab = "ratio of mean value per protien with KD mean", main = "KD thermophile mean / mean(KD)" )
lines(rep(1,540))

dev.new()
plot(ratio_hw_meso, xlab = "indexing", ylab = "ratio of mean value per protien with HW mean", main = "HW mesophile mean/ mean(HW)")
lines(rep(1,540))

dev.new()
plot(ratio_hw_therm, xlab = "indexing",ylab = "ratio of mean value per protien with HW mean",main = "HW thermophile mean/ mean(HW)" )
lines(rep(1,540))


dev.new()
plot(ratio_ES_meso, xlab = "indexing",ylab = "ratio of mean value per protien with ES mean",main = "ES mesohpile mean/ mean(ES)")
lines(rep(1,540))

dev.new()
plot(ratio_ES_therm,xlab = "indexing",ylab = "ratio of mean value per protien with ES mean",main = " ES thermophile mean/ mean(ES)" )
lines(rep(1,540))


## KS Test (two sample)------------------------------------------------------------------------------------------------------------------------
#Kyte-Doolittle 

KS_test_KD <- matrix(list(), 2, 540)

for (j in 1:540){ 
	
	pop = ks.test(unlist(KD_therm[j]), unlist(KD_meso[j]))
	KS_test_KD[1,j] = pop$statistic
	KS_test_KD[2,j] = pop$p.value
	}

p_values <- KS_test_KD[2,]

p_0.5_KD_KS = with(p_values, c(sum(p_values <=0.05))) # p values from kyte-doolittle for KS test - number of stat sig values 
# with applies an expression to a dataset. with(data, expression), in our case we go through and check if a value is greater than or less than
# 0.05, this returns true (1) or false(0). Then we sum up all the values and get the total number of p-values that are stat sig. 

#Hopp-Woods 

KS_test_HW <- matrix(list(), 2, 540)

for (j in 1:540){ 
	
	pop = ks.test(unlist(HW_therm[j]), unlist(HW_meso[j]))
	KS_test_HW[1,j] = pop$statistic
	KS_test_HW[2,j] = pop$p.value
	}

p_values_HW_KS <- KS_test_HW[2,]

p_0.5_HW_KS = with(p_values_HW_KS, c(sum(p_values_HW_KS <=0.05)))

# Engleman-Steitz 

KS_test_ES <- matrix(list(), 2, 540)

for (j in 1:540){ 
	
	pop = ks.test(unlist(ES_therm[j]), unlist(ES_meso[j]))
	KS_test_ES[1,j] = pop$statistic
	KS_test_ES[2,j] = pop$p.value
	}

p_values_ES_KS <- KS_test_ES[2,]

p_0.5_ES_KS = with(p_values_ES_KS, c(sum(p_values_ES_KS <=0.05)))
#----------------------------------------------------------------------------------------------------------------------------------------------
# Anderson-Darling
# Do normality first, this makes sense
library(nortest)
ad_meso_kd <- sum(sapply(KD_meso, function(x) ad.test(unlist(x))$p.value)>= 0.05)
ad_therm_kd <- sum(sapply(KD_therm, function(x) ad.test(unlist(x))$p.value)>= 0.05)

ad_meso_hw <- sum(sapply(HW_meso, function(x) ad.test(unlist(x))$p.value)>=0.05)
ad_therm_hw <- sum(sapply(HW_therm, function(x) ad.test(unlist(x))$p.value)>=0.05)

ad_meso_es <- sum(sapply(ES_meso, function(x) ad.test(unlist(x))$p.value)>=0.05)
ad_therm_es <- sum(sapply(ES_therm, function(x) ad.test(unlist(x))$p.value)>0.05)

# None of these are normally distributed 
#---------------------------------------------------------------------------------------------------------------------------------------------
# Lets generate some Histograms!!!!! 


	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
