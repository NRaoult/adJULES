#library('R.utils')
## code for replacing xoff with x1 in control file     
#
#site <- 'US-Ha1'
#pft <- 1
#alpha <- 1

#adJULES_results <- '/scratch/nr278/adJULES/Results/old_results/'
#imogen_results  <- '/scratch/nr278/imogen/results/Investigating_offset/'
#imogen          <- '/scratch/nr278/imogen/'

# default
#PFT_type <- c('BT','NT','C3','C4','Sh')
#p_code   <- c('nl0','alpha','f0','tlow','tupp','rd','dcdl','dq') #could change to makde more generic
#nPFT <- length(PFT_type)
#nz   <- length(p_code)

## Specific zs
#opt_z0 <- array(NA,dim=c(nPFT,nz))
#opt_z1 <- array(NA,dim=c(nPFT,nz))

#for (iPFT in 1:nPFT){
#	opt_z1[iPFT,] <- round(loadToEnv(paste0(adJULES_results,'zsamples_',PFT_type[iPFT],'_all_',alpha,'_GPP.RData'))[['z1']],3)
#	opt_z0[iPFT,] <- loadToEnv(paste0(adJULES_results,'zsamples_',PFT_type[iPFT],'_all_',alpha,'_GPP.RData'))[['z0']]
#}


#xoff <- loadToEnv(load(paste0(adJULES_results,'zsamples_BT_single14_0_1_GPP.RData')))[['xoff']]

#parameterchoice <- c(2,7,12,17,22,46,68,73) + pft -1
#xvary <- rep(FALSE,length(xoff))
#xvary[parameterchoice] <-TRUE 

#x1          <- xoff       #
#x1[xvary]   <- opt_z1[pft,] #ip       

#jinfile <- '/scratch/nr278/adJULES/ADJULES_code/jules_verify_data/BT/control/US-Ha1_1996.jin'
#scan original control file into R
control_file  <- scan(jinfile,character(0),sep="\n",comment.char="#")

new_control_file <- control_file # create new control file object


q10_leaf_line <- grep(">INIT_MISC",control_file) + 4
oldstr <- unlist(strsplit(control_file[q10_leaf_line],"  "))[1]
newstr <- as.character(x1[1])
new_control_file[q10_leaf_line] <- sub(oldstr,newstr,new_control_file[q10_leaf_line])


nl0_line <- grep(">INIT_VEG_PFT",control_file) + 40
oldstr <- unlist(strsplit(control_file[nl0_line],split="!")   ) 
newstr <- c(paste(x1[2],x1[3],x1[4],x1[5],x1[6], sep=" , "), oldstr[2])
new_control_file[nl0_line] <- paste(newstr[1],"!",oldstr[2])

alpha_line <- grep(">INIT_VEG_PFT",control_file) + 22
oldstr <- unlist(strsplit(control_file[alpha_line],split="!")   ) 
newstr <- c(paste(x1[7],x1[8],x1[9],x1[10],x1[11], sep=" , "), oldstr[2])
new_control_file[alpha_line] <- paste(newstr[1],"!",oldstr[2])

f0_line <- grep(">INIT_VEG_PFT",control_file) + 37
oldstr <- unlist(strsplit(control_file[f0_line],split="!")   ) 
newstr <- c(paste(x1[12],x1[13],x1[14],x1[15],x1[16], sep=" , "), oldstr[2])
new_control_file[f0_line] <- paste(newstr[1],"!",oldstr[2])

tlow_line <- grep(">INIT_VEG_PFT",control_file) + 46
oldstr <- unlist(strsplit(control_file[tlow_line],split="!")   ) 
newstr <- c(paste(x1[17],x1[18],x1[19],x1[20],x1[21], sep=" , "), oldstr[2])
new_control_file[tlow_line] <- paste(newstr[1],"!",oldstr[2])

tupp_line <- grep(">INIT_VEG_PFT",control_file) + 47
oldstr <- unlist(strsplit(control_file[tupp_line],split="!")   ) 
newstr <- c(paste(x1[22],x1[23],x1[24],x1[25],x1[26], sep=" , "), oldstr[2])
new_control_file[tupp_line] <- paste(newstr[1],"!",oldstr[2])


lai_line <- grep(">INIT_VEG_PFT",control_file) + 8
oldstr <- unlist(strsplit(control_file[lai_line],split="!")   ) 
newstr <- c(paste(x1[27],x1[28],x1[29],x1[30],x1[31], sep=" , "), oldstr[2])
new_control_file[lai_line] <- paste(newstr[1],"!",oldstr[2])

canht_ft_line <- grep(">INIT_VEG_PFT",control_file) + 7
oldstr <- unlist(strsplit(control_file[canht_ft_line],split="!")   ) 
newstr <- c(paste(x1[32],x1[33],x1[34],x1[35],x1[36], sep=" , "), oldstr[2])
new_control_file[canht_ft_line] <- paste(newstr[1],"!",oldstr[2])

satcon_line <- grep(">DATA_DZSOIL",control_file) + 6
oldstr <- unlist(strsplit(control_file[satcon_line],split="!")   ) 
newstr <- c(paste(x1[37],x1[38],x1[39],x1[40], sep=" , "), oldstr[2])
new_control_file[satcon_line] <- paste(newstr[1],"!",oldstr[2])

b_line <- grep(">DATA_DZSOIL",control_file) + 4
oldstr <- unlist(strsplit(control_file[b_line],split="!")   ) 
newstr <- c(paste(x1[41],x1[42],x1[43],x1[44], sep=" , "), oldstr[2])
new_control_file[b_line] <- paste(newstr[1],"!",oldstr[2])

print("check that cs works.....")
cs_line <- grep(">INIT_IC",control_file) + 36 
oldstr <- unlist(strsplit(control_file[cs_line],split="!")   ) 
newstr <- c(x1[45], oldstr[2])
new_control_file[cs_line] <- paste(newstr[1],"!",oldstr[2])
print(control_file[cs_line])
print(new_control_file[cs_line])


rootd_ft_line <- grep(">INIT_VEG_PFT",control_file) + 14
oldstr <- unlist(strsplit(control_file[rootd_ft_line],split="!")   ) 
newstr <- c(paste(x1[46],x1[47],x1[48],x1[49],x1[50], sep=" , "), oldstr[2])
new_control_file[rootd_ft_line] <- paste(newstr[1],"!",oldstr[2])

albsoil_line <- grep(">DATA_DZSOIL",control_file) + 12
oldstr <- unlist(strsplit(control_file[albsoil_line],split="!")   ) 
newstr <- c(x1[59], oldstr[2])
new_control_file[albsoil_line] <- paste(newstr[1],"!",oldstr[2])

sm_crit_line <- grep(">DATA_DZSOIL",control_file) + 8
oldstr <- unlist(strsplit(control_file[sm_crit_line],split="!")   ) 
newstr <- c(paste(x1[60]*x1[82],x1[61]*x1[83],x1[62]*x1[84],x1[63]*x1[85], sep=" , "), oldstr[2])
new_control_file[sm_crit_line] <- paste(newstr[1],"!",oldstr[2])

sm_wilt_line <- grep(">DATA_DZSOIL",control_file) + 9
oldstr <- unlist(strsplit(control_file[sm_wilt_line],split="!")   ) 
newstr <- c(paste(x1[64]*x1[60]*x1[82],x1[65]*x1[61]*x1[83],x1[66]*x1[62]*x1[84],x1[67]*x1[63]*x1[85], sep=" , "), oldstr[2])
new_control_file[sm_wilt_line] <- paste(newstr[1],"!",oldstr[2])

dcatch_dlai_line <- grep(">INIT_VEG_PFT",control_file) + 10
oldstr <- unlist(strsplit(control_file[dcatch_dlai_line],split="!")   ) 
newstr <- c(paste(x1[68],x1[69],x1[70],x1[71],x1[72], sep=" , "), oldstr[2])
new_control_file[dcatch_dlai_line] <- paste(newstr[1],"!",oldstr[2])

dqcrit_line <- grep(">INIT_VEG_PFT",control_file) + 35
oldstr <- unlist(strsplit(control_file[dqcrit_line],split="!")   ) 
newstr <- c(paste(x1[73],x1[74],x1[75],x1[76],x1[77], sep=" , "), oldstr[2])
new_control_file[dqcrit_line] <- paste(newstr[1],"!",oldstr[2])

sathh_line <- grep(">DATA_DZSOIL",control_file) + 5
oldstr <- unlist(strsplit(control_file[sathh_line],split="!")   ) 
newstr <- c(paste(x1[78],x1[79],x1[80],x1[81], sep=" , "), oldstr[2])
new_control_file[sathh_line] <- paste(newstr[1],"!",oldstr[2])

sm_sat_line <- grep(">DATA_DZSOIL",control_file) + 7
oldstr <- unlist(strsplit(control_file[sm_sat_line],split="!")   ) 
newstr <- c(paste(x1[82],x1[83],x1[84],x1[85], sep=" , "), oldstr[2])
new_control_file[sm_sat_line] <- paste(newstr[1],"!",oldstr[2])

hcap_line <- grep(">DATA_DZSOIL",control_file) + 10
oldstr <- unlist(strsplit(control_file[hcap_line],split="!")   ) 
newstr <- c(paste(x1[86],x1[87],x1[88],x1[89], sep=" , "), oldstr[2])
new_control_file[hcap_line] <- paste(newstr[1],"!",oldstr[2])

hcon_line <- grep(">DATA_DZSOIL",control_file) + 11
oldstr <- unlist(strsplit(control_file[hcon_line],split="!")   ) 
newstr <- c(paste(x1[90],x1[91],x1[92],x1[92], sep=" , "), oldstr[2])
new_control_file[hcon_line] <- paste(newstr[1],"!",oldstr[2])

q10_soil_line <- grep(">INIT_MISC",control_file) + 7
oldstr <- unlist(strsplit(control_file[q10_soil_line],split="!")   ) 
newstr <- c(x1[94], oldstr[2])
new_control_file[q10_soil_line] <- paste(newstr[1],"!",oldstr[2])


write(new_control_file,paste(verifydir,"/control/optimised_",jinname,".jin",sep=""))

