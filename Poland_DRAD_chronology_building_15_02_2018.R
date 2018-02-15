#install.packages("tgram")
# Load Packages -----------------------------------------------------------
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tgram)
# library()

rm(list = ls())
setwd("/Users/danielbalanzategui/Documents/1_Active_study/new_poland_pine_clsm_analysis/DRAD_poland")
tree_one_DRAD_poland<-read.table("104_tracheid_in.txt",header=TRUE, sep="\t")
tree_two_DRAD_poland<-read.table("106_tracheid_in.txt",header=TRUE, sep="\t")
tree_three_DRAD_poland<-read.table("107_tracheid_in.txt",header=TRUE, sep="\t")
tree_four_DRAD_poland<-read.table("113_tracheid_in.txt",header=TRUE, sep="\t")
tree_five_DRAD_poland<-read.table("baby119b_cellb_output.txt",header=TRUE, sep="\t")

#####
# setup probabilities 
#####
probs <- seq(.0, 1,length.out = 21) # percentiles to base stds on
#####
# Multiple parameters
#####
#params <- c("CWTRAD", "DRAD_poland", "DTAN")
#multi_strippers <- lapply(params, function(x){interp.standardise(dat = tree_one_DRAD_poland,probs = probs,listed = F,param = x)}) %>% 
      #setNames(., nm = params) 

#####
# interpolate and standardise function
#####
#dat<-tree_one
interp.standardise_DRAD_poland <- function(dat, probs, param = "DRAD", listed = F){
      # function inputs: 
      # dat = data frame from RAPTOR (includes ROW and POSITION)
      # probs = vectors of percentiles to be computed
      # ordered = if FALSE, function extracts quantiles based on param size distribution (not along radii)
      #           if TRUE, function extracts quantiles of cell position (along a radii),
      #           then the corresponding cell characteristics.
      # listed = if F, returns data.frame, if T, returns list for each year
      require(tidyr)
      require(tgram)
      
      
      
      # check if class = data.frame
      stopifnot(is.data.frame(dat), 
                grepl(pattern = 'ROW', 
                      colnames(dat)) %>% sum() >= 1)
      
      
      
      # ID years
      year_col <- grepl(pattern = 'YEAR',
                        colnames(dat))
      
      years <- dat[ ,year_col] %>% unique()
      n_years <- years %>% length()
      
      # ID Rows
      row_col <- grepl(pattern = 'ROW',
                       colnames(dat))
      
      pos_col <- grepl(pattern = 'POSITION',
                       colnames(dat))
      
      # Drop NA ROWS or POSITIONs
      dat <- dat %>% filter(!is.na(ROW), !is.na(POSITION))
      
      res_list <- list()
      # for each year extract values and compute param percentiles
      for(y in 1:n_years){
            #current_year <- 2011
            current_year <- years[y]
            
            # filter for current year
            year_dat <- dat %>% filter(YEAR == current_year)

            
            
            # identify max number of cell positions, compute pctls and store in list
            
            max_pos <- year_dat[ ,pos_col] %>% max(na.rm = T)


            
            # force ordering of rows and cells
            year_dat_ord <- year_dat[order(year_dat[ ,row_col],
                                           year_dat[ ,pos_col]), ]
            
            # standardise to max_rows (based on tgram::standz.all)
            
            year_std <- tgram::standz.all(traq = year_dat_ord[ ,param],
                                          series = year_dat_ord[ ,"ROW"],
                                          G = max_pos)
            # Extract standardised values, make into df and change column names
            
            year_std <- year_std[[1]] 
            
            
            
            year_pctl <- apply(year_std, MAR = 1, 
                               FUN = function(x) x[length(x) - ceiling(quantile(1:length(x), 
                                                              probs = probs)) + 1])
            
            year_pctl <- year_pctl %>% t()
            
            
            col_names <- (probs * 100) %>% as.character()
            colnames(year_pctl) <- col_names
            
            year_pctl <- year_pctl %>% as.data.frame()
            

            # bind year and row column to each list element
            
            # 
            # year_pctl$Percentile <- colnames(year_pctl) %>% sub(pattern = '%',
                                                                # replacement = '')

            
            year_pctl$YEAR <- current_year
            year_pctl$ROW <- as.numeric(rownames(year_pctl))
            
            # # 
            year_res <- gather_(data = year_pctl,
                               key_col = 'Percentiles',
                               value_col = param,
                               gather_cols = col_names)
            

            year_res$Percentiles <- year_res$Percentiles %>% sub(pattern = '%',
                                                               replacement = '') %>%
                  as.numeric()
            # year_res <- year_res[ order(year_res$R), ]

            # split resulting year DF by ROW and compute percentiles according to 
            # probs argument            
      #       split_row_dat <- year_dat %>% split(f = .[ ,row_col])
      #       
      #       if(ordered == F){
      #             
      #             
      #             row_pctls <- sapply(split_row_dat, 
      #                                 FUN = function(x) quantile(x[,'DRAD_poland'],
      #                                                            probs = probs))
      #             
      #             
      #             # bind a year column and pctl column to the data set
      #             
      #             
      #             pctls <- row.names(row_pctls) %>% 
      #                   sub(pattern = '%', 
      #                       replacement = '')
      #             
      #             pctls <- cbind(row_pctls,
      #                            data.frame('Percentiles' = as.numeric(pctls),
      #                                       'YEAR' = as.character(current_year),
      #                                       stringsAsFactors = F))
      #             
      #             
      #             
      #             
      #       } else if(ordered == T) {
      #             
      #             
      #             # If Ordered == T use cell position rather than absolute DRAD_poland size
      #             # to choose quantiles (i.e. standardize by position along radii rather than
      #             # cell size)
      #             
      #             
      #             pos_pctls <- split_row_dat %>% 
      #                   sapply(., function(x){quantile(x[,'POSITION'],
      #                                                  probs = probs)
      #                         
      #                   })
      #             
      #             pos_pctls <- pos_pctls %>% ceiling()
      #             
      #             row_pctls <- sapply(seq_along(split_row_dat), 
      #                                 FUN = function(x) {split_row_dat[[x]][pos_pctls[,x],'DRAD_poland']})
      #             
      #             
      #             # bind a year column and pctl column to the data set
      #             pctls <- row.names(pos_pctls) %>% 
      #                   sub(pattern = '%', 
      #                       replacement = '')
      #             
      #             
      #             pctls <- cbind(row_pctls,
      #                            data.frame('Percentiles' = as.numeric(pctls),
      #                                       'YEAR' = as.character(current_year),
      #                                       stringsAsFactors = F))
      #             
      #             
      #             
      #       }
      #       
      #       pctls_gathered <- tidyr::gather(data = as.data.frame(pctls), 
      #                                       key = 'ROW',
      #                                       value = 'DRAD_poland',
      #                                       -Percentiles, -YEAR)
      #       
      #       res_list[[y]] <- pctls_gathered
      #       
      #       
      # }
            
            
            
            year_res$ORIGIN <- 'Interp. Percentile Estimate'
            res_list[[y]] <- year_res[order(year_res[ ,"ROW"], year_res[ ,"Percentiles"]), ]
            # res_list[[y]] <- year_res
      }
      
      if(listed == T){
            return(res_list)
            
      } else {
            do.call(rbind, res_list) %>% return()
      }
}
# interp_stand_tracheid_tree_one <- interp.standardise(dat = tree_one,probs = probs,listed = F)
#interp_new_data <- interp.standardise(dat = tree_one,probs = probs,listed = F,param = "DRAD_poland")

interp_stand_tracheid_tree_one_DRAD_poland <- interp.standardise_DRAD_poland(dat = tree_one_DRAD_poland,probs = probs,listed = F,param = "DRAD")
interp_stand_tracheid_tree_two_DRAD_poland <- interp.standardise_DRAD_poland(dat = tree_two_DRAD_poland,probs = probs,listed = F,param = "DRAD")
interp_stand_tracheid_tree_three_DRAD_poland <- interp.standardise_DRAD_poland(dat = tree_three_DRAD_poland,probs = probs,listed = F,param = "DRAD")
interp_stand_tracheid_tree_four_DRAD_poland <- interp.standardise_DRAD_poland(dat = tree_four_DRAD_poland,probs = probs,listed = F,param = "DRAD")
interp_stand_tracheid_tree_five_DRAD_poland <- interp.standardise_DRAD_poland(dat = tree_five_DRAD_poland,probs = probs,listed = F,param = "DRAD")


#####
# tree_one
#####
colnames(interp_stand_tracheid_tree_one_DRAD_poland)<-c("Year_tree_one_DRAD_poland","Row_tree_one_DRAD_poland","Percentiles_tree_one_DRAD_poland", "tree_one_DRAD_poland","ORIGIN_tree_one_DRAD_poland")
attach(interp_stand_tracheid_tree_one_DRAD_poland)
interp_stand_tracheid_tree_one_annual_DRAD_poland <- ddply(interp_stand_tracheid_tree_one_DRAD_poland, c("Year_tree_one_DRAD_poland","Percentiles_tree_one_DRAD_poland"), summarise,tree_one_DRAD_poland=mean(tree_one_DRAD_poland))
zero_tree_one_DRAD_poland<-interp_stand_tracheid_tree_one_annual_DRAD_poland[with(interp_stand_tracheid_tree_one_annual_DRAD_poland, Percentiles_tree_one_DRAD_poland==0),]
colnames(zero_tree_one_DRAD_poland)<-c("Year","Percentiles_tree_one","zero_DRAD_poland_tree_one")
five_tree_one_DRAD_poland<-interp_stand_tracheid_tree_one_annual_DRAD_poland[with(interp_stand_tracheid_tree_one_annual_DRAD_poland, Percentiles_tree_one_DRAD_poland==5),]
colnames(five_tree_one_DRAD_poland)<-c("Year","Percentiles_tree_one","five_DRAD_poland_tree_one")
ten_tree_one_DRAD_poland<-interp_stand_tracheid_tree_one_annual_DRAD_poland[with(interp_stand_tracheid_tree_one_annual_DRAD_poland, Percentiles_tree_one_DRAD_poland==10),]
colnames(ten_tree_one_DRAD_poland)<-c("Year","Percentiles_tree_one","ten_DRAD_poland_tree_one")
fifteen_tree_one_DRAD_poland<-interp_stand_tracheid_tree_one_annual_DRAD_poland[with(interp_stand_tracheid_tree_one_annual_DRAD_poland, Percentiles_tree_one_DRAD_poland==15),]
colnames(fifteen_tree_one_DRAD_poland)<-c("Year","Percentiles_tree_one","fifteen_DRAD_poland_tree_one")
twenty_tree_one_DRAD_poland<-interp_stand_tracheid_tree_one_annual_DRAD_poland[with(interp_stand_tracheid_tree_one_annual_DRAD_poland, Percentiles_tree_one_DRAD_poland==20),]
colnames(twenty_tree_one_DRAD_poland)<-c("Year","Percentiles_tree_one","twenty_DRAD_poland_tree_one")
twentyfive_tree_one_DRAD_poland<-interp_stand_tracheid_tree_one_annual_DRAD_poland[with(interp_stand_tracheid_tree_one_annual_DRAD_poland, Percentiles_tree_one_DRAD_poland==25),]
colnames(twentyfive_tree_one_DRAD_poland)<-c("Year","Percentiles_tree_one","twentyfive_DRAD_poland_tree_one")
thirty_tree_one_DRAD_poland<-interp_stand_tracheid_tree_one_annual_DRAD_poland[with(interp_stand_tracheid_tree_one_annual_DRAD_poland, Percentiles_tree_one_DRAD_poland==30),]
colnames(thirty_tree_one_DRAD_poland)<-c("Year","Percentiles_tree_one","thirty_DRAD_poland_tree_one")
thirtyfive_tree_one_DRAD_poland<-interp_stand_tracheid_tree_one_annual_DRAD_poland[with(interp_stand_tracheid_tree_one_annual_DRAD_poland, Percentiles_tree_one_DRAD_poland==35),]
colnames(thirtyfive_tree_one_DRAD_poland)<-c("Year","Percentiles_tree_one","thirtyfive_DRAD_poland_tree_one")
fourty_tree_one_DRAD_poland<-interp_stand_tracheid_tree_one_annual_DRAD_poland[with(interp_stand_tracheid_tree_one_annual_DRAD_poland, Percentiles_tree_one_DRAD_poland==40),]
colnames(fourty_tree_one_DRAD_poland)<-c("Year","Percentiles_tree_one","fourty_DRAD_poland_tree_one")
fourtyfive_tree_one_DRAD_poland<-interp_stand_tracheid_tree_one_annual_DRAD_poland[with(interp_stand_tracheid_tree_one_annual_DRAD_poland, Percentiles_tree_one_DRAD_poland==45),]
colnames(fourtyfive_tree_one_DRAD_poland)<-c("Year","Percentiles_tree_one","fourtyfive_DRAD_poland_tree_one")
fifty_tree_one_DRAD_poland<-interp_stand_tracheid_tree_one_annual_DRAD_poland[with(interp_stand_tracheid_tree_one_annual_DRAD_poland, Percentiles_tree_one_DRAD_poland==50),]
colnames(fifty_tree_one_DRAD_poland)<-c("Year","Percentiles_tree_one","fifty_DRAD_poland_tree_one")
fiftyfive_tree_one_DRAD_poland<-interp_stand_tracheid_tree_one_annual_DRAD_poland[with(interp_stand_tracheid_tree_one_annual_DRAD_poland, Percentiles_tree_one_DRAD_poland==55),]
colnames(fiftyfive_tree_one_DRAD_poland)<-c("Year","Percentiles_tree_one","fiftyfive_DRAD_poland_tree_one")
sixty_tree_one_DRAD_poland<-interp_stand_tracheid_tree_one_annual_DRAD_poland[with(interp_stand_tracheid_tree_one_annual_DRAD_poland, Percentiles_tree_one_DRAD_poland==60),]
colnames(sixty_tree_one_DRAD_poland)<-c("Year","Percentiles_tree_one","sixty_DRAD_poland_tree_one")
sixtyfive_tree_one_DRAD_poland<-interp_stand_tracheid_tree_one_annual_DRAD_poland[with(interp_stand_tracheid_tree_one_annual_DRAD_poland, Percentiles_tree_one_DRAD_poland==65),]
colnames(sixtyfive_tree_one_DRAD_poland)<-c("Year","Percentiles_tree_one","sixtyfive_DRAD_poland_tree_one")
seventy_tree_one_DRAD_poland<-interp_stand_tracheid_tree_one_annual_DRAD_poland[with(interp_stand_tracheid_tree_one_annual_DRAD_poland, Percentiles_tree_one_DRAD_poland==70),]
colnames(seventy_tree_one_DRAD_poland)<-c("Year","Percentiles_tree_one","seventy_DRAD_poland_tree_one")
seventyfive_tree_one_DRAD_poland<-interp_stand_tracheid_tree_one_annual_DRAD_poland[with(interp_stand_tracheid_tree_one_annual_DRAD_poland, Percentiles_tree_one_DRAD_poland==75),]
colnames(seventyfive_tree_one_DRAD_poland)<-c("Year","Percentiles_tree_one","seventyfive_DRAD_poland_tree_one")
eighty_tree_one_DRAD_poland<-interp_stand_tracheid_tree_one_annual_DRAD_poland[with(interp_stand_tracheid_tree_one_annual_DRAD_poland, Percentiles_tree_one_DRAD_poland==80),]
colnames(eighty_tree_one_DRAD_poland)<-c("Year","Percentiles_tree_one","eighty_DRAD_poland_tree_one")
eightyfive_tree_one_DRAD_poland<-interp_stand_tracheid_tree_one_annual_DRAD_poland[with(interp_stand_tracheid_tree_one_annual_DRAD_poland, Percentiles_tree_one_DRAD_poland==85),]
colnames(eightyfive_tree_one_DRAD_poland)<-c("Year","Percentiles_tree_one","eightyfive_DRAD_poland_tree_one")
ninety_tree_one_DRAD_poland<-interp_stand_tracheid_tree_one_annual_DRAD_poland[with(interp_stand_tracheid_tree_one_annual_DRAD_poland, Percentiles_tree_one_DRAD_poland==90),]
colnames(ninety_tree_one_DRAD_poland)<-c("Year","Percentiles_tree_one","ninety_DRAD_poland_tree_one")
ninetyfive_tree_one_DRAD_poland<-interp_stand_tracheid_tree_one_annual_DRAD_poland[with(interp_stand_tracheid_tree_one_annual_DRAD_poland, Percentiles_tree_one_DRAD_poland==95),]
colnames(ninetyfive_tree_one_DRAD_poland)<-c("Year","Percentiles_tree_one","ninetyfive_DRAD_poland_tree_one")
hundred_tree_one_DRAD_poland<-interp_stand_tracheid_tree_one_annual_DRAD_poland[with(interp_stand_tracheid_tree_one_annual_DRAD_poland, Percentiles_tree_one_DRAD_poland==100),]
colnames(hundred_tree_one_DRAD_poland)<-c("Year","Percentiles_tree_one","hundred_DRAD_poland_tree_one")

#####
# tree_two_DRAD_poland
#####
colnames(interp_stand_tracheid_tree_two_DRAD_poland)<-c("Year_tree_two_DRAD_poland","Row_tree_two_DRAD_poland","Percentiles_tree_two_DRAD_poland", "tree_two_DRAD_poland","ORIGIN_tree_two_DRAD_poland")
attach(interp_stand_tracheid_tree_two_DRAD_poland)
interp_stand_tracheid_tree_two_annual_DRAD_poland <- ddply(interp_stand_tracheid_tree_two_DRAD_poland, c("Year_tree_two_DRAD_poland","Percentiles_tree_two_DRAD_poland"), summarise,tree_two_DRAD_poland=mean(tree_two_DRAD_poland))
zero_tree_two_DRAD_poland<-interp_stand_tracheid_tree_two_annual_DRAD_poland[with(interp_stand_tracheid_tree_two_annual_DRAD_poland, Percentiles_tree_two_DRAD_poland==0),]
colnames(zero_tree_two_DRAD_poland)<-c("Year","Percentiles_tree_two","zero_DRAD_poland_tree_two")
five_tree_two_DRAD_poland<-interp_stand_tracheid_tree_two_annual_DRAD_poland[with(interp_stand_tracheid_tree_two_annual_DRAD_poland, Percentiles_tree_two_DRAD_poland==5),]
colnames(five_tree_two_DRAD_poland)<-c("Year","Percentiles_tree_two","five_DRAD_poland_tree_two")
ten_tree_two_DRAD_poland<-interp_stand_tracheid_tree_two_annual_DRAD_poland[with(interp_stand_tracheid_tree_two_annual_DRAD_poland, Percentiles_tree_two_DRAD_poland==10),]
colnames(ten_tree_two_DRAD_poland)<-c("Year","Percentiles_tree_two","ten_DRAD_poland_tree_two")
fifteen_tree_two_DRAD_poland<-interp_stand_tracheid_tree_two_annual_DRAD_poland[with(interp_stand_tracheid_tree_two_annual_DRAD_poland, Percentiles_tree_two_DRAD_poland==15),]
colnames(fifteen_tree_two_DRAD_poland)<-c("Year","Percentiles_tree_two","fifteen_DRAD_poland_tree_two")
twenty_tree_two_DRAD_poland<-interp_stand_tracheid_tree_two_annual_DRAD_poland[with(interp_stand_tracheid_tree_two_annual_DRAD_poland, Percentiles_tree_two_DRAD_poland==20),]
colnames(twenty_tree_two_DRAD_poland)<-c("Year","Percentiles_tree_two","twenty_DRAD_poland_tree_two")
twentyfive_tree_two_DRAD_poland<-interp_stand_tracheid_tree_two_annual_DRAD_poland[with(interp_stand_tracheid_tree_two_annual_DRAD_poland, Percentiles_tree_two_DRAD_poland==25),]
colnames(twentyfive_tree_two_DRAD_poland)<-c("Year","Percentiles_tree_two","twentyfive_DRAD_poland_tree_two")
thirty_tree_two_DRAD_poland<-interp_stand_tracheid_tree_two_annual_DRAD_poland[with(interp_stand_tracheid_tree_two_annual_DRAD_poland, Percentiles_tree_two_DRAD_poland==30),]
colnames(thirty_tree_two_DRAD_poland)<-c("Year","Percentiles_tree_two","thirty_DRAD_poland_tree_two")
thirtyfive_tree_two_DRAD_poland<-interp_stand_tracheid_tree_two_annual_DRAD_poland[with(interp_stand_tracheid_tree_two_annual_DRAD_poland, Percentiles_tree_two_DRAD_poland==35),]
colnames(thirtyfive_tree_two_DRAD_poland)<-c("Year","Percentiles_tree_two","thirtyfive_DRAD_poland_tree_two")
fourty_tree_two_DRAD_poland<-interp_stand_tracheid_tree_two_annual_DRAD_poland[with(interp_stand_tracheid_tree_two_annual_DRAD_poland, Percentiles_tree_two_DRAD_poland==40),]
colnames(fourty_tree_two_DRAD_poland)<-c("Year","Percentiles_tree_two","fourty_DRAD_poland_tree_two")
fourtyfive_tree_two_DRAD_poland<-interp_stand_tracheid_tree_two_annual_DRAD_poland[with(interp_stand_tracheid_tree_two_annual_DRAD_poland, Percentiles_tree_two_DRAD_poland==45),]
colnames(fourtyfive_tree_two_DRAD_poland)<-c("Year","Percentiles_tree_two","fourtyfive_DRAD_poland_tree_two")
fifty_tree_two_DRAD_poland<-interp_stand_tracheid_tree_two_annual_DRAD_poland[with(interp_stand_tracheid_tree_two_annual_DRAD_poland, Percentiles_tree_two_DRAD_poland==50),]
colnames(fifty_tree_two_DRAD_poland)<-c("Year","Percentiles_tree_two","fifty_DRAD_poland_tree_two")
fiftyfive_tree_two_DRAD_poland<-interp_stand_tracheid_tree_two_annual_DRAD_poland[with(interp_stand_tracheid_tree_two_annual_DRAD_poland, Percentiles_tree_two_DRAD_poland==55),]
colnames(fiftyfive_tree_two_DRAD_poland)<-c("Year","Percentiles_tree_two","fiftyfive_DRAD_poland_tree_two")
sixty_tree_two_DRAD_poland<-interp_stand_tracheid_tree_two_annual_DRAD_poland[with(interp_stand_tracheid_tree_two_annual_DRAD_poland, Percentiles_tree_two_DRAD_poland==60),]
colnames(sixty_tree_two_DRAD_poland)<-c("Year","Percentiles_tree_two","sixty_DRAD_poland_tree_two")
sixtyfive_tree_two_DRAD_poland<-interp_stand_tracheid_tree_two_annual_DRAD_poland[with(interp_stand_tracheid_tree_two_annual_DRAD_poland, Percentiles_tree_two_DRAD_poland==65),]
colnames(sixtyfive_tree_two_DRAD_poland)<-c("Year","Percentiles_tree_two","sixtyfive_DRAD_poland_tree_two")
seventy_tree_two_DRAD_poland<-interp_stand_tracheid_tree_two_annual_DRAD_poland[with(interp_stand_tracheid_tree_two_annual_DRAD_poland, Percentiles_tree_two_DRAD_poland==70),]
colnames(seventy_tree_two_DRAD_poland)<-c("Year","Percentiles_tree_two","seventy_DRAD_poland_tree_two")
seventyfive_tree_two_DRAD_poland<-interp_stand_tracheid_tree_two_annual_DRAD_poland[with(interp_stand_tracheid_tree_two_annual_DRAD_poland, Percentiles_tree_two_DRAD_poland==75),]
colnames(seventyfive_tree_two_DRAD_poland)<-c("Year","Percentiles_tree_two","seventyfive_DRAD_poland_tree_two")
eighty_tree_two_DRAD_poland<-interp_stand_tracheid_tree_two_annual_DRAD_poland[with(interp_stand_tracheid_tree_two_annual_DRAD_poland, Percentiles_tree_two_DRAD_poland==80),]
colnames(eighty_tree_two_DRAD_poland)<-c("Year","Percentiles_tree_two","eighty_DRAD_poland_tree_two")
eightyfive_tree_two_DRAD_poland<-interp_stand_tracheid_tree_two_annual_DRAD_poland[with(interp_stand_tracheid_tree_two_annual_DRAD_poland, Percentiles_tree_two_DRAD_poland==85),]
colnames(eightyfive_tree_two_DRAD_poland)<-c("Year","Percentiles_tree_two","eightyfive_DRAD_poland_tree_two")
ninety_tree_two_DRAD_poland<-interp_stand_tracheid_tree_two_annual_DRAD_poland[with(interp_stand_tracheid_tree_two_annual_DRAD_poland, Percentiles_tree_two_DRAD_poland==90),]
colnames(ninety_tree_two_DRAD_poland)<-c("Year","Percentiles_tree_two","ninety_DRAD_poland_tree_two")
ninetyfive_tree_two_DRAD_poland<-interp_stand_tracheid_tree_two_annual_DRAD_poland[with(interp_stand_tracheid_tree_two_annual_DRAD_poland, Percentiles_tree_two_DRAD_poland==95),]
colnames(ninetyfive_tree_two_DRAD_poland)<-c("Year","Percentiles_tree_two","ninetyfive_DRAD_poland_tree_two")
hundred_tree_two_DRAD_poland<-interp_stand_tracheid_tree_two_annual_DRAD_poland[with(interp_stand_tracheid_tree_two_annual_DRAD_poland, Percentiles_tree_two_DRAD_poland==100),]
colnames(hundred_tree_two_DRAD_poland)<-c("Year","Percentiles_tree_two","hundred_DRAD_poland_tree_two")

#####
# tree_three_DRAD_poland
#####
colnames(interp_stand_tracheid_tree_three_DRAD_poland)<-c("Year_tree_three_DRAD_poland","Row_tree_three_DRAD_poland","Percentiles_tree_three_DRAD_poland", "tree_three_DRAD_poland","ORIGIN_tree_three_DRAD_poland")
attach(interp_stand_tracheid_tree_three_DRAD_poland)
interp_stand_tracheid_tree_three_annual_DRAD_poland <- ddply(interp_stand_tracheid_tree_three_DRAD_poland, c("Year_tree_three_DRAD_poland","Percentiles_tree_three_DRAD_poland"), summarise,tree_three_DRAD_poland=mean(tree_three_DRAD_poland))
zero_tree_three_DRAD_poland<-interp_stand_tracheid_tree_three_annual_DRAD_poland[with(interp_stand_tracheid_tree_three_annual_DRAD_poland, Percentiles_tree_three_DRAD_poland==0),]
colnames(zero_tree_three_DRAD_poland)<-c("Year","Percentiles_tree_three","zero_DRAD_poland_tree_three")
five_tree_three_DRAD_poland<-interp_stand_tracheid_tree_three_annual_DRAD_poland[with(interp_stand_tracheid_tree_three_annual_DRAD_poland, Percentiles_tree_three_DRAD_poland==5),]
colnames(five_tree_three_DRAD_poland)<-c("Year","Percentiles_tree_three","five_DRAD_poland_tree_three")
ten_tree_three_DRAD_poland<-interp_stand_tracheid_tree_three_annual_DRAD_poland[with(interp_stand_tracheid_tree_three_annual_DRAD_poland, Percentiles_tree_three_DRAD_poland==10),]
colnames(ten_tree_three_DRAD_poland)<-c("Year","Percentiles_tree_three","ten_DRAD_poland_tree_three")
fifteen_tree_three_DRAD_poland<-interp_stand_tracheid_tree_three_annual_DRAD_poland[with(interp_stand_tracheid_tree_three_annual_DRAD_poland, Percentiles_tree_three_DRAD_poland==15),]
colnames(fifteen_tree_three_DRAD_poland)<-c("Year","Percentiles_tree_three","fifteen_DRAD_poland_tree_three")
twenty_tree_three_DRAD_poland<-interp_stand_tracheid_tree_three_annual_DRAD_poland[with(interp_stand_tracheid_tree_three_annual_DRAD_poland, Percentiles_tree_three_DRAD_poland==20),]
colnames(twenty_tree_three_DRAD_poland)<-c("Year","Percentiles_tree_three","twenty_DRAD_poland_tree_three")
twentyfive_tree_three_DRAD_poland<-interp_stand_tracheid_tree_three_annual_DRAD_poland[with(interp_stand_tracheid_tree_three_annual_DRAD_poland, Percentiles_tree_three_DRAD_poland==25),]
colnames(twentyfive_tree_three_DRAD_poland)<-c("Year","Percentiles_tree_three","twentyfive_DRAD_poland_tree_three")
thirty_tree_three_DRAD_poland<-interp_stand_tracheid_tree_three_annual_DRAD_poland[with(interp_stand_tracheid_tree_three_annual_DRAD_poland, Percentiles_tree_three_DRAD_poland==30),]
colnames(thirty_tree_three_DRAD_poland)<-c("Year","Percentiles_tree_three","thirty_DRAD_poland_tree_three")
thirtyfive_tree_three_DRAD_poland<-interp_stand_tracheid_tree_three_annual_DRAD_poland[with(interp_stand_tracheid_tree_three_annual_DRAD_poland, Percentiles_tree_three_DRAD_poland==35),]
colnames(thirtyfive_tree_three_DRAD_poland)<-c("Year","Percentiles_tree_three","thirtyfive_DRAD_poland_tree_three")
fourty_tree_three_DRAD_poland<-interp_stand_tracheid_tree_three_annual_DRAD_poland[with(interp_stand_tracheid_tree_three_annual_DRAD_poland, Percentiles_tree_three_DRAD_poland==40),]
colnames(fourty_tree_three_DRAD_poland)<-c("Year","Percentiles_tree_three","fourty_DRAD_poland_tree_three")
fourtyfive_tree_three_DRAD_poland<-interp_stand_tracheid_tree_three_annual_DRAD_poland[with(interp_stand_tracheid_tree_three_annual_DRAD_poland, Percentiles_tree_three_DRAD_poland==45),]
colnames(fourtyfive_tree_three_DRAD_poland)<-c("Year","Percentiles_tree_three","fourtyfive_DRAD_poland_tree_three")
fifty_tree_three_DRAD_poland<-interp_stand_tracheid_tree_three_annual_DRAD_poland[with(interp_stand_tracheid_tree_three_annual_DRAD_poland, Percentiles_tree_three_DRAD_poland==50),]
colnames(fifty_tree_three_DRAD_poland)<-c("Year","Percentiles_tree_three","fifty_DRAD_poland_tree_three")
fiftyfive_tree_three_DRAD_poland<-interp_stand_tracheid_tree_three_annual_DRAD_poland[with(interp_stand_tracheid_tree_three_annual_DRAD_poland, Percentiles_tree_three_DRAD_poland==55),]
colnames(fiftyfive_tree_three_DRAD_poland)<-c("Year","Percentiles_tree_three","fiftyfive_DRAD_poland_tree_three")
sixty_tree_three_DRAD_poland<-interp_stand_tracheid_tree_three_annual_DRAD_poland[with(interp_stand_tracheid_tree_three_annual_DRAD_poland, Percentiles_tree_three_DRAD_poland==60),]
colnames(sixty_tree_three_DRAD_poland)<-c("Year","Percentiles_tree_three","sixty_DRAD_poland_tree_three")
sixtyfive_tree_three_DRAD_poland<-interp_stand_tracheid_tree_three_annual_DRAD_poland[with(interp_stand_tracheid_tree_three_annual_DRAD_poland, Percentiles_tree_three_DRAD_poland==65),]
colnames(sixtyfive_tree_three_DRAD_poland)<-c("Year","Percentiles_tree_three","sixtyfive_DRAD_poland_tree_three")
seventy_tree_three_DRAD_poland<-interp_stand_tracheid_tree_three_annual_DRAD_poland[with(interp_stand_tracheid_tree_three_annual_DRAD_poland, Percentiles_tree_three_DRAD_poland==70),]
colnames(seventy_tree_three_DRAD_poland)<-c("Year","Percentiles_tree_three","seventy_DRAD_poland_tree_three")
seventyfive_tree_three_DRAD_poland<-interp_stand_tracheid_tree_three_annual_DRAD_poland[with(interp_stand_tracheid_tree_three_annual_DRAD_poland, Percentiles_tree_three_DRAD_poland==75),]
colnames(seventyfive_tree_three_DRAD_poland)<-c("Year","Percentiles_tree_three","seventyfive_DRAD_poland_tree_three")
eighty_tree_three_DRAD_poland<-interp_stand_tracheid_tree_three_annual_DRAD_poland[with(interp_stand_tracheid_tree_three_annual_DRAD_poland, Percentiles_tree_three_DRAD_poland==80),]
colnames(eighty_tree_three_DRAD_poland)<-c("Year","Percentiles_tree_three","eighty_DRAD_poland_tree_three")
eightyfive_tree_three_DRAD_poland<-interp_stand_tracheid_tree_three_annual_DRAD_poland[with(interp_stand_tracheid_tree_three_annual_DRAD_poland, Percentiles_tree_three_DRAD_poland==85),]
colnames(eightyfive_tree_three_DRAD_poland)<-c("Year","Percentiles_tree_three","eightyfive_DRAD_poland_tree_three")
ninety_tree_three_DRAD_poland<-interp_stand_tracheid_tree_three_annual_DRAD_poland[with(interp_stand_tracheid_tree_three_annual_DRAD_poland, Percentiles_tree_three_DRAD_poland==90),]
colnames(ninety_tree_three_DRAD_poland)<-c("Year","Percentiles_tree_three","ninety_DRAD_poland_tree_three")
ninetyfive_tree_three_DRAD_poland<-interp_stand_tracheid_tree_three_annual_DRAD_poland[with(interp_stand_tracheid_tree_three_annual_DRAD_poland, Percentiles_tree_three_DRAD_poland==95),]
colnames(ninetyfive_tree_three_DRAD_poland)<-c("Year","Percentiles_tree_three","ninetyfive_DRAD_poland_tree_three")
hundred_tree_three_DRAD_poland<-interp_stand_tracheid_tree_three_annual_DRAD_poland[with(interp_stand_tracheid_tree_three_annual_DRAD_poland, Percentiles_tree_three_DRAD_poland==100),]
colnames(hundred_tree_three_DRAD_poland)<-c("Year","Percentiles_tree_three","hundred_DRAD_poland_tree_three")

#####
# tree_four_DRAD_poland
#####
colnames(interp_stand_tracheid_tree_four_DRAD_poland)<-c("Year_tree_four_DRAD_poland","Row_tree_four_DRAD_poland","Percentiles_tree_four_DRAD_poland", "tree_four_DRAD_poland","ORIGIN_tree_four_DRAD_poland")
attach(interp_stand_tracheid_tree_four_DRAD_poland)
interp_stand_tracheid_tree_four_annual_DRAD_poland <- ddply(interp_stand_tracheid_tree_four_DRAD_poland, c("Year_tree_four_DRAD_poland","Percentiles_tree_four_DRAD_poland"), summarise,tree_four_DRAD_poland=mean(tree_four_DRAD_poland))
zero_tree_four_DRAD_poland<-interp_stand_tracheid_tree_four_annual_DRAD_poland[with(interp_stand_tracheid_tree_four_annual_DRAD_poland, Percentiles_tree_four_DRAD_poland==0),]
colnames(zero_tree_four_DRAD_poland)<-c("Year","Percentiles_tree_four","zero_DRAD_poland_tree_four")
five_tree_four_DRAD_poland<-interp_stand_tracheid_tree_four_annual_DRAD_poland[with(interp_stand_tracheid_tree_four_annual_DRAD_poland, Percentiles_tree_four_DRAD_poland==5),]
colnames(five_tree_four_DRAD_poland)<-c("Year","Percentiles_tree_four","five_DRAD_poland_tree_four")
ten_tree_four_DRAD_poland<-interp_stand_tracheid_tree_four_annual_DRAD_poland[with(interp_stand_tracheid_tree_four_annual_DRAD_poland, Percentiles_tree_four_DRAD_poland==10),]
colnames(ten_tree_four_DRAD_poland)<-c("Year","Percentiles_tree_four","ten_DRAD_poland_tree_four")
fifteen_tree_four_DRAD_poland<-interp_stand_tracheid_tree_four_annual_DRAD_poland[with(interp_stand_tracheid_tree_four_annual_DRAD_poland, Percentiles_tree_four_DRAD_poland==15),]
colnames(fifteen_tree_four_DRAD_poland)<-c("Year","Percentiles_tree_four","fifteen_DRAD_poland_tree_four")
twenty_tree_four_DRAD_poland<-interp_stand_tracheid_tree_four_annual_DRAD_poland[with(interp_stand_tracheid_tree_four_annual_DRAD_poland, Percentiles_tree_four_DRAD_poland==20),]
colnames(twenty_tree_four_DRAD_poland)<-c("Year","Percentiles_tree_four","twenty_DRAD_poland_tree_four")
twentyfive_tree_four_DRAD_poland<-interp_stand_tracheid_tree_four_annual_DRAD_poland[with(interp_stand_tracheid_tree_four_annual_DRAD_poland, Percentiles_tree_four_DRAD_poland==25),]
colnames(twentyfive_tree_four_DRAD_poland)<-c("Year","Percentiles_tree_four","twentyfive_DRAD_poland_tree_four")
thirty_tree_four_DRAD_poland<-interp_stand_tracheid_tree_four_annual_DRAD_poland[with(interp_stand_tracheid_tree_four_annual_DRAD_poland, Percentiles_tree_four_DRAD_poland==30),]
colnames(thirty_tree_four_DRAD_poland)<-c("Year","Percentiles_tree_four","thirty_DRAD_poland_tree_four")
thirtyfive_tree_four_DRAD_poland<-interp_stand_tracheid_tree_four_annual_DRAD_poland[with(interp_stand_tracheid_tree_four_annual_DRAD_poland, Percentiles_tree_four_DRAD_poland==35),]
colnames(thirtyfive_tree_four_DRAD_poland)<-c("Year","Percentiles_tree_four","thirtyfive_DRAD_poland_tree_four")
fourty_tree_four_DRAD_poland<-interp_stand_tracheid_tree_four_annual_DRAD_poland[with(interp_stand_tracheid_tree_four_annual_DRAD_poland, Percentiles_tree_four_DRAD_poland==40),]
colnames(fourty_tree_four_DRAD_poland)<-c("Year","Percentiles_tree_four","fourty_DRAD_poland_tree_four")
fourtyfive_tree_four_DRAD_poland<-interp_stand_tracheid_tree_four_annual_DRAD_poland[with(interp_stand_tracheid_tree_four_annual_DRAD_poland, Percentiles_tree_four_DRAD_poland==45),]
colnames(fourtyfive_tree_four_DRAD_poland)<-c("Year","Percentiles_tree_four","fourtyfive_DRAD_poland_tree_four")
fifty_tree_four_DRAD_poland<-interp_stand_tracheid_tree_four_annual_DRAD_poland[with(interp_stand_tracheid_tree_four_annual_DRAD_poland, Percentiles_tree_four_DRAD_poland==50),]
colnames(fifty_tree_four_DRAD_poland)<-c("Year","Percentiles_tree_four","fifty_DRAD_poland_tree_four")
fiftyfive_tree_four_DRAD_poland<-interp_stand_tracheid_tree_four_annual_DRAD_poland[with(interp_stand_tracheid_tree_four_annual_DRAD_poland, Percentiles_tree_four_DRAD_poland==55),]
colnames(fiftyfive_tree_four_DRAD_poland)<-c("Year","Percentiles_tree_four","fiftyfive_DRAD_poland_tree_four")
sixty_tree_four_DRAD_poland<-interp_stand_tracheid_tree_four_annual_DRAD_poland[with(interp_stand_tracheid_tree_four_annual_DRAD_poland, Percentiles_tree_four_DRAD_poland==60),]
colnames(sixty_tree_four_DRAD_poland)<-c("Year","Percentiles_tree_four","sixty_DRAD_poland_tree_four")
sixtyfive_tree_four_DRAD_poland<-interp_stand_tracheid_tree_four_annual_DRAD_poland[with(interp_stand_tracheid_tree_four_annual_DRAD_poland, Percentiles_tree_four_DRAD_poland==65),]
colnames(sixtyfive_tree_four_DRAD_poland)<-c("Year","Percentiles_tree_four","sixtyfive_DRAD_poland_tree_four")
seventy_tree_four_DRAD_poland<-interp_stand_tracheid_tree_four_annual_DRAD_poland[with(interp_stand_tracheid_tree_four_annual_DRAD_poland, Percentiles_tree_four_DRAD_poland==70),]
colnames(seventy_tree_four_DRAD_poland)<-c("Year","Percentiles_tree_four","seventy_DRAD_poland_tree_four")
seventyfive_tree_four_DRAD_poland<-interp_stand_tracheid_tree_four_annual_DRAD_poland[with(interp_stand_tracheid_tree_four_annual_DRAD_poland, Percentiles_tree_four_DRAD_poland==75),]
colnames(seventyfive_tree_four_DRAD_poland)<-c("Year","Percentiles_tree_four","seventyfive_DRAD_poland_tree_four")
eighty_tree_four_DRAD_poland<-interp_stand_tracheid_tree_four_annual_DRAD_poland[with(interp_stand_tracheid_tree_four_annual_DRAD_poland, Percentiles_tree_four_DRAD_poland==80),]
colnames(eighty_tree_four_DRAD_poland)<-c("Year","Percentiles_tree_four","eighty_DRAD_poland_tree_four")
eightyfive_tree_four_DRAD_poland<-interp_stand_tracheid_tree_four_annual_DRAD_poland[with(interp_stand_tracheid_tree_four_annual_DRAD_poland, Percentiles_tree_four_DRAD_poland==85),]
colnames(eightyfive_tree_four_DRAD_poland)<-c("Year","Percentiles_tree_four","eightyfive_DRAD_poland_tree_four")
ninety_tree_four_DRAD_poland<-interp_stand_tracheid_tree_four_annual_DRAD_poland[with(interp_stand_tracheid_tree_four_annual_DRAD_poland, Percentiles_tree_four_DRAD_poland==90),]
colnames(ninety_tree_four_DRAD_poland)<-c("Year","Percentiles_tree_four","ninety_DRAD_poland_tree_four")
ninetyfive_tree_four_DRAD_poland<-interp_stand_tracheid_tree_four_annual_DRAD_poland[with(interp_stand_tracheid_tree_four_annual_DRAD_poland, Percentiles_tree_four_DRAD_poland==95),]
colnames(ninetyfive_tree_four_DRAD_poland)<-c("Year","Percentiles_tree_four","ninetyfive_DRAD_poland_tree_four")
hundred_tree_four_DRAD_poland<-interp_stand_tracheid_tree_four_annual_DRAD_poland[with(interp_stand_tracheid_tree_four_annual_DRAD_poland, Percentiles_tree_four_DRAD_poland==100),]
colnames(hundred_tree_four_DRAD_poland)<-c("Year","Percentiles_tree_four","hundred_DRAD_poland_tree_four")

#####
# tree_five_DRAD_poland
#####
colnames(interp_stand_tracheid_tree_five_DRAD_poland)<-c("Year_tree_five_DRAD_poland","Row_tree_five_DRAD_poland","Percentiles_tree_five_DRAD_poland", "tree_five_DRAD_poland","ORIGIN_tree_five_DRAD_poland")
attach(interp_stand_tracheid_tree_five_DRAD_poland)
interp_stand_tracheid_tree_five_annual_DRAD_poland <- ddply(interp_stand_tracheid_tree_five_DRAD_poland, c("Year_tree_five_DRAD_poland","Percentiles_tree_five_DRAD_poland"), summarise,tree_five_DRAD_poland=mean(tree_five_DRAD_poland))
zero_tree_five_DRAD_poland<-interp_stand_tracheid_tree_five_annual_DRAD_poland[with(interp_stand_tracheid_tree_five_annual_DRAD_poland, Percentiles_tree_five_DRAD_poland==0),]
colnames(zero_tree_five_DRAD_poland)<-c("Year","Percentiles_tree_five","zero_DRAD_poland_tree_five")
five_tree_five_DRAD_poland<-interp_stand_tracheid_tree_five_annual_DRAD_poland[with(interp_stand_tracheid_tree_five_annual_DRAD_poland, Percentiles_tree_five_DRAD_poland==5),]
colnames(five_tree_five_DRAD_poland)<-c("Year","Percentiles_tree_five","five_DRAD_poland_tree_five")
ten_tree_five_DRAD_poland<-interp_stand_tracheid_tree_five_annual_DRAD_poland[with(interp_stand_tracheid_tree_five_annual_DRAD_poland, Percentiles_tree_five_DRAD_poland==10),]
colnames(ten_tree_five_DRAD_poland)<-c("Year","Percentiles_tree_five","ten_DRAD_poland_tree_five")
fifteen_tree_five_DRAD_poland<-interp_stand_tracheid_tree_five_annual_DRAD_poland[with(interp_stand_tracheid_tree_five_annual_DRAD_poland, Percentiles_tree_five_DRAD_poland==15),]
colnames(fifteen_tree_five_DRAD_poland)<-c("Year","Percentiles_tree_five","fifteen_DRAD_poland_tree_five")
twenty_tree_five_DRAD_poland<-interp_stand_tracheid_tree_five_annual_DRAD_poland[with(interp_stand_tracheid_tree_five_annual_DRAD_poland, Percentiles_tree_five_DRAD_poland==20),]
colnames(twenty_tree_five_DRAD_poland)<-c("Year","Percentiles_tree_five","twenty_DRAD_poland_tree_five")
twentyfive_tree_five_DRAD_poland<-interp_stand_tracheid_tree_five_annual_DRAD_poland[with(interp_stand_tracheid_tree_five_annual_DRAD_poland, Percentiles_tree_five_DRAD_poland==25),]
colnames(twentyfive_tree_five_DRAD_poland)<-c("Year","Percentiles_tree_five","twentyfive_DRAD_poland_tree_five")
thirty_tree_five_DRAD_poland<-interp_stand_tracheid_tree_five_annual_DRAD_poland[with(interp_stand_tracheid_tree_five_annual_DRAD_poland, Percentiles_tree_five_DRAD_poland==30),]
colnames(thirty_tree_five_DRAD_poland)<-c("Year","Percentiles_tree_five","thirty_DRAD_poland_tree_five")
thirtyfive_tree_five_DRAD_poland<-interp_stand_tracheid_tree_five_annual_DRAD_poland[with(interp_stand_tracheid_tree_five_annual_DRAD_poland, Percentiles_tree_five_DRAD_poland==35),]
colnames(thirtyfive_tree_five_DRAD_poland)<-c("Year","Percentiles_tree_five","thirtyfive_DRAD_poland_tree_five")
fourty_tree_five_DRAD_poland<-interp_stand_tracheid_tree_five_annual_DRAD_poland[with(interp_stand_tracheid_tree_five_annual_DRAD_poland, Percentiles_tree_five_DRAD_poland==40),]
colnames(fourty_tree_five_DRAD_poland)<-c("Year","Percentiles_tree_five","fivety_DRAD_poland_tree_five")
fourtyfive_tree_five_DRAD_poland<-interp_stand_tracheid_tree_five_annual_DRAD_poland[with(interp_stand_tracheid_tree_five_annual_DRAD_poland, Percentiles_tree_five_DRAD_poland==45),]
colnames(fourtyfive_tree_five_DRAD_poland)<-c("Year","Percentiles_tree_five","fivetyfive_DRAD_poland_tree_five")
fifty_tree_five_DRAD_poland<-interp_stand_tracheid_tree_five_annual_DRAD_poland[with(interp_stand_tracheid_tree_five_annual_DRAD_poland, Percentiles_tree_five_DRAD_poland==50),]
colnames(fifty_tree_five_DRAD_poland)<-c("Year","Percentiles_tree_five","fifty_DRAD_poland_tree_five")
fiftyfive_tree_five_DRAD_poland<-interp_stand_tracheid_tree_five_annual_DRAD_poland[with(interp_stand_tracheid_tree_five_annual_DRAD_poland, Percentiles_tree_five_DRAD_poland==55),]
colnames(fiftyfive_tree_five_DRAD_poland)<-c("Year","Percentiles_tree_five","fiftyfive_DRAD_poland_tree_five")
sixty_tree_five_DRAD_poland<-interp_stand_tracheid_tree_five_annual_DRAD_poland[with(interp_stand_tracheid_tree_five_annual_DRAD_poland, Percentiles_tree_five_DRAD_poland==60),]
colnames(sixty_tree_five_DRAD_poland)<-c("Year","Percentiles_tree_five","sixty_DRAD_poland_tree_five")
sixtyfive_tree_five_DRAD_poland<-interp_stand_tracheid_tree_five_annual_DRAD_poland[with(interp_stand_tracheid_tree_five_annual_DRAD_poland, Percentiles_tree_five_DRAD_poland==65),]
colnames(sixtyfive_tree_five_DRAD_poland)<-c("Year","Percentiles_tree_five","sixtyfive_DRAD_poland_tree_five")
seventy_tree_five_DRAD_poland<-interp_stand_tracheid_tree_five_annual_DRAD_poland[with(interp_stand_tracheid_tree_five_annual_DRAD_poland, Percentiles_tree_five_DRAD_poland==70),]
colnames(seventy_tree_five_DRAD_poland)<-c("Year","Percentiles_tree_five","seventy_DRAD_poland_tree_five")
seventyfive_tree_five_DRAD_poland<-interp_stand_tracheid_tree_five_annual_DRAD_poland[with(interp_stand_tracheid_tree_five_annual_DRAD_poland, Percentiles_tree_five_DRAD_poland==75),]
colnames(seventyfive_tree_five_DRAD_poland)<-c("Year","Percentiles_tree_five","seventyfive_DRAD_poland_tree_five")
eighty_tree_five_DRAD_poland<-interp_stand_tracheid_tree_five_annual_DRAD_poland[with(interp_stand_tracheid_tree_five_annual_DRAD_poland, Percentiles_tree_five_DRAD_poland==80),]
colnames(eighty_tree_five_DRAD_poland)<-c("Year","Percentiles_tree_five","eighty_DRAD_poland_tree_five")
eightyfive_tree_five_DRAD_poland<-interp_stand_tracheid_tree_five_annual_DRAD_poland[with(interp_stand_tracheid_tree_five_annual_DRAD_poland, Percentiles_tree_five_DRAD_poland==85),]
colnames(eightyfive_tree_five_DRAD_poland)<-c("Year","Percentiles_tree_five","eightyfive_DRAD_poland_tree_five")
ninety_tree_five_DRAD_poland<-interp_stand_tracheid_tree_five_annual_DRAD_poland[with(interp_stand_tracheid_tree_five_annual_DRAD_poland, Percentiles_tree_five_DRAD_poland==90),]
colnames(ninety_tree_five_DRAD_poland)<-c("Year","Percentiles_tree_five","ninety_DRAD_poland_tree_five")
ninetyfive_tree_five_DRAD_poland<-interp_stand_tracheid_tree_five_annual_DRAD_poland[with(interp_stand_tracheid_tree_five_annual_DRAD_poland, Percentiles_tree_five_DRAD_poland==95),]
colnames(ninetyfive_tree_five_DRAD_poland)<-c("Year","Percentiles_tree_five","ninetyfive_DRAD_poland_tree_five")
hundred_tree_five_DRAD_poland<-interp_stand_tracheid_tree_five_annual_DRAD_poland[with(interp_stand_tracheid_tree_five_annual_DRAD_poland, Percentiles_tree_five_DRAD_poland==100),]
colnames(hundred_tree_five_DRAD_poland)<-c("Year","Percentiles_tree_five","hundred_DRAD_poland_tree_five")


#####
# zero percentile, latewood section
#####
#remove precentile column
zero_tree_one_DRAD_poland_2<-as.data.frame(zero_tree_one_DRAD_poland[,-2])
zero_tree_two_DRAD_poland_2<-as.data.frame(zero_tree_two_DRAD_poland[,-2])
zero_tree_three_DRAD_poland_2<-as.data.frame(zero_tree_three_DRAD_poland[,-2])
zero_tree_four_DRAD_poland_2<-as.data.frame(zero_tree_four_DRAD_poland[,-2])
zero_tree_five_DRAD_poland_2<-as.data.frame(zero_tree_five_DRAD_poland[,-2])

#merge individual chronologies
zero_chron_DRAD_poland.2 <- merge(zero_tree_one_DRAD_poland_2, zero_tree_two_DRAD_poland_2, all = TRUE)
zero_chron_DRAD_poland.2[is.na(zero_chron_DRAD_poland.2)] <- NA
zero_chron_DRAD_poland.2<- merge(zero_chron_DRAD_poland.2, zero_tree_three_DRAD_poland_2, all = TRUE)
zero_chron_DRAD_poland.2[is.na(zero_chron_DRAD_poland.2)] <- NA
zero_chron_DRAD_poland.2<- merge(zero_chron_DRAD_poland.2, zero_tree_four_DRAD_poland_2, all = TRUE)
zero_chron_DRAD_poland.2[is.na(zero_chron_DRAD_poland.2)] <- NA
zero_chron_DRAD_poland.2<- merge(zero_chron_DRAD_poland.2, zero_tree_five_DRAD_poland_2, all = TRUE)
zero_chron_DRAD_poland.2[is.na(zero_chron_DRAD_poland.2)] <- NA
zero_chron_DRAD_poland.3<-zero_chron_DRAD_poland.2[-1,]

#####
#RESULTS
#raw, standardised and mean chronologies
#####
zero_chron_DRAD_poland.4 <- as.data.frame(transform(zero_chron_DRAD_poland.3, zero_mean = rowMeans(zero_chron_DRAD_poland.3[,-1], na.rm = TRUE)))
zero_chron_DRAD_poland.scaled<-zero_chron_DRAD_poland.3[,-1]
zero_chron_DRAD_poland.scaled.2<- as.data.frame(scale(zero_chron_DRAD_poland.scaled))
zero_chron_DRAD_poland.scaled.3<- as.data.frame(transform(zero_chron_DRAD_poland.scaled.2, z_zero_mean = rowMeans(zero_chron_DRAD_poland.scaled.2, na.rm = TRUE)))
zero_chron_DRAD_poland_final<-cbind(zero_chron_DRAD_poland.4,zero_chron_DRAD_poland.scaled.3)
colnames(zero_chron_DRAD_poland_final)=c("year","zero_tree_one_DRAD_poland","zero_tree_two_DRAD_poland","zero_tree_three_DRAD_poland","zero_tree_four_DRAD_poland","zero_tree_five_DRAD_poland","zero_mean_DRAD_poland","z_zero_tree_one_DRAD_poland","z_zero_tree_two_DRAD_poland","z_zero_tree_three_DRAD_poland","z_zero_tree_four_DRAD_poland","z_zero_tree_five_DRAD_poland","z_zero_mean_DRAD_poland")
write.table(zero_chron_DRAD_poland_final, "zero_chron_DRAD_poland.txt")

#standardised plots
plot(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1969,2015), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
lines(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_one_DRAD_poland,col="#008597",lwd=0.75)
lines(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75)
lines(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75)
lines(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_four_DRAD_poland,col="#FF6093",lwd=0.75)
lines(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_five_DRAD_poland,col="#008597",lwd=0.75)
ylabs_1<-c("-2","-1","0","1","2","3","4")
cticks_1<-seq(-2, 4, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)

#####
# ten percentile, latewood section
#####
#remove precentile column
ten_tree_one_DRAD_poland_2<-as.data.frame(ten_tree_one_DRAD_poland[,-2])
ten_tree_two_DRAD_poland_2<-as.data.frame(ten_tree_two_DRAD_poland[,-2])
ten_tree_three_DRAD_poland_2<-as.data.frame(ten_tree_three_DRAD_poland[,-2])
ten_tree_four_DRAD_poland_2<-as.data.frame(ten_tree_four_DRAD_poland[,-2])
ten_tree_five_DRAD_poland_2<-as.data.frame(ten_tree_five_DRAD_poland[,-2])

#merge individual chronologies
ten_chron_DRAD_poland.2 <- merge(ten_tree_one_DRAD_poland_2, ten_tree_two_DRAD_poland_2, all = TRUE)
ten_chron_DRAD_poland.2[is.na(ten_chron_DRAD_poland.2)] <- NA
ten_chron_DRAD_poland.2<- merge(ten_chron_DRAD_poland.2, ten_tree_three_DRAD_poland_2, all = TRUE)
ten_chron_DRAD_poland.2[is.na(ten_chron_DRAD_poland.2)] <- NA
ten_chron_DRAD_poland.2<- merge(ten_chron_DRAD_poland.2, ten_tree_four_DRAD_poland_2, all = TRUE)
ten_chron_DRAD_poland.2[is.na(ten_chron_DRAD_poland.2)] <- NA
ten_chron_DRAD_poland.2<- merge(ten_chron_DRAD_poland.2, ten_tree_five_DRAD_poland_2, all = TRUE)
ten_chron_DRAD_poland.2[is.na(ten_chron_DRAD_poland.2)] <- NA
ten_chron_DRAD_poland.3<-ten_chron_DRAD_poland.2[-1,]

#####
#RESULTS
#raw, standardised and mean chronologies
#####
ten_chron_DRAD_poland.4 <- as.data.frame(transform(ten_chron_DRAD_poland.3, ten_mean = rowMeans(ten_chron_DRAD_poland.3[,-1], na.rm = TRUE)))
ten_chron_DRAD_poland.scaled<-ten_chron_DRAD_poland.3[,-1]
ten_chron_DRAD_poland.scaled.2<- as.data.frame(scale(ten_chron_DRAD_poland.scaled))
ten_chron_DRAD_poland.scaled.3<- as.data.frame(transform(ten_chron_DRAD_poland.scaled.2, z_ten_mean = rowMeans(ten_chron_DRAD_poland.scaled.2, na.rm = TRUE)))
ten_chron_DRAD_poland_final<-cbind(ten_chron_DRAD_poland.4,ten_chron_DRAD_poland.scaled.3)
colnames(ten_chron_DRAD_poland_final)=c("year","ten_tree_one_DRAD_poland","ten_tree_two_DRAD_poland","ten_tree_three_DRAD_poland","ten_tree_four_DRAD_poland","ten_tree_five_DRAD_poland","ten_mean_DRAD_poland","z_ten_tree_one_DRAD_poland","z_ten_tree_two_DRAD_poland","z_ten_tree_three_DRAD_poland","z_ten_tree_four_DRAD_poland","z_ten_tree_five_DRAD_poland","z_ten_mean_DRAD_poland")
write.table(ten_chron_DRAD_poland_final, "ten_chron_DRAD_poland.txt")

#standardised plots
plot(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1969,2015), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_one_DRAD_poland,col="#008597",lwd=0.75)
lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75)
lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75)
lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_four_DRAD_poland,col="#FF6093",lwd=0.75)
lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_five_DRAD_poland,col="#008597",lwd=0.75)
ylabs_1<-c("-2","-1","0","1","2","3","4")
cticks_1<-seq(-2, 4, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)

#####
# twenty percentile, latewood section
#####
#remove precentile column
twenty_tree_one_DRAD_poland_2<-as.data.frame(twenty_tree_one_DRAD_poland[,-2])
twenty_tree_two_DRAD_poland_2<-as.data.frame(twenty_tree_two_DRAD_poland[,-2])
twenty_tree_three_DRAD_poland_2<-as.data.frame(twenty_tree_three_DRAD_poland[,-2])
twenty_tree_four_DRAD_poland_2<-as.data.frame(twenty_tree_four_DRAD_poland[,-2])
twenty_tree_five_DRAD_poland_2<-as.data.frame(twenty_tree_five_DRAD_poland[,-2])

#merge individual chronologies
twenty_chron_DRAD_poland.2 <- merge(twenty_tree_one_DRAD_poland_2, twenty_tree_two_DRAD_poland_2, all = TRUE)
twenty_chron_DRAD_poland.2[is.na(twenty_chron_DRAD_poland.2)] <- NA
twenty_chron_DRAD_poland.2<- merge(twenty_chron_DRAD_poland.2, twenty_tree_three_DRAD_poland_2, all = TRUE)
twenty_chron_DRAD_poland.2[is.na(twenty_chron_DRAD_poland.2)] <- NA
twenty_chron_DRAD_poland.2<- merge(twenty_chron_DRAD_poland.2, twenty_tree_four_DRAD_poland_2, all = TRUE)
twenty_chron_DRAD_poland.2[is.na(twenty_chron_DRAD_poland.2)] <- NA
twenty_chron_DRAD_poland.2<- merge(twenty_chron_DRAD_poland.2, twenty_tree_five_DRAD_poland_2, all = TRUE)
twenty_chron_DRAD_poland.2[is.na(twenty_chron_DRAD_poland.2)] <- NA
twenty_chron_DRAD_poland.3<-twenty_chron_DRAD_poland.2[-1,]

#####
#RESULTS
#raw, standardised and mean chronologies
#####
twenty_chron_DRAD_poland.4 <- as.data.frame(transform(twenty_chron_DRAD_poland.3, twenty_mean = rowMeans(twenty_chron_DRAD_poland.3[,-1], na.rm = TRUE)))
twenty_chron_DRAD_poland.scaled<-twenty_chron_DRAD_poland.3[,-1]
twenty_chron_DRAD_poland.scaled.2<- as.data.frame(scale(twenty_chron_DRAD_poland.scaled))
twenty_chron_DRAD_poland.scaled.3<- as.data.frame(transform(twenty_chron_DRAD_poland.scaled.2, z_twenty_mean = rowMeans(twenty_chron_DRAD_poland.scaled.2, na.rm = TRUE)))
twenty_chron_DRAD_poland_final<-cbind(twenty_chron_DRAD_poland.4,twenty_chron_DRAD_poland.scaled.3)
colnames(twenty_chron_DRAD_poland_final)=c("year","twenty_tree_one_DRAD_poland","twenty_tree_two_DRAD_poland","twenty_tree_three_DRAD_poland","twenty_tree_four_DRAD_poland","twenty_tree_five_DRAD_poland","twenty_mean_DRAD_poland","z_twenty_tree_one_DRAD_poland","z_twenty_tree_two_DRAD_poland","z_twenty_tree_three_DRAD_poland","z_twenty_tree_four_DRAD_poland","z_twenty_tree_five_DRAD_poland","z_twenty_mean_DRAD_poland")
write.table(twenty_chron_DRAD_poland_final, "twenty_chron_DRAD_poland.txt")

#standardised plots
plot(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1969,2015), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_one_DRAD_poland,col="#008597",lwd=0.75)
lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75)
lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75)
lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75)
lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_five_DRAD_poland,col="#008597",lwd=0.75)
ylabs_1<-c("-2","-1","0","1","2","3","4")
cticks_1<-seq(-2, 4, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)

#####
# thirty percentile, latewood section
#####
#remove precentile column
thirty_tree_one_DRAD_poland_2<-as.data.frame(thirty_tree_one_DRAD_poland[,-2])
thirty_tree_two_DRAD_poland_2<-as.data.frame(thirty_tree_two_DRAD_poland[,-2])
thirty_tree_three_DRAD_poland_2<-as.data.frame(thirty_tree_three_DRAD_poland[,-2])
thirty_tree_four_DRAD_poland_2<-as.data.frame(thirty_tree_four_DRAD_poland[,-2])
thirty_tree_five_DRAD_poland_2<-as.data.frame(thirty_tree_five_DRAD_poland[,-2])

#merge individual chronologies
thirty_chron_DRAD_poland.2 <- merge(thirty_tree_one_DRAD_poland_2, thirty_tree_two_DRAD_poland_2, all = TRUE)
thirty_chron_DRAD_poland.2[is.na(thirty_chron_DRAD_poland.2)] <- NA
thirty_chron_DRAD_poland.2<- merge(thirty_chron_DRAD_poland.2, thirty_tree_three_DRAD_poland_2, all = TRUE)
thirty_chron_DRAD_poland.2[is.na(thirty_chron_DRAD_poland.2)] <- NA
thirty_chron_DRAD_poland.2<- merge(thirty_chron_DRAD_poland.2, thirty_tree_four_DRAD_poland_2, all = TRUE)
thirty_chron_DRAD_poland.2[is.na(thirty_chron_DRAD_poland.2)] <- NA
thirty_chron_DRAD_poland.2<- merge(thirty_chron_DRAD_poland.2, thirty_tree_five_DRAD_poland_2, all = TRUE)
thirty_chron_DRAD_poland.2[is.na(thirty_chron_DRAD_poland.2)] <- NA
thirty_chron_DRAD_poland.3<-thirty_chron_DRAD_poland.2[-1,]

#####
#RESULTS
#raw, standardised and mean chronologies
#####
thirty_chron_DRAD_poland.4 <- as.data.frame(transform(thirty_chron_DRAD_poland.3, thirty_mean = rowMeans(thirty_chron_DRAD_poland.3[,-1], na.rm = TRUE)))
thirty_chron_DRAD_poland.scaled<-thirty_chron_DRAD_poland.3[,-1]
thirty_chron_DRAD_poland.scaled.2<- as.data.frame(scale(thirty_chron_DRAD_poland.scaled))
thirty_chron_DRAD_poland.scaled.3<- as.data.frame(transform(thirty_chron_DRAD_poland.scaled.2, z_thirty_mean = rowMeans(thirty_chron_DRAD_poland.scaled.2, na.rm = TRUE)))
thirty_chron_DRAD_poland_final<-cbind(thirty_chron_DRAD_poland.4,thirty_chron_DRAD_poland.scaled.3)
colnames(thirty_chron_DRAD_poland_final)=c("year","thirty_tree_one_DRAD_poland","thirty_tree_two_DRAD_poland","thirty_tree_three_DRAD_poland","thirty_tree_four_DRAD_poland","thirty_tree_five_DRAD_poland","thirty_mean_DRAD_poland","z_thirty_tree_one_DRAD_poland","z_thirty_tree_two_DRAD_poland","z_thirty_tree_three_DRAD_poland","z_thirty_tree_four_DRAD_poland","z_thirty_tree_five_DRAD_poland","z_thirty_mean_DRAD_poland")
write.table(thirty_chron_DRAD_poland_final, "thirty_chron_DRAD_poland.txt")

#standardised plots
plot(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1969,2015), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_one_DRAD_poland,col="#008597",lwd=0.75)
lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75)
lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75)
lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75)
lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_five_DRAD_poland,col="#008597",lwd=0.75)
ylabs_1<-c("-2","-1","0","1","2","3","4")
cticks_1<-seq(-2, 4, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)

#####
# fourty percentile, latewood section
#####
#remove precentile column
fourty_tree_one_DRAD_poland_2<-as.data.frame(fourty_tree_one_DRAD_poland[,-2])
fourty_tree_two_DRAD_poland_2<-as.data.frame(fourty_tree_two_DRAD_poland[,-2])
fourty_tree_three_DRAD_poland_2<-as.data.frame(fourty_tree_three_DRAD_poland[,-2])
fourty_tree_four_DRAD_poland_2<-as.data.frame(fourty_tree_four_DRAD_poland[,-2])
fourty_tree_five_DRAD_poland_2<-as.data.frame(fourty_tree_five_DRAD_poland[,-2])

#merge individual chronologies
fourty_chron_DRAD_poland.2 <- merge(fourty_tree_one_DRAD_poland_2, fourty_tree_two_DRAD_poland_2, all = TRUE)
fourty_chron_DRAD_poland.2[is.na(fourty_chron_DRAD_poland.2)] <- NA
fourty_chron_DRAD_poland.2<- merge(fourty_chron_DRAD_poland.2, fourty_tree_three_DRAD_poland_2, all = TRUE)
fourty_chron_DRAD_poland.2[is.na(fourty_chron_DRAD_poland.2)] <- NA
fourty_chron_DRAD_poland.2<- merge(fourty_chron_DRAD_poland.2, fourty_tree_four_DRAD_poland_2, all = TRUE)
fourty_chron_DRAD_poland.2[is.na(fourty_chron_DRAD_poland.2)] <- NA
fourty_chron_DRAD_poland.2<- merge(fourty_chron_DRAD_poland.2, fourty_tree_five_DRAD_poland_2, all = TRUE)
fourty_chron_DRAD_poland.2[is.na(fourty_chron_DRAD_poland.2)] <- NA
fourty_chron_DRAD_poland.3<-fourty_chron_DRAD_poland.2[-1,]

#####
#RESULTS
#raw, standardised and mean chronologies
#####
fourty_chron_DRAD_poland.4 <- as.data.frame(transform(fourty_chron_DRAD_poland.3, fourty_mean = rowMeans(fourty_chron_DRAD_poland.3[,-1], na.rm = TRUE)))
fourty_chron_DRAD_poland.scaled<-fourty_chron_DRAD_poland.3[,-1]
fourty_chron_DRAD_poland.scaled.2<- as.data.frame(scale(fourty_chron_DRAD_poland.scaled))
fourty_chron_DRAD_poland.scaled.3<- as.data.frame(transform(fourty_chron_DRAD_poland.scaled.2, z_fourty_mean = rowMeans(fourty_chron_DRAD_poland.scaled.2, na.rm = TRUE)))
fourty_chron_DRAD_poland_final<-cbind(fourty_chron_DRAD_poland.4,fourty_chron_DRAD_poland.scaled.3)
colnames(fourty_chron_DRAD_poland_final)=c("year","fourty_tree_one_DRAD_poland","fourty_tree_two_DRAD_poland","fourty_tree_three_DRAD_poland","fourty_tree_four_DRAD_poland","fourty_tree_five_DRAD_poland","fourty_mean_DRAD_poland","z_fourty_tree_one_DRAD_poland","z_fourty_tree_two_DRAD_poland","z_fourty_tree_three_DRAD_poland","z_fourty_tree_four_DRAD_poland","z_fourty_tree_five_DRAD_poland","z_fourty_mean_DRAD_poland")
write.table(fourty_chron_DRAD_poland_final, "fourty_chron_DRAD_poland.txt")

#standardised plots
plot(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1969,2015), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_one_DRAD_poland,col="#008597",lwd=0.75)
lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75)
lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75)
lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75)
lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_five_DRAD_poland,col="#008597",lwd=0.75)
ylabs_1<-c("-2","-1","0","1","2","3","4")
cticks_1<-seq(-2, 4, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)

#####
# fifty percentile, latewood section
#####
#remove precentile column
fifty_tree_one_DRAD_poland_2<-as.data.frame(fifty_tree_one_DRAD_poland[,-2])
fifty_tree_two_DRAD_poland_2<-as.data.frame(fifty_tree_two_DRAD_poland[,-2])
fifty_tree_three_DRAD_poland_2<-as.data.frame(fifty_tree_three_DRAD_poland[,-2])
fifty_tree_four_DRAD_poland_2<-as.data.frame(fifty_tree_four_DRAD_poland[,-2])
fifty_tree_five_DRAD_poland_2<-as.data.frame(fifty_tree_five_DRAD_poland[,-2])

#merge individual chronologies
fifty_chron_DRAD_poland.2 <- merge(fifty_tree_one_DRAD_poland_2, fifty_tree_two_DRAD_poland_2, all = TRUE)
fifty_chron_DRAD_poland.2[is.na(fifty_chron_DRAD_poland.2)] <- NA
fifty_chron_DRAD_poland.2<- merge(fifty_chron_DRAD_poland.2, fifty_tree_three_DRAD_poland_2, all = TRUE)
fifty_chron_DRAD_poland.2[is.na(fifty_chron_DRAD_poland.2)] <- NA
fifty_chron_DRAD_poland.2<- merge(fifty_chron_DRAD_poland.2, fifty_tree_four_DRAD_poland_2, all = TRUE)
fifty_chron_DRAD_poland.2[is.na(fifty_chron_DRAD_poland.2)] <- NA
fifty_chron_DRAD_poland.2<- merge(fifty_chron_DRAD_poland.2, fifty_tree_five_DRAD_poland_2, all = TRUE)
fifty_chron_DRAD_poland.2[is.na(fifty_chron_DRAD_poland.2)] <- NA
fifty_chron_DRAD_poland.3<-fifty_chron_DRAD_poland.2[-1,]

#####
#RESULTS
#raw, standardised and mean chronologies
#####
fifty_chron_DRAD_poland.4 <- as.data.frame(transform(fifty_chron_DRAD_poland.3, fifty_mean = rowMeans(fifty_chron_DRAD_poland.3[,-1], na.rm = TRUE)))
fifty_chron_DRAD_poland.scaled<-fifty_chron_DRAD_poland.3[,-1]
fifty_chron_DRAD_poland.scaled.2<- as.data.frame(scale(fifty_chron_DRAD_poland.scaled))
fifty_chron_DRAD_poland.scaled.3<- as.data.frame(transform(fifty_chron_DRAD_poland.scaled.2, z_fifty_mean = rowMeans(fifty_chron_DRAD_poland.scaled.2, na.rm = TRUE)))
fifty_chron_DRAD_poland_final<-cbind(fifty_chron_DRAD_poland.4,fifty_chron_DRAD_poland.scaled.3)
colnames(fifty_chron_DRAD_poland_final)=c("year","fifty_tree_one_DRAD_poland","fifty_tree_two_DRAD_poland","fifty_tree_three_DRAD_poland","fifty_tree_four_DRAD_poland","fifty_tree_five_DRAD_poland","fifty_mean_DRAD_poland","z_fifty_tree_one_DRAD_poland","z_fifty_tree_two_DRAD_poland","z_fifty_tree_three_DRAD_poland","z_fifty_tree_four_DRAD_poland","z_fifty_tree_five_DRAD_poland","z_fifty_mean_DRAD_poland")
write.table(fifty_chron_DRAD_poland_final, "fifty_chron_DRAD_poland.txt")

#standardised plots
plot(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1969,2015), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_one_DRAD_poland,col="#008597",lwd=0.75)
lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75)
lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75)
lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75)
lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_five_DRAD_poland,col="#008597",lwd=0.75)
ylabs_1<-c("-2","-1","0","1","2","3","4")
cticks_1<-seq(-2, 4, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)

#####
# sixty percentile, latewood section
#####
#remove precentile column
sixty_tree_one_DRAD_poland_2<-as.data.frame(sixty_tree_one_DRAD_poland[,-2])
sixty_tree_two_DRAD_poland_2<-as.data.frame(sixty_tree_two_DRAD_poland[,-2])
sixty_tree_three_DRAD_poland_2<-as.data.frame(sixty_tree_three_DRAD_poland[,-2])
sixty_tree_four_DRAD_poland_2<-as.data.frame(sixty_tree_four_DRAD_poland[,-2])
sixty_tree_five_DRAD_poland_2<-as.data.frame(sixty_tree_five_DRAD_poland[,-2])

#merge individual chronologies
sixty_chron_DRAD_poland.2 <- merge(sixty_tree_one_DRAD_poland_2, sixty_tree_two_DRAD_poland_2, all = TRUE)
sixty_chron_DRAD_poland.2[is.na(sixty_chron_DRAD_poland.2)] <- NA
sixty_chron_DRAD_poland.2<- merge(sixty_chron_DRAD_poland.2, sixty_tree_three_DRAD_poland_2, all = TRUE)
sixty_chron_DRAD_poland.2[is.na(sixty_chron_DRAD_poland.2)] <- NA
sixty_chron_DRAD_poland.2<- merge(sixty_chron_DRAD_poland.2, sixty_tree_four_DRAD_poland_2, all = TRUE)
sixty_chron_DRAD_poland.2[is.na(sixty_chron_DRAD_poland.2)] <- NA
sixty_chron_DRAD_poland.2<- merge(sixty_chron_DRAD_poland.2, sixty_tree_five_DRAD_poland_2, all = TRUE)
sixty_chron_DRAD_poland.2[is.na(sixty_chron_DRAD_poland.2)] <- NA
sixty_chron_DRAD_poland.3<-sixty_chron_DRAD_poland.2[-1,]

#####
#RESULTS
#raw, standardised and mean chronologies
#####
sixty_chron_DRAD_poland.4 <- as.data.frame(transform(sixty_chron_DRAD_poland.3, sixty_mean = rowMeans(sixty_chron_DRAD_poland.3[,-1], na.rm = TRUE)))
sixty_chron_DRAD_poland.scaled<-sixty_chron_DRAD_poland.3[,-1]
sixty_chron_DRAD_poland.scaled.2<- as.data.frame(scale(sixty_chron_DRAD_poland.scaled))
sixty_chron_DRAD_poland.scaled.3<- as.data.frame(transform(sixty_chron_DRAD_poland.scaled.2, z_sixty_mean = rowMeans(sixty_chron_DRAD_poland.scaled.2, na.rm = TRUE)))
sixty_chron_DRAD_poland_final<-cbind(sixty_chron_DRAD_poland.4,sixty_chron_DRAD_poland.scaled.3)
colnames(sixty_chron_DRAD_poland_final)=c("year","sixty_tree_one_DRAD_poland","sixty_tree_two_DRAD_poland","sixty_tree_three_DRAD_poland","sixty_tree_four_DRAD_poland","sixty_tree_five_DRAD_poland","sixty_mean_DRAD_poland","z_sixty_tree_one_DRAD_poland","z_sixty_tree_two_DRAD_poland","z_sixty_tree_three_DRAD_poland","z_sixty_tree_four_DRAD_poland","z_sixty_tree_five_DRAD_poland","z_sixty_mean_DRAD_poland")
write.table(sixty_chron_DRAD_poland_final, "sixty_chron_DRAD_poland.txt")

#standardised plots
plot(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1969,2015), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_one_DRAD_poland,col="#008597",lwd=0.75)
lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75)
lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75)
lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75)
lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_five_DRAD_poland,col="#008597",lwd=0.75)
ylabs_1<-c("-2","-1","0","1","2","3","4")
cticks_1<-seq(-2, 4, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
#####
# seventy percentile, latewood section
#####
#remove precentile column
seventy_tree_one_DRAD_poland_2<-as.data.frame(seventy_tree_one_DRAD_poland[,-2])
seventy_tree_two_DRAD_poland_2<-as.data.frame(seventy_tree_two_DRAD_poland[,-2])
seventy_tree_three_DRAD_poland_2<-as.data.frame(seventy_tree_three_DRAD_poland[,-2])
seventy_tree_four_DRAD_poland_2<-as.data.frame(seventy_tree_four_DRAD_poland[,-2])
seventy_tree_five_DRAD_poland_2<-as.data.frame(seventy_tree_five_DRAD_poland[,-2])

#merge individual chronologies
seventy_chron_DRAD_poland.2 <- merge(seventy_tree_one_DRAD_poland_2, seventy_tree_two_DRAD_poland_2, all = TRUE)
seventy_chron_DRAD_poland.2[is.na(seventy_chron_DRAD_poland.2)] <- NA
seventy_chron_DRAD_poland.2<- merge(seventy_chron_DRAD_poland.2, seventy_tree_three_DRAD_poland_2, all = TRUE)
seventy_chron_DRAD_poland.2[is.na(seventy_chron_DRAD_poland.2)] <- NA
seventy_chron_DRAD_poland.2<- merge(seventy_chron_DRAD_poland.2, seventy_tree_four_DRAD_poland_2, all = TRUE)
seventy_chron_DRAD_poland.2[is.na(seventy_chron_DRAD_poland.2)] <- NA
seventy_chron_DRAD_poland.2<- merge(seventy_chron_DRAD_poland.2, seventy_tree_five_DRAD_poland_2, all = TRUE)
seventy_chron_DRAD_poland.2[is.na(seventy_chron_DRAD_poland.2)] <- NA
seventy_chron_DRAD_poland.3<-seventy_chron_DRAD_poland.2[-1,]

#####
#RESULTS
#raw, standardised and mean chronologies
#####
seventy_chron_DRAD_poland.4 <- as.data.frame(transform(seventy_chron_DRAD_poland.3, seventy_mean = rowMeans(seventy_chron_DRAD_poland.3[,-1], na.rm = TRUE)))
seventy_chron_DRAD_poland.scaled<-seventy_chron_DRAD_poland.3[,-1]
seventy_chron_DRAD_poland.scaled.2<- as.data.frame(scale(seventy_chron_DRAD_poland.scaled))
seventy_chron_DRAD_poland.scaled.3<- as.data.frame(transform(seventy_chron_DRAD_poland.scaled.2, z_seventy_mean = rowMeans(seventy_chron_DRAD_poland.scaled.2, na.rm = TRUE)))
seventy_chron_DRAD_poland_final<-cbind(seventy_chron_DRAD_poland.4,seventy_chron_DRAD_poland.scaled.3)
colnames(seventy_chron_DRAD_poland_final)=c("year","seventy_tree_one_DRAD_poland","seventy_tree_two_DRAD_poland","seventy_tree_three_DRAD_poland","seventy_tree_four_DRAD_poland","seventy_tree_five_DRAD_poland","seventy_mean_DRAD_poland","z_seventy_tree_one_DRAD_poland","z_seventy_tree_two_DRAD_poland","z_seventy_tree_three_DRAD_poland","z_seventy_tree_four_DRAD_poland","z_seventy_tree_five_DRAD_poland","z_seventy_mean_DRAD_poland")
write.table(seventy_chron_DRAD_poland_final, "seventy_chron_DRAD_poland.txt")

#standardised plots
plot(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1969,2015), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_one_DRAD_poland,col="#008597",lwd=0.75)
lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75)
lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75)
lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_four_DRAD_poland,col="#FF6093",lwd=0.75)
lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_five_DRAD_poland,col="#008597",lwd=0.75)
ylabs_1<-c("-2","-1","0","1","2","3","4")
cticks_1<-seq(-2, 4, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)

#####
# seventyfive percentile, latewood section
#####
#remove precentile column
seventyfive_tree_one_DRAD_poland_2<-as.data.frame(seventyfive_tree_one_DRAD_poland[,-2])
seventyfive_tree_two_DRAD_poland_2<-as.data.frame(seventyfive_tree_two_DRAD_poland[,-2])
seventyfive_tree_three_DRAD_poland_2<-as.data.frame(seventyfive_tree_three_DRAD_poland[,-2])
seventyfive_tree_four_DRAD_poland_2<-as.data.frame(seventyfive_tree_four_DRAD_poland[,-2])
seventyfive_tree_five_DRAD_poland_2<-as.data.frame(seventyfive_tree_five_DRAD_poland[,-2])

#merge individual chronologies
seventyfive_chron_DRAD_poland.2 <- merge(seventyfive_tree_one_DRAD_poland_2, seventyfive_tree_two_DRAD_poland_2, all = TRUE)
seventyfive_chron_DRAD_poland.2[is.na(seventyfive_chron_DRAD_poland.2)] <- NA
seventyfive_chron_DRAD_poland.2<- merge(seventyfive_chron_DRAD_poland.2, seventyfive_tree_three_DRAD_poland_2, all = TRUE)
seventyfive_chron_DRAD_poland.2[is.na(seventyfive_chron_DRAD_poland.2)] <- NA
seventyfive_chron_DRAD_poland.2<- merge(seventyfive_chron_DRAD_poland.2, seventyfive_tree_four_DRAD_poland_2, all = TRUE)
seventyfive_chron_DRAD_poland.2[is.na(seventyfive_chron_DRAD_poland.2)] <- NA
seventyfive_chron_DRAD_poland.2<- merge(seventyfive_chron_DRAD_poland.2, seventyfive_tree_five_DRAD_poland_2, all = TRUE)
seventyfive_chron_DRAD_poland.2[is.na(seventyfive_chron_DRAD_poland.2)] <- NA
seventyfive_chron_DRAD_poland.3<-seventyfive_chron_DRAD_poland.2[-1,]

#####
#RESULTS
#raw, standardised and mean chronologies
#####
seventyfive_chron_DRAD_poland.4 <- as.data.frame(transform(seventyfive_chron_DRAD_poland.3, seventyfive_mean = rowMeans(seventyfive_chron_DRAD_poland.3[,-1], na.rm = TRUE)))
seventyfive_chron_DRAD_poland.scaled<-seventyfive_chron_DRAD_poland.3[,-1]
seventyfive_chron_DRAD_poland.scaled.2<- as.data.frame(scale(seventyfive_chron_DRAD_poland.scaled))
seventyfive_chron_DRAD_poland.scaled.3<- as.data.frame(transform(seventyfive_chron_DRAD_poland.scaled.2, z_seventyfive_mean = rowMeans(seventyfive_chron_DRAD_poland.scaled.2, na.rm = TRUE)))
seventyfive_chron_DRAD_poland_final<-cbind(seventyfive_chron_DRAD_poland.4,seventyfive_chron_DRAD_poland.scaled.3)
colnames(seventyfive_chron_DRAD_poland_final)=c("year","seventyfive_tree_one_DRAD_poland","seventyfive_tree_two_DRAD_poland","seventyfive_tree_three_DRAD_poland","seventyfive_tree_four_DRAD_poland","seventyfive_tree_five_DRAD_poland","seventyfive_mean_DRAD_poland","z_seventyfive_tree_one_DRAD_poland","z_seventyfive_tree_two_DRAD_poland","z_seventyfive_tree_three_DRAD_poland","z_seventyfive_tree_four_DRAD_poland","z_seventyfive_tree_five_DRAD_poland","z_seventyfive_mean_DRAD_poland")
write.table(seventyfive_chron_DRAD_poland_final, "seventyfive_chron_DRAD_poland.txt")

#standardised plots
plot(seventyfive_chron_DRAD_poland_final$year, seventyfive_chron_DRAD_poland_final$z_seventyfive_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1969,2015), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
lines(seventyfive_chron_DRAD_poland_final$year, seventyfive_chron_DRAD_poland_final$z_seventyfive_tree_one_DRAD_poland,col="#008597",lwd=0.75)
lines(seventyfive_chron_DRAD_poland_final$year, seventyfive_chron_DRAD_poland_final$z_seventyfive_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75)
lines(seventyfive_chron_DRAD_poland_final$year, seventyfive_chron_DRAD_poland_final$z_seventyfive_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75)
lines(seventyfive_chron_DRAD_poland_final$year, seventyfive_chron_DRAD_poland_final$z_seventyfive_tree_four_DRAD_poland,col="#FF6093",lwd=0.75)
lines(seventyfive_chron_DRAD_poland_final$year, seventyfive_chron_DRAD_poland_final$z_seventyfive_tree_five_DRAD_poland,col="#008597",lwd=0.75)
ylabs_1<-c("-2","-1","0","1","2","3","4")
cticks_1<-seq(-2, 4, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)


#####
# eighty percentile, latewood section
#####
#remove precentile column
eighty_tree_one_DRAD_poland_2<-as.data.frame(eighty_tree_one_DRAD_poland[,-2])
eighty_tree_two_DRAD_poland_2<-as.data.frame(eighty_tree_two_DRAD_poland[,-2])
eighty_tree_three_DRAD_poland_2<-as.data.frame(eighty_tree_three_DRAD_poland[,-2])
eighty_tree_four_DRAD_poland_2<-as.data.frame(eighty_tree_four_DRAD_poland[,-2])
eighty_tree_five_DRAD_poland_2<-as.data.frame(eighty_tree_five_DRAD_poland[,-2])

#merge individual chronologies
eighty_chron_DRAD_poland.2 <- merge(eighty_tree_one_DRAD_poland_2, eighty_tree_two_DRAD_poland_2, all = TRUE)
eighty_chron_DRAD_poland.2[is.na(eighty_chron_DRAD_poland.2)] <- NA
eighty_chron_DRAD_poland.2<- merge(eighty_chron_DRAD_poland.2, eighty_tree_three_DRAD_poland_2, all = TRUE)
eighty_chron_DRAD_poland.2[is.na(eighty_chron_DRAD_poland.2)] <- NA
eighty_chron_DRAD_poland.2<- merge(eighty_chron_DRAD_poland.2, eighty_tree_four_DRAD_poland_2, all = TRUE)
eighty_chron_DRAD_poland.2[is.na(eighty_chron_DRAD_poland.2)] <- NA
eighty_chron_DRAD_poland.2<- merge(eighty_chron_DRAD_poland.2, eighty_tree_five_DRAD_poland_2, all = TRUE)
eighty_chron_DRAD_poland.2[is.na(eighty_chron_DRAD_poland.2)] <- NA
eighty_chron_DRAD_poland.3<-eighty_chron_DRAD_poland.2[-1,]

#####
#RESULTS
#raw, standardised and mean chronologies
#####
eighty_chron_DRAD_poland.4 <- as.data.frame(transform(eighty_chron_DRAD_poland.3, eighty_mean = rowMeans(eighty_chron_DRAD_poland.3[,-1], na.rm = TRUE)))
eighty_chron_DRAD_poland.scaled<-eighty_chron_DRAD_poland.3[,-1]
eighty_chron_DRAD_poland.scaled.2<- as.data.frame(scale(eighty_chron_DRAD_poland.scaled))
eighty_chron_DRAD_poland.scaled.3<- as.data.frame(transform(eighty_chron_DRAD_poland.scaled.2, z_eighty_mean = rowMeans(eighty_chron_DRAD_poland.scaled.2, na.rm = TRUE)))
eighty_chron_DRAD_poland_final<-cbind(eighty_chron_DRAD_poland.4,eighty_chron_DRAD_poland.scaled.3)
colnames(eighty_chron_DRAD_poland_final)=c("year","eighty_tree_one_DRAD_poland","eighty_tree_two_DRAD_poland","eighty_tree_three_DRAD_poland","eighty_tree_four_DRAD_poland","eighty_tree_five_DRAD_poland","eighty_mean_DRAD_poland","z_eighty_tree_one_DRAD_poland","z_eighty_tree_two_DRAD_poland","z_eighty_tree_three_DRAD_poland","z_eighty_tree_four_DRAD_poland","z_eighty_tree_five_DRAD_poland","z_eighty_mean_DRAD_poland")
write.table(eighty_chron_DRAD_poland_final, "eighty_chron_DRAD_poland.txt")

#standardised plots
plot(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1969,2015), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_one_DRAD_poland,col="#008597",lwd=0.75)
lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75)
lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75)
lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75)
lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_five_DRAD_poland,col="#008597",lwd=0.75)
ylabs_1<-c("-2","-1","0","1","2","3","4")
cticks_1<-seq(-2, 4, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
#####
# ninety percentile, latewood section
#####
#remove precentile column
ninety_tree_one_DRAD_poland_2<-as.data.frame(ninety_tree_one_DRAD_poland[,-2])
ninety_tree_two_DRAD_poland_2<-as.data.frame(ninety_tree_two_DRAD_poland[,-2])
ninety_tree_three_DRAD_poland_2<-as.data.frame(ninety_tree_three_DRAD_poland[,-2])
ninety_tree_four_DRAD_poland_2<-as.data.frame(ninety_tree_four_DRAD_poland[,-2])
ninety_tree_five_DRAD_poland_2<-as.data.frame(ninety_tree_five_DRAD_poland[,-2])

#merge individual chronologies
ninety_chron_DRAD_poland.2 <- merge(ninety_tree_one_DRAD_poland_2, ninety_tree_two_DRAD_poland_2, all = TRUE)
ninety_chron_DRAD_poland.2[is.na(ninety_chron_DRAD_poland.2)] <- NA
ninety_chron_DRAD_poland.2<- merge(ninety_chron_DRAD_poland.2, ninety_tree_three_DRAD_poland_2, all = TRUE)
ninety_chron_DRAD_poland.2[is.na(ninety_chron_DRAD_poland.2)] <- NA
ninety_chron_DRAD_poland.2<- merge(ninety_chron_DRAD_poland.2, ninety_tree_four_DRAD_poland_2, all = TRUE)
ninety_chron_DRAD_poland.2[is.na(ninety_chron_DRAD_poland.2)] <- NA
ninety_chron_DRAD_poland.2<- merge(ninety_chron_DRAD_poland.2, ninety_tree_five_DRAD_poland_2, all = TRUE)
ninety_chron_DRAD_poland.2[is.na(ninety_chron_DRAD_poland.2)] <- NA
ninety_chron_DRAD_poland.3<-ninety_chron_DRAD_poland.2[-1,]

#####
#RESULTS
#raw, standardised and mean chronologies
#####
ninety_chron_DRAD_poland.4 <- as.data.frame(transform(ninety_chron_DRAD_poland.3, ninety_mean = rowMeans(ninety_chron_DRAD_poland.3[,-1], na.rm = TRUE)))
ninety_chron_DRAD_poland.scaled<-ninety_chron_DRAD_poland.3[,-1]
ninety_chron_DRAD_poland.scaled.2<- as.data.frame(scale(ninety_chron_DRAD_poland.scaled))
ninety_chron_DRAD_poland.scaled.3<- as.data.frame(transform(ninety_chron_DRAD_poland.scaled.2, z_ninety_mean = rowMeans(ninety_chron_DRAD_poland.scaled.2, na.rm = TRUE)))
ninety_chron_DRAD_poland_final<-cbind(ninety_chron_DRAD_poland.4,ninety_chron_DRAD_poland.scaled.3)
colnames(ninety_chron_DRAD_poland_final)=c("year","ninety_tree_one_DRAD_poland","ninety_tree_two_DRAD_poland","ninety_tree_three_DRAD_poland","ninety_tree_four_DRAD_poland","ninety_tree_five_DRAD_poland","ninety_mean_DRAD_poland","z_ninety_tree_one_DRAD_poland","z_ninety_tree_two_DRAD_poland","z_ninety_tree_three_DRAD_poland","z_ninety_tree_four_DRAD_poland","z_ninety_tree_five_DRAD_poland","z_ninety_mean_DRAD_poland")
write.table(ninety_chron_DRAD_poland_final, "ninety_chron_DRAD_poland.txt")

#standardised plots
plot(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1969,2015), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_one_DRAD_poland,col="#008597",lwd=0.75)
lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75)
lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75)
lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_four_DRAD_poland,col="#FF6093",lwd=0.75)
lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_five_DRAD_poland,col="#008597",lwd=0.75)
ylabs_1<-c("-2","-1","0","1","2","3","4")
cticks_1<-seq(-2, 4, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
#####
# hundred percentile, latewood section
#####
#remove precentile column
hundred_tree_one_DRAD_poland_2<-as.data.frame(hundred_tree_one_DRAD_poland[,-2])
hundred_tree_two_DRAD_poland_2<-as.data.frame(hundred_tree_two_DRAD_poland[,-2])
hundred_tree_three_DRAD_poland_2<-as.data.frame(hundred_tree_three_DRAD_poland[,-2])
hundred_tree_four_DRAD_poland_2<-as.data.frame(hundred_tree_four_DRAD_poland[,-2])
hundred_tree_five_DRAD_poland_2<-as.data.frame(hundred_tree_five_DRAD_poland[,-2])

#merge individual chronologies
hundred_chron_DRAD_poland.2 <- merge(hundred_tree_one_DRAD_poland_2, hundred_tree_two_DRAD_poland_2, all = TRUE)
hundred_chron_DRAD_poland.2[is.na(hundred_chron_DRAD_poland.2)] <- NA
hundred_chron_DRAD_poland.2<- merge(hundred_chron_DRAD_poland.2, hundred_tree_three_DRAD_poland_2, all = TRUE)
hundred_chron_DRAD_poland.2[is.na(hundred_chron_DRAD_poland.2)] <- NA
hundred_chron_DRAD_poland.2<- merge(hundred_chron_DRAD_poland.2, hundred_tree_four_DRAD_poland_2, all = TRUE)
hundred_chron_DRAD_poland.2[is.na(hundred_chron_DRAD_poland.2)] <- NA
hundred_chron_DRAD_poland.2<- merge(hundred_chron_DRAD_poland.2, hundred_tree_five_DRAD_poland_2, all = TRUE)
hundred_chron_DRAD_poland.2[is.na(hundred_chron_DRAD_poland.2)] <- NA
hundred_chron_DRAD_poland.3<-hundred_chron_DRAD_poland.2[-1,]

#####
#RESULTS
#raw, standardised and mean chronologies
#####
hundred_chron_DRAD_poland.4 <- as.data.frame(transform(hundred_chron_DRAD_poland.3, hundred_mean = rowMeans(hundred_chron_DRAD_poland.3[,-1], na.rm = TRUE)))
hundred_chron_DRAD_poland.scaled<-hundred_chron_DRAD_poland.3[,-1]
hundred_chron_DRAD_poland.scaled.2<- as.data.frame(scale(hundred_chron_DRAD_poland.scaled))
hundred_chron_DRAD_poland.scaled.3<- as.data.frame(transform(hundred_chron_DRAD_poland.scaled.2, z_hundred_mean = rowMeans(hundred_chron_DRAD_poland.scaled.2, na.rm = TRUE)))
hundred_chron_DRAD_poland_final<-cbind(hundred_chron_DRAD_poland.4,hundred_chron_DRAD_poland.scaled.3)
colnames(hundred_chron_DRAD_poland_final)=c("year","hundred_tree_one_DRAD_poland","hundred_tree_two_DRAD_poland","hundred_tree_three_DRAD_poland","hundred_tree_four_DRAD_poland","hundred_tree_five_DRAD_poland","hundred_mean_DRAD_poland","z_hundred_tree_one_DRAD_poland","z_hundred_tree_two_DRAD_poland","z_hundred_tree_three_DRAD_poland","z_hundred_tree_four_DRAD_poland","z_hundred_tree_five_DRAD_poland","z_hundred_mean_DRAD_poland")
write.table(hundred_chron_DRAD_poland_final, "hundred_chron_DRAD_poland.txt")

#standardised plots
plot(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1969,2015), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_one_DRAD_poland,col="#008597",lwd=0.75)
lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75)
lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75)
lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_four_DRAD_poland,col="#FF6093",lwd=0.75)
lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_five_DRAD_poland,col="#008597",lwd=0.75)
ylabs_1<-c("-2","-1","0","1","2","3","4")
cticks_1<-seq(-2, 4, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)

#####
# chronologies output
#####
year<-read.table("year_poland.txt", header=TRUE, sep="\t")
z_zero_mean_DRAD_poland = rowMeans(zero_chron_DRAD_poland.scaled.2, na.rm = TRUE)
z_ten_mean_DRAD_poland = rowMeans(ten_chron_DRAD_poland.scaled.2, na.rm = TRUE)
z_twenty_mean_DRAD_poland = rowMeans(twenty_chron_DRAD_poland.scaled.2, na.rm = TRUE)
z_thirty_mean_DRAD_poland = rowMeans(thirty_chron_DRAD_poland.scaled.2, na.rm = TRUE)
z_fourty_mean_DRAD_poland = rowMeans(fourty_chron_DRAD_poland.scaled.2, na.rm = TRUE)
z_fifty_mean_DRAD_poland = rowMeans(fifty_chron_DRAD_poland.scaled.2, na.rm = TRUE)
z_sixty_mean_DRAD_poland = rowMeans(sixty_chron_DRAD_poland.scaled.2, na.rm = TRUE)
z_seventy_mean_DRAD_poland = rowMeans(seventy_chron_DRAD_poland.scaled.2, na.rm = TRUE)
z_eighty_mean_DRAD_poland = rowMeans(eighty_chron_DRAD_poland.scaled.2, na.rm = TRUE)
z_ninety_mean_DRAD_poland = rowMeans(ninety_chron_DRAD_poland.scaled.2, na.rm = TRUE)
z_hundred_mean_DRAD_poland = rowMeans(hundred_chron_DRAD_poland.scaled.2, na.rm = TRUE)

DRAD_mean_chrons_poland<-cbind(year$year,z_zero_mean_DRAD_poland,z_ten_mean_DRAD_poland,z_twenty_mean_DRAD_poland,z_thirty_mean_DRAD_poland,z_fourty_mean_DRAD_poland,z_fifty_mean_DRAD_poland,z_sixty_mean_DRAD_poland,z_seventy_mean_DRAD_poland,z_eighty_mean_DRAD_poland,z_ninety_mean_DRAD_poland,z_hundred_mean_DRAD_poland)
colnames(DRAD_mean_chrons_poland)<-c("year","zero_mean_DRAD","ten_mean_mean_DRAD","twenty_mean_DRAD","thirty_mean_DRAD","fourty_mean_DRAD","fifty_mean_DRAD","sixty_mean_DRAD","seventy_mean_DRAD","eighty_mean_DRAD","ninety_mean_DRAD","hundred_mean_DRAD")
write.table(DRAD_mean_chrons_poland, "DRAD_mean_chrons_poland.txt")

###########
#PLOTS
###########
dev.off()
pdf("DRAD_poland_plots.pdf",width=3)
par(mar=c(0.3,1.5,0.1,1.5), oma=c(3,0.5,0.5,0.5),mfcol=c(11,1))

plot(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(zero_chron_final$year, zero_chron_final$z_zero_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"zero_percentile_chron",cex=0.5,adj=0)  
  
plot(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"ten_percentile_chron",cex=0.5,adj=0)  

plot(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"twenty_percentile_chron",cex=0.5,adj=0)  

plot(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"thirty_percentile_chron",cex=0.5,adj=0)  

plot(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"fourty_percentile_chron",cex=0.5,adj=0)  

plot(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"fifty_percentile_chron",cex=0.5,adj=0)  

plot(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"sixty_percentile_chron",cex=0.5,adj=0)  

plot(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"seventy_percentile_chron",cex=0.5,adj=0)  

plot(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"eighty_percentile_chron",cex=0.5,adj=0)  

plot(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"ninety_percentile_chron",cex=0.5,adj=0)  

plot(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
xticks_2<-seq(1950, 2010, by=5)
axis(1, at=xticks_2, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"hundred_percentile_chron",cex=0.5,adj=0)  
dev.off()


###########
#PLOTS
###########
dev.off()
pdf("DRAD_poland_plots.pdf",width=3)
par(mar=c(0.3,1.5,0.1,1.5), oma=c(3,0.5,0.5,0.5),mfcol=c(11,1))

plot(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(zero_chron_final$year, zero_chron_final$z_zero_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"zero_percentile_chron",cex=0.5,adj=0)  
  
plot(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"ten_percentile_chron",cex=0.5,adj=0)  

plot(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"twenty_percentile_chron",cex=0.5,adj=0)  

plot(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"thirty_percentile_chron",cex=0.5,adj=0)  

plot(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"fourty_percentile_chron",cex=0.5,adj=0)  

plot(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"fifty_percentile_chron",cex=0.5,adj=0)  

plot(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"sixty_percentile_chron",cex=0.5,adj=0)  

plot(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"seventy_percentile_chron",cex=0.5,adj=0)  

plot(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"eighty_percentile_chron",cex=0.5,adj=0)  

plot(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"ninety_percentile_chron",cex=0.5,adj=0)  

plot(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
xticks_2<-seq(1950, 2010, by=5)
axis(1, at=xticks_2, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"hundred_percentile_chron",cex=0.5,adj=0)  
dev.off()


###########
#PLOTS
###########
dev.off()
pdf("DRAD_poland_plots.pdf",width=3)
par(mar=c(0.3,1.5,0.1,1.5), oma=c(3,0.5,0.5,0.5),mfcol=c(11,1))

plot(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(zero_chron_DRAD_poland_final$year, zero_chron_DRAD_poland_final$z_zero_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(zero_chron_final$year, zero_chron_final$z_zero_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"zero_percentile_chron",cex=0.5,adj=0)  
  
plot(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(ten_chron_DRAD_poland_final$year, ten_chron_DRAD_poland_final$z_ten_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"ten_percentile_chron",cex=0.5,adj=0)  

plot(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(twenty_chron_DRAD_poland_final$year, twenty_chron_DRAD_poland_final$z_twenty_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"twenty_percentile_chron",cex=0.5,adj=0)  

plot(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(thirty_chron_DRAD_poland_final$year, thirty_chron_DRAD_poland_final$z_thirty_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"thirty_percentile_chron",cex=0.5,adj=0)  

plot(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(fourty_chron_DRAD_poland_final$year, fourty_chron_DRAD_poland_final$z_fourty_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"fourty_percentile_chron",cex=0.5,adj=0)  

plot(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(fifty_chron_DRAD_poland_final$year, fifty_chron_DRAD_poland_final$z_fifty_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"fifty_percentile_chron",cex=0.5,adj=0)  

plot(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(sixty_chron_DRAD_poland_final$year, sixty_chron_DRAD_poland_final$z_sixty_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"sixty_percentile_chron",cex=0.5,adj=0)  

plot(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(seventy_chron_DRAD_poland_final$year, seventy_chron_DRAD_poland_final$z_seventy_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"seventy_percentile_chron",cex=0.5,adj=0)  

plot(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(eighty_chron_DRAD_poland_final$year, eighty_chron_DRAD_poland_final$z_eighty_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"eighty_percentile_chron",cex=0.5,adj=0)  

plot(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(ninety_chron_DRAD_poland_final$year, ninety_chron_DRAD_poland_final$z_ninety_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"ninety_percentile_chron",cex=0.5,adj=0)  

plot(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_one_DRAD_poland, type="n", ylim=c(-3.5,4),xlim=c(1950,2010), xlab="", ylab="",xaxt="n",axes=FALSE,xaxs="i")
abline(h=0, lty=2, lwd=0.5)
lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_one_DRAD_poland,col="#008597",lwd=0.75,lty=1)
lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_two_DRAD_poland,col="#7EFBE3",lwd=0.75,lty=1)
lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_three_DRAD_poland,col="#FF95AE",lwd=0.75,lty=1)
lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_four_DRAD_poland,col="#FF6093",lwd=0.75,lty=1)
lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_tree_five_DRAD_poland,col="#423828",lwd=0.75,lty=1)
#lines(hundred_chron_DRAD_poland_final$year, hundred_chron_DRAD_poland_final$z_hundred_mean,col="#3F1500", lwd=1.25)
ylabs_1<-c("-3","-2","-1","0","1","2","3")
cticks_1<-seq(-3, 3, by=1)
axis(2, at=cticks_1, ylabs_1, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
xticks_2<-seq(1950, 2010, by=5)
axis(1, at=xticks_2, line=0.5,cex.axis=0.5,las=2,col="black",col.axis="black", mgp=c(0, 0.6, 0),lwd=0.5,tck=-0.05)
text(1950,3.5,"hundred_percentile_chron",cex=0.5,adj=0)  
dev.off()
