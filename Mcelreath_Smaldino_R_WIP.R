library(rethinking)

rm(list=ls())  

#parameters for the world
beta <- 0.2
power <- 1 - beta
alpha <- c(0.05)
trueneg <- 1 - alpha
baserates <- c(0.005)
p_novel <- 0.8

n <- 20 #hypotheses tested each round
rounds <- 5000
total_h <- n * rounds

#communication parameters
pos_nov_com <- 1
pos_rep_com <- c(1)
neg_nov_com <- 0
neg_rep_com <- 1

#FOR DIFFERENT BASERATES and communication parameters
for(alpha in alpha){
for(baserate in baserates) {
  for(PR_COM in pos_rep_com) {
  
  novel_h <- rbinom(total_h, 1, baserate)
  novel_h_id <- 1:length(novel_h)
  
  #tested hypotheses
  r_df <- data.frame(id = rep(NA, total_h),
                     true = rep(NA, total_h),
                     result = rep(NA, total_h), 
                     comm = rep(NA, total_h), 
                     tally = rep(0, total_h))
  
  #data frame to store what happens each round
  q_df <- data.frame(id = rep(NA, n),
                     true = rep(NA, n),
                     novel = rep(NA, n), 
                     result = rep(NA, n), 
                     comm = rep(NA, n), 
                     tally = rep(0, n))
  
  #FOR MANY RUNS
  for(runs in 1:rounds){
    
    #if no hypotheses to be replicated, then samples only novel hypotheses
    if(is.na(r_df$true[r_df$comm==1][1])){
      q_df$true <- novel_h[1:n]
      #remove already sampled hypotheses
      novel_h <- novel_h[-c(1:n)]
      #update novel column in q_df
      q_df$novel <- 1 
    } else{
      
      #data frame of communicated hypotheses
      communicated_df <- subset(r_df, comm == 1)
      
      for(i in 1:n) {
        if(runif(1, 0, 1) < p_novel){
          q_df$true[i] <- novel_h[1]
          q_df$novel[i] <- 1 
          #remove already sampled hypothesis
          novel_h <- novel_h[-1]
        } else {
          #generates position from those positions that contain tested, communicated hypotheses
          position <- sample(which(!is.na(communicated_df$id)), 1, replace = TRUE)
          #adds previously tested hypotheses to the q_df data frame
          q_df$id[i] <- communicated_df$id[position]
          q_df$true[i] <- communicated_df$true[position]
          q_df$novel[i] <- 0
        } }
    } 
    
    #assign id's to novel hypotheses
    xx <- length(q_df$id[q_df$novel == 1])
    q_df$id[q_df$novel == 1] <- novel_h_id[1:xx]
    #remove already used id's from novel_h_id
    novel_h_id <- novel_h_id[-c(1:xx)]
    
    #TEST HYPOTHESES
    #results for false hypotheses: positive result with probability alpha
    q_df <- within(q_df, result[true == 0] <- rbinom(length(result[true == 0]), 1, alpha))
    #results for true hypotheses: positive result with probability power
    q_df <- within(q_df, result[true == 1] <- rbinom(length(result[true == 1]), 1, power))
    
    #DETERMINE WHETHER HYPOTHESES ARE COMMUNICATED
    #novel results
    q_df <- within(q_df, comm[novel==1 & result == 1] <- rbinom(length(comm[novel==1 & result == 1]), 1, pos_nov_com))
    q_df <- within(q_df, comm[novel==1 & result == 0] <- rbinom(length(comm[novel==1 & result == 0]), 1, neg_nov_com))
    #replications results
    q_df <- within(q_df, comm[novel==0 & result == 1] <- rbinom(length(comm[novel==0 & result == 1]), 1, PR_COM))
    q_df <- within(q_df, comm[novel==0 & result == 0] <- rbinom(length(comm[novel==0 & result == 0]), 1, neg_rep_com))
    
    #UPDATE TALLY
    q_df$tally[q_df$comm == 1 & q_df$result == 1] <- 1
    q_df$tally[q_df$comm == 1 & q_df$result == 0] <- -1
    q_df$tally[q_df$comm == 0] <- 0
    
    #ADD HYPOTHESES TO POOL OF COMMUNICATED HYPOTHESES
    
    #update data in r_df data frame
    r_spot <- match(NA, r_df$id)
    r_spot2 <- r_spot + nrow(q_df) - 1
    r_df$id[r_spot:r_spot2] <- q_df$id
    r_df$true[r_spot:r_spot2] <- q_df$true
    r_df$result[r_spot:r_spot2] <- q_df$result
    r_df$comm[r_spot:r_spot2] <- q_df$comm
    r_df$tally[r_spot:r_spot2] <- q_df$tally
  }
  
  ##aggregate data frame for each unique hypothesis id, for communicated results only, because only these influence the tally
  r_df_comm <- r_df[r_df$comm==1,]
  final_df <- aggregate(tally ~ id + true, data=r_df_comm, FUN=sum)
  
  ##plot
  simplehist(final_df$tally[final_df$true==0], ylab = "Number of Hypotheses", xlab = "Tally", 
             main = paste("False Hyp., +nov:", pos_nov_com, ", +rep", PR_COM, ", -nov", neg_nov_com, 
                          ", -rep", neg_rep_com, ", Base Rate = ", baserate, "alpha:", alpha))
  simplehist(final_df$tally[final_df$true==1], ylab = "Number of Hypotheses", xlab = "Tally", 
             main = paste("True Hyp., +nov:", pos_nov_com, ", +rep", PR_COM, ", -nov", neg_nov_com, 
                          ", -rep", neg_rep_com, ", Base Rate = ", baserate, "alpha:", alpha))
  
  ##precision: proportion of hypotheses at each tally that are true
  precision_df <- aggregate(true ~ tally, data=final_df, FUN = function(x) sum(x) / length(x))
  
  gg <- ggplot(precision_df, aes(x = tally)) + geom_histogram(aes(y = true), stat = "identity") + 
    scale_y_continuous(name = "Proportion of True Hypotheses at Tally") + 
    scale_x_continuous(limits=c(-5, 10)) +
    ggtitle(paste("Comm Param = +nov:", pos_nov_com, ", +rep", PR_COM, ", -nov", neg_nov_com, 
                  ", -rep", neg_rep_com, ", Base Rate = ", baserate, "alpha:", alpha)) + theme_bw()
  
  plot(gg)
  } } }

#in a column called tally, stores the number of observations for each hypothesis
n_studies_per_id <- aggregate(tally ~ id + true, data=r_df_comm, FUN = function(x) length(x))
final_df$n_tests <- n_studies_per_id$tally

#unique tally values, and limit to between the tally bound for all tallies in final_df
tally_bound <- c(-10, 10)
u.tally <- unique(final_df$tally)

#limiting tallies
final_df$tally[final_df$tally < tally_bound[1]] <- tally_bound[1]
final_df$tally[final_df$tally > tally_bound[2]] <- tally_bound[2]
u.tally <- u.tally[u.tally >= tally_bound[1] & u.tally <= tally_bound[2]]
u.tally <- sort(u.tally)

#for T/F hypotheses, for each unique tally value, creating a separate data frame and storing it in a list
TH_tally_list <- list()
FH_tally_list <- list()

#False Hypotheses
for(i in u.tally){
  df_1_tally <- final_df[final_df$true==0 & final_df$tally==i,]
  
  if(nrow(df_1_tally) > 0) {
    
  simplehist(df_1_tally$n_tests, ylab = "Frequency", xlab = "Number of Hypothesis Tests", 
             main = paste("FH., Tally:",i, "+nov:", pos_nov_com, ", +rep", PR_COM, ", -nov", neg_nov_com, 
                          ", -rep", neg_rep_com, ", Base Rate = ", baserate, "alpha:", alpha))
          
  
  df_1_tally$pos_nov_com <- pos_nov_com
  df_1_tally$neg_nov_com <- neg_nov_com
  df_1_tally$pos_rep_com <- PR_COM
  df_1_tally$neg_nov_com <- neg_rep_com
    
  FH_tally_list[[length(FH_tally_list)+1]] <- df_1_tally }
}

#True Hypotheses
for(i in u.tally){
  df_1_tally <- final_df[final_df$true==1 & final_df$tally==i,]
  
  if(nrow(df_1_tally) > 0) {
  
  simplehist(df_1_tally$n_tests, ylab = "Frequency", xlab = "Number of Hypothesis Tests", 
             main = paste("TH., Tally:",i, "+nov:", pos_nov_com, ", +rep", PR_COM, ", -nov", neg_nov_com, 
                          ", -rep", neg_rep_com, ", Base Rate = ", baserate, "alpha:", alpha))
  
  df_1_tally$pos_nov_com <- pos_nov_com
  df_1_tally$neg_nov_com <- neg_nov_com
  df_1_tally$pos_rep_com <- PR_COM
  df_1_tally$neg_nov_com <- neg_rep_com
  
  TH_tally_list[[length(TH_tally_list)+1]] <- df_1_tally }
}

#For each element in FH_tally_list and TH_tally_list, we can see a table for the number of times that x tests 
# led to that tally as follows

unique(TH_tally_list[[5]]$tally) #the tally
table(TH_tally_list[[5]]$n_tests) #the distribution of tests




