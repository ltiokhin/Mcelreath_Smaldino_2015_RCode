library(rethinking)

rm(list=ls())  

#parameters for the world
beta <- 0.2
power <- 1 - beta
alpha <- c(0.05)
trueneg <- 1 - alpha
baserates <- c(0.01)
p_novel <- 0.8

n <- 40 #hypotheses tested each round
rounds <- 1000
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
  
  #empty list to store the frequency of the number of hypothesis tests to reached a given tally
  n_test_freq_list <- list()
  
  #FOR MANY RUNS
  for(runs in 1:rounds){
    
    #if no hypotheses to be replicated, then samples only novel hypotheses
    if(is.na(r_df$true[1])) {
      q_df$true <- novel_h[1:n]
      #remove already sampled hypotheses
      novel_h <- novel_h[-c(1:n)]
      #update novel column in q_df
      q_df$novel <- 1 
    } else{
      
      for(i in 1:n) {
        if(runif(1, 0, 1) < p_novel){
          q_df$true[i] <- novel_h[1]
          q_df$novel[i] <- 1 
          #remove already sampled hypothesis
          novel_h <- novel_h[-1]
        } else {
          #generates position from those positions that contain tested hypotheses
          position <- sample(which(!is.na(r_df$id)), 1, replace = TRUE)
          #adds previously tested hypotheses to the q_df data frame
          q_df$id[i] <- r_df$id[position]
          q_df$true[i] <- r_df$true[position]
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
  
  ##aggregate data frame for each unique hypothesis id
  final_df_t <- r_df[r_df$true==1,]
  final_df_t <- aggregate(tally ~ id, data=final_df_t, FUN=sum)
  final_df_f <- r_df[r_df$true==0,]
  final_df_f <- aggregate(tally ~ id, data=final_df_f, FUN=sum)
  final_df <- rbind(final_df_t, final_df_f)
  
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

#addd a column to final_df that represents how many times each hypothesis was tested ####kaka
length(unique(r_df$id))
length(unique(final_df$id))


final_df$n_tests <- NA
n_studies_per_id <- as.data.frame(table(r_df$id))
z <- cbind(final_df$id, as.numeric(n_studies_per_id$Var1))
options("max.print" = 2000)
tail(z)

for(i in n_studies_per_id$Var1){
  final_df$n_tests[final_df$id==i] <- n_studies_per_id$Var1[i]
}


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

final_df[final_df$id==6132,]

#6132 appears 5 times - -what is its tally

