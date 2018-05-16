library(rethinking)

rm(list=ls())  

#parameters for the world
beta <- 0.2
power <- 1 - beta
alpha <- 0.05
trueneg <- 1 - alpha
baserates <- c(0.001, 0.01, 0.1, 0.5)
p_novel <- 0.8
#novel hypotheses 
n <- 20 #hypotheses tested each round

#FOR DIFFERENT BASERATES
for(baserate in baserates) {
  
novel_h <- rbinom(1e5, 1, baserate)
                    
#tested hypotheses
r_df <- data.frame(id = rep(NA, 1e5),
                  true = rep(NA, 1e5),
                  result = rep(NA, 1e5), 
                  comm = rep(NA, 1e5), 
                  tally = rep(0, 1e5))

#communication parameters
pos_nov_com <- 1
pos_rep_com <- 0
neg_nov_com <- 0
neg_rep_com <- 1

#data frame to store what happens each round
q_df <- data.frame(id = rep(NA, n),
                  true = rep(NA, n),
                  novel = rep(NA, n), 
                  result = rep(NA, n), 
                  comm = rep(NA, n), 
                  tally = rep(0, n))

#FOR MANY RUNS
for(runs in 1:1000){
  
#GENERATE HYPOTHESES TO BE TESTED
#If there are no hypotheses to be replicated, then all hypotheses are novel. If there is at least 1 hypothesis to be replicated, 
# then, for each of n possible hypotheses, we check whether an runif number is less than p_novel. If so, then hypothesis tested is novel. if not, then it's a replication. 

  #need to assign ids

#if no hypotheses to be replicated, then samples only novel hypotheses
if(is.na(r_df$true[1])) {
  q_df$true <- novel_h[1:n]
  #remove already sampled hypotheses
  novel_h <- novel_h[-1:n]
  #update novel column in q_df
  q_df$novel <- 1 
} else{
  
  for(i in 1:n) {
    if(runif(1, 0, 1) < p_novel){
      q_df$true[i] <- novel_h[1]
      q_df$novel[i] <- 1 
      #remove the rows associated with tested hypotheses from novel_h
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
q_df <- within(q_df, comm[novel==0 & result == 1] <- rbinom(length(comm[novel==0 & result == 1]), 1, pos_rep_com))
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

} #end loop

##aggregate data frame for each unique hypothesis id
final_df <- aggregate(tally ~ id + true, data=r_df, FUN=sum)
simplehist(final_df$tally[final_df$true==0], ylab = "Number of Hypotheses", xlab = "Tally", 
           main = paste("False Hypotheses, Base Rate = ", baserate))
simplehist(final_df$tally[final_df$true==1], ylab = "Number of Hypotheses", xlab = "Tally", 
           main = paste("True Hypotheses, Base Rate = ", baserate))

##precision: proportion of hypotheses at each tally that are true
precision_df <- aggregate(true ~ tally, data=final_df, FUN = function(x) sum(x) / length(x))

gg <- ggplot(precision_df, aes(x = tally)) + geom_histogram(aes(y = true), stat = "identity") + 
  scale_y_continuous(name = "Proportion of True Hypotheses at Tally") + 
  ggtitle(paste("Only + Novel Results and - Replications Comm; Base Rate =", baserate)) +
  theme_bw()

plot(gg)
}




