##### Code to analyses the simulations & plots op characteristics of models
##### Requires inputs from simulations.R 

# Function to extract op characteristics from a list of simulations
# sims = list of simulations
op_char <- function(sims, alpha=0.05) {
  
  # initialise output
  input         <- array(0, dim=c(length(sims),7))
  results       <- array(0, dim=c(length(sims),6))
  results.early <- array(0, dim=c(length(sims),6))
  results.late  <- array(0, dim=c(length(sims),6))
  
  for (i in 1:length(sims)) {
    
    ##### input parameters of the simulations
    if( length(sims[[i]]$r)<2 ){
      input[i,1] <- input[i,2] <- sims[[i]]$r
    }else{
      input[i,1] <- sims[[i]]$r[1]
      input[i,2] <- sims[[i]]$r[2]
    }
    input[i,3] <- sims[[i]]$nSim
    input[i,4] <- sims[[i]]$N
    input[i,5] <- sims[[i]]$tau
    input[i,6] <- sims[[i]]$tend
    input[i,7] <- mean(sims[[i]]$cox.t[, "convergence"], na.rm = TRUE)
    
    # true HR r^(-1/relapse_scale)
    delta.e <- input[i,1]^(-1/relapse_scale)
    delta.l <- input[i,2]^(-1/relapse_scale)
    
    # true HR overall
    lambda.0 = exp(relapse_coefficients)
    lambda.1 = lambda.0*input[i,1]
    S0 <- pweibull( input[i,5], shape=relapse_scale, scale=lambda.0, lower.tail=FALSE )
    S1 <- pweibull( input[i,5], shape=relapse_scale, scale=lambda.1, lower.tail=FALSE )
    omega_bar <- 1/2*(S0+S1)
    
    delta <- exp((1-omega_bar)*log(delta.e)+omega_bar*log(delta.l))
    
    
    ##### Cox univariate model
    # mean of estimate (HR)
    results[i,1] <- mean(sims[[i]]$cox[,"HR"], na.rm = TRUE)
    # sd of estimate (HR)
    results[i,2] <- sd(sims[[i]]$cox[,"HR"], na.rm = TRUE)
    # mse
    results[i,3] <- sqrt(mean((sims[[i]]$cox[,"HR"]-delta)^2))
    # bias
    results[i,4] <- mean(sims[[i]]$cox[,"HR"], na.rm = TRUE)-delta
    # power
    results[i,5] <- mean(sims[[i]]$cox[,"p-value"]<alpha, na.rm = TRUE)
    # coverage
    results[i,6] <- 100*mean((sims[[i]]$cox[,"lower CI"]<delta) 
                             & (sims[[i]]$cox[,"higher CI"]>delta), na.rm = TRUE)
    
    
    ##### Cox model with period
    # mean of estimate (HR.early)
    results.early[i,1] <- mean(sims[[i]]$cox.t[,"HR.early"], na.rm = TRUE)
    # sd of estimate (HR)
    results.early[i,2] <- sd(sims[[i]]$cox.t[,"HR.early"], na.rm = TRUE)
    # mse
    results.early[i,3] <- sqrt(mean((sims[[i]]$cox.t[,"HR.early"]-delta.e)^2))
    # bias
    results.early[i,4] <- mean(sims[[i]]$cox.t[,"HR.early"], na.rm = TRUE)-delta.e
    # power
    results.early[i,5] <- mean(sims[[i]]$cox.t[,"p-value.early"]<alpha, na.rm = TRUE)
    # coverage
    results.early[i,6] <- 100*mean((sims[[i]]$cox.t[,"lower CI.early"]<delta.e) 
                                   & (sims[[i]]$cox.t[,"higher CI.early"]>delta.e), na.rm = TRUE)
    
    
    # mean of estimate (HR.late)
    results.late[i,1] <- mean(sims[[i]]$cox.t[,"HR.late"], na.rm = TRUE)
    # sd of estimate (HR)
    results.late[i,2] <- sd(sims[[i]]$cox.t[,"HR.late"], na.rm = TRUE)
    # mse
    results.late[i,3] <- sqrt(mean((sims[[i]]$cox.t[,"HR.late"]-delta.l)^2))
    # bias
    results.late[i,4] <- mean(sims[[i]]$cox.t[,"HR.late"], na.rm = TRUE)-delta.l
    # power
    results.late[i,5] <- mean(sims[[i]]$cox.t[,"p-value.late"]<alpha, na.rm = TRUE)
    # coverage
    results.late[i,6] <- 100*mean((sims[[i]]$cox.t[,"lower CI.late"]<delta.l) 
                                  & (sims[[i]]$cox.t[,"higher CI.late"]>delta.l), na.rm = TRUE)
    
  }
  
  # change arrays to dataframe & rename columns
    # bind three array together
  output <- cbind(input, results, results.early, results.late)
  output <- as.data.frame(output)
  
  #colnames
  colnames(output) <- c("r.e","r.l","nSim","sample","tau","tend","Convergence",
                        "Mean","SD","MSE","Bias","Power","Coverage",
                        "Mean.early","SD.early","MSE.early","Bias.early","Power.early","Coverage.early",
                        "Mean.late","SD.late","MSE.late","Bias.late","Power.late","Coverage.late")
  
  # round results to 3 decimals
  output <- apply(output, 2, round, digits=3)
  output <- as.data.frame(output)
  
  return(output)
  
}



# Function to plot operational characteristics of a set of simulations
# filename = name of the .RData file containing the saved simulations
# NOTE: files are assumed to be located in the working directory 
#       figures are saved to folder /Figures located inside the working directory
analyse_results <- function(filename=NULL){
  
  # load simulations
  load(file=filename)
  label <- gsub(".RData", "_", filename)
  
  # calculate the operationa characteristics from a set of simulations
  results <- as.data.frame(op_char(sims))
  
  # re-organise results, binding "HR overall", "HR early" and "HR late" vertically
  tmp <- cbind(HR="HR overall", 
               results[,c("r.e","r.l","nSim","sample","tau","tend","Convergence")], 
               results[,c("Mean","SD","MSE","Bias","Power","Coverage")])
  
  tmp.early <- cbind(HR="HR early", 
                     results[,c("r.e","r.l","nSim","sample","tau","tend","Convergence")], 
                     results[,grep("early", names(results))])
  
  tmp.late  <- cbind(HR="HR late", 
                     results[,c("r.e","r.l","nSim","sample","tau","tend","Convergence")],
                     results[,grep("late", names(results))])
  
  names(tmp.early) <- names(tmp)
  names(tmp.late) <- names(tmp)
  
  table <- rbind(tmp, tmp.early, tmp.late)
  
  # define delta (true HR) from values of r
  table$delta <- paste0(round(table$r.e^(-1/relapse_scale), digits=2), " - ",
                        round(table$r.l^(-1/relapse_scale), digits=2))
  
  # plot Mean & SE (with preview)
  fig <- ggplot(data=table, 
                aes(x=tend, y=Mean, group=delta, colour=delta)) + 
    geom_errorbar(aes(ymin=Mean-SD/sqrt(nSim), ymax=Mean+SD/sqrt(nSim)), width=.1) +
    ggtitle("Mean & SE") +
    geom_line()+geom_point() + 
    facet_grid(sample+nSim~HR, scales="free_x", labeller=label_both)
  ggsave(filename = paste0(PATH,"/Figures/",label,"Mean.png"),
         width = 20, height = 14, units = "cm")
  print(fig)
  
  # re-organise table in long form
  table.long <- melt(table, id=c("HR","r.e","r.l","nSim","sample","tau","tend","delta"))
  
  # function to plot figures of op characteristics (figures are saved)
  op_fig <- function(scale_char, data=table.long){
    data %>% 
      filter(variable==scale_char) %>%
      ggplot(aes(x=tend, y=value, colour=delta)) +
      ggtitle(scale_char) +
      geom_line()+geom_point() + 
      facet_grid(sample+nSim~HR, scales="free_x", labeller=label_both)
    ggsave(filename = paste0(PATH,"/Figures/",label,scale_char,".png"),
           width = 20, height = 14, units = "cm")
  }
  
  # plot each of the op characteristics
  op_fig("MSE")
  op_fig("Bias")
  op_fig("Power")
  op_fig("Coverage")
  op_fig("Convergence")
  
  return(table)
  
} 


# analyse simulations
table <- analyse_results("sims_1.RData")
table <- analyse_results("sims_2.RData")
table <- analyse_results("sims_3.RData")
table <- analyse_results("sims_4.RData")
table <- analyse_results("sims_5.RData")
table <- analyse_results("sims_6.RData")



####### Model 1: simulate the null hypothesis #######
table <- analyse_results("sims_1.RData")

# Use only sample = 167
table.1 <- melt(table, id = c("HR","r.e","r.l","nSim","sample","tau","tend","delta"))

# plot MSE, Bias, Power, Coverage
ggplot(data=table.1[table.1$variable %in% c("MSE","Bias","Power","Coverage"), ],
       aes(x=tend, y=value, group=factor(sample), colour=factor(sample))) +
  geom_line()+geom_point() +
  geom_vline(xintercept = 20, linetype="dashed") +
  xlab("Analysis change point") + ylab("") +
  scale_color_discrete("sample") +
  facet_grid(variable~HR, scales="free_y") + 
  theme_bw(base_size = 20)
ggsave(filename = paste0(PATH,"/Figures/fig1.png"), width = 24, height = 18, units = "cm")




####### Model 2: constant overall treatment effect #######
table <- analyse_results("sims_2.RData")

# Use only sample = 167
table.2 <- melt(table[table$sample==167, ], id = c("HR","r.e","r.l","nSim","sample","tau","tend","delta"))

# plot MSE, Bias, Power, Coverage
ggplot(data=table.2[table.2$variable %in% c("MSE","Bias","Power","Coverage"), ],
       aes(x=tend, y=value, group=delta, colour=delta)) +
  geom_line()+geom_point() +
  geom_vline(xintercept = 20, linetype="dashed") +
  xlab("Analysis change point") + ylab("") +
  scale_color_discrete(expression(HR[e] - HR[l])) +
  facet_grid(variable~HR, scales="free_y") + 
  theme_bw(base_size = 20)
ggsave(filename = paste0(PATH,"/Figures/fig2.png"), width = 24, height = 18, units = "cm")




####### Model 3: late efficacy (HRearly = 1), treatment effect only after month 20) ######
table <- analyse_results("sims_3.RData")

# Use only sample = 167
table.3 <- melt(table[table$sample==167, ], id = c("HR","r.e","r.l","nSim","sample","tau","tend","delta"))

# plot MSE, Bias, Power, Coverage
ggplot(data=table.3[table.3$variable %in% c("MSE","Bias","Power","Coverage"), ],
       aes(x=tend, y=value, group=delta, colour=delta)) +
  geom_line()+geom_point() +
  geom_vline(xintercept = 20, linetype="dashed") +
  xlab("Analysis change point") + ylab("") +
  scale_color_discrete(expression(HR[e] - HR[l])) +
  facet_grid(variable~HR, scales="free_y") + 
  theme_bw(base_size = 20)
ggsave(filename = paste0(PATH,"/Figures/fig3.png"), width = 24, height = 18, units = "cm")




####### Model 4: early efficacy (HRlate = 1), treatment effect up until month 20) ######
table <- analyse_results("sims_4.RData")

# Use only sample = 167
table.4 <- melt(table[table$sample==167, ], id = c("HR","r.e","r.l","nSim","sample","tau","tend","delta"))

# plot MSE, Bias, Power, Coverage
ggplot(data=table.4[table.4$variable %in% c("MSE","Bias","Power","Coverage"), ],
       aes(x=tend, y=value, group=delta, colour=delta)) +
  geom_line()+geom_point() +
  geom_vline(xintercept = 20, linetype="dashed") +
  xlab("Analysis change point") + ylab("") +
  scale_color_discrete(expression(HR[e] - HR[l])) +
  facet_grid(variable~HR, scales="free_y") + 
  theme_bw(base_size = 20)
ggsave(filename = paste0(PATH,"/Figures/fig4.png"), width = 24, height = 18, units = "cm")




####### ####### Model 5: different treatment effects before/after month 20  ######
table <- analyse_results("sims_5.RData")

# Use only sample = 167
table.5 <- melt(table[table$sample==167, ], id = c("HR","r.e","r.l","nSim","sample","tau","tend","delta"))

# plot MSE, Bias, Power, Coverage
ggplot(data=table.5[table.5$variable %in% c("MSE","Bias","Power","Coverage"), ],
       aes(x=tend, y=value, group=delta, colour=delta)) +
  geom_line()+geom_point() +
  geom_vline(xintercept = 20, linetype="dashed") +
  xlab("Analysis change point") + ylab("") +
  scale_color_discrete(expression(HR[e] - HR[l])) +
  facet_grid(variable~HR, scales="free_y") + 
  theme_bw(base_size = 20)
ggsave(filename = paste0(PATH,"/Figures/fig5.png"), width = 24, height = 18, units = "cm")




####### Model 6: different treatment effects, wrong change point tau ######
table <- analyse_results("sims_6.RData")

# Use only sample = 167
table.6 <- melt(table[table$sample==167, ], id = c("HR","r.e","r.l","nSim","sample","tau","tend","delta"))

# plot MSE, Bias, Power, Coverage
ggplot(data=table.6[table.6$variable %in% c("MSE","Bias","Power","Coverage"), ],
       aes(x=tau, y=value, group=delta, colour=delta)) +
  geom_line()+geom_point() +
  geom_vline(xintercept = 20, linetype="dashed") +
  xlab("True change point") + ylab("") +
  scale_color_discrete(expression(HR[e] - HR[l])) +
  facet_grid(variable~HR, scales="free_y") + 
  theme_bw(base_size = 20)
ggsave(filename = paste0(PATH,"/Figures/fig6.png"), width = 24, height = 18, units = "cm")




##################################################
# save all data
save.image(file="save_results_Oct03.Rdata")



