# input is the combined result table from runSVM
plotSVM <- function(res) {
  # aggregate data
  dat.plot <- res %>% group_by(run_name, train_prop) %>% summarise_all(funs(mean, sd))
  # make plot
  p1 <- ggplot(res_plot, aes(x=props, y=Avg, group = gp)) + 
    geom_errorbar(aes(ymin=Avg-Sd, ymax=Avg+Sd), width=.02) +
    geom_line( aes(x=props, y=Avg, linetype = method)) + 
    geom_point(size=4) + 
    xlab("Traning Data Proportion") +
    ylab("Prediction Accuracy") +
    theme(legend.position = "none") +
    ggtitle("A")
  
}