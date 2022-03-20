# visualization function - weights of assets according to the lambda sequences
# real_data
real_viz <- function(ww, eta){
  plot(ww[1,], ylim=c(min(ww),max(ww)),type='l', lwd=1, main=paste("data = ", data_set[number], " eta = ",eta))
  for(i in 2:nrow(ww)) lines(ww[i,], type='l') #, col=ddd[i])
}

# simulation_data
simul_viz <- function(ww, eta) {
  plot(ww[1,], ylim=c(min(ww),max(ww)), type='l', col=1, lwd=2, main=paste("data = simulation, eta = ", eta))
  for(i in 2:4) lines(ww[i,], type='l', col=1, lwd=2)
  for(i in 5:8) lines(ww[i,], type='l', col=2, lwd=2)
  for(i in 9:12) lines(ww[i,], type='l', col=3, lwd=2)
} 
