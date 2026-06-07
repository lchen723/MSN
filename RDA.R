library(MASS)  
library(stats)
library(survival)
library(glasso)
library(igraph)
library(Matrix)
library(ggplot2)
library(dplyr)
library(viridis)
library(lattice)

## data setup
Xstar[[9]] <- Xstar0
Y[[9]] <- Y0
###
X = NULL
for(k in 1:9) {
X[[k]] = Xstar[[k]][,1:200]
}
###

event_time = NULL
delta = NULL
for(k in 1:9) {
event_time[[k]] = Y[[k]][,1] + 0.1
delta[[k]] = (Y[[k]][,1] == Y[[k]][,2])*1
}

event_names_9 <- c(
  "GES58812", "GES48390", "GES42568", "GES31448", 
  "GES21653", "GES20711", "GES20685", "GES16446", 
  "GSE88770"
)
gene_names = colnames(X[[1]])
## Analysis

est_beta = NULL 
for(k in 1:9) {
est_beta[[k]] = boosting_beta(event_time[[k]], delta[[k]], X[[k]])$beta
}

est_sur = NULL; est_b = NULL
for(k in 1:9) {
semi_est = semi_estimation(y_vec = seq(0,max(event_time[[k]]),length=100), X_target = X[[k]], 
     Y_train = event_time[[k]], Delta_train = delta[[k]], X_train = X[[k]],kappa=0) 

est_sur[[k]] = semi_est$est_F
est_b[[k]] = semi_est$est_beta
}

seq_network = NULL
se = seq(0,1,length=50); se = se[-1]
for(k in 1:length(se)) {
seq_network[[k]] = network_estimate(Y = event_time, delta = delta, X = X, rho = se[k],kappa=0)

}

est_network = network_estimate_BIC(Y = event_time, delta = delta, X = X, rho_seq = seq(0,0.8,length=50),kappa=0) 

### draw graphs
  net = network::network(est_network$pre_matrix, directed = FALSE)
  network::network.vertex.names(net)=paste0("Y",network::network.vertex.names(net))
  GGally::ggnet2(net,size=10,node.color = "lightgray",label=event_names_9,label.size = 2.7,mode = "circle")+ 
  labs(title = "Scenario II-0.8")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

### summarize beta

Beta = NULL
for(k in 1:9) {
Beta = cbind(Beta, est_b[[k]])
}
freq = c()
for(j in 1:dim(Beta)[1]) {
freq = c(freq, length(which(Beta[j,]!=0)))
}
gene_summary = data.frame(
gene_level = factor(gene_names[which(freq>=2)], levels = gene_names[which(freq>=2)]),
freq_count = freq[which(freq>=2)]
)
ggplot(gene_summary, aes(x = gene_level, y = freq_count, fill = gene_level)) +
  geom_bar(stat = "identity", color = "black", width = 0.5) + 
  geom_text(aes(label = gene_level), vjust = -0.5, fontface = "bold") +
  theme_minimal(base_size = 14) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(title = "X-Axis Rotated Bar Chart", x = "", y = "Count Number") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))



#####################

## merge data
ID = which(est_b[[4]]!=0)
event_time_new = c(event_time[[4]], event_time[[5]])
delta_new = c(delta[[4]],delta[[5]])
X_new = rbind(X[[4]],X[[5]])
X_new = X_new[,ID]

TL = semi_estimation(y_vec = seq(0,20,length=100), X_target = X_new, 
     Y_train = event_time_new, Delta_train = delta_new, X_train = X_new,kappa=0)

GES31448 = semi_estimation(y_vec = seq(0,20,length=100), X_target = X[[4]], 
           Y_train = event_time[[4]], Delta_train = delta[[4]], X_train = X[[4]],kappa=0)

GES21653 = semi_estimation(y_vec = seq(0,20,length=100), X_target = X[[5]], 
           Y_train = event_time[[5]], Delta_train = delta[[5]], X_train = X[[5]],kappa=0)

plot(seq(0,20,length=100), 1-TL$est_F[,250],type="l", lwd=3,xlab="t", ylab="S")
points(seq(0,20,length=100), 1-GES31448$est_F[,250],col=2,type="l", lwd=3,lty=2)
points(seq(0,20,length=100), 1-GES21653$est_F[,250],col=4,type="l", lwd=3,lty=3)









