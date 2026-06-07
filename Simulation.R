#==========library==========
library(MASS)  
library(stats)
library(survival)
library(glasso)
library(igraph)
library(Matrix)
library(energy)
library(ggplot2)
library(pROC)
library(dplyr)
library(viridis)
library(tidyverse)
library(PAGE)
#========== p-dimention ==========
set.seed(168)
simGraph <-function(p=50, type="scale-free"){
  if(type=="scale-free"){
    g <- igraph::barabasi.game(n=p, power=0.01, zero.appeal=p, directed=F)
    theta <- as.matrix(igraph::get.adjacency(g, type="both"))
  }
  #
  if(type=="hub"){
    if(p > 40) {g = ceiling(p/20)}
    if(p <= 40) {g = 2} 
    
    g.list = c(rep(floor(p/g), (g - p%%g)), rep(floor(p/g)+1, p%%g))
    g.ind = rep(c(1:g), g.list)
    u = 0.1
    v = 0.3
    theta = matrix(0, p, p)
    for (i in 1:g) {
      tmp = which(g.ind == i)
      theta[tmp[1], tmp] = 1
      theta[tmp, tmp[1]] = 1
    }
  }
  #
  if(type=="lattice"){
    
    latticeStructure <- function(p){
      nnodes <- sqrt(p)
      if((nnodes %% 1) == 0){
        rnodes = nnodes
        cnodes = nnodes
      }
      if((nnodes %% 1) != 0){
        cnodes = 3
        rnodes = p / cnodes
        while((rnodes %% 1) != 0){
          cnodes = cnodes + 1
          rnodes = p / cnodes
        }
      }
      return(c(rnodes, cnodes))
    }
    
    temp = latticeStructure(p)
    rnodes = temp[1]
    cnodes = temp[2]
    
    theta=diag(1,nrow=rnodes*cnodes,ncol=cnodes*rnodes)
    alterSign = FALSE
    for(i in 1:(rnodes*cnodes)){
      #print(i);
      if( i%% rnodes !=0 ){
        theta[i,i+1]=1
      }
      if((i-1)%%rnodes!=0){
        theta[i,i-1]=1
      }
      if((i-rnodes)>0){
        theta[i,i-rnodes]=1
        if(alterSign) 	theta[i,i-rnodes] = -1
        
      }			
      if((i+rnodes)<=rnodes*cnodes){
        theta[i,i+rnodes]=1
        if(alterSign) 	theta[i,i+rnodes] = -1
        
      }			
    }
  }
  diag(theta) <- 0
  return(theta)
}
XMRF.Sim <-function(n=100, p=50, model="LPGM", graph.type="hub"){
  mydata <- list()
  
  # Generate a network first
  B <- simGraph(p=p, type=graph.type)
  mydata$B <- B
  
  if(model == "LPGM" || model == "SPGM" || model == "TPGM"){
    ## implement the speed.pois="fast" in mp.generator
    lambda = 2
    lambda.c = 0.5
    
    A = adj2A(B,type=graph.type)
    sigma=lambda*B
    ltri.sigma=sigma[lower.tri(sigma)]
    nonzero.sigma = ltri.sigma[which(ltri.sigma !=0 )]
    Y.lambda = c(rep(lambda,nrow(sigma)), nonzero.sigma)
    
    Y = rmpois(n,Y.lambda)
    X =A%*%Y
    # add the labmda.c to all the nodes. 
    X = X + rmpois(ncol(X), rep(lambda.c,nrow(X)))
    mydata$X=X
    ##mydata$A=A
    ##mydata$sigma=sigma
  }
  #
  if(model == "GGM"){
    # Simmuate a Gaussian multivariate data matrix
    u = 0.1
    v = 0.3
    diag(B) = 0
    omega = B * v
    diag(omega) = abs(min(eigen(omega)$values)) + 0.1 + u
    sigma = cov2cor(solve(omega))
    omega = solve(sigma)
    X = MASS::mvrnorm(n, rep(0, p), sigma)
    
    mydata$X <- t(X)
    ##mydata$A=NULL
    mydata$sigma=sigma
  }
  #
  if(model == "ISM"){
    theta = 0.4
    maxit = 1000
    
    X <- matrix(rbinom(n*p, 1, 0.6), nrow=n, ncol=p)
    X[X==0] <- -1
    # Changing sign of the parts of adjacency matrix 
    ####dims = 1:size(x,2);Ising0.
    Theta <- -theta*B
    Theta[,1:(p/2)] <- -Theta[,1:(p/2)]
    Iterat <- array(0, dim=c(n,p,maxit))
    t = 1
    # RUN GIBBS SAMPLER
    while (t < maxit){
      t = t+1
      # Loop over dimensions  
      for (iD in 1:p){
        odds <- exp(2*X[, -iD]%*%Theta[-iD, iD])
        prob <- odds /(1+odds)
        X[, iD] <- rbinom(length(prob), 1, prob)*2 -1  
      }
      Iterat[,,t] <- X
    }
    xDat = (X+1)/2
    mydata$X <- t(xDat)
    ##mydata$sigma=NULL
  }
  #
  if(model == "PGM"){
    mu = 0.5
    sigma = sqrt(0.01)
    
    mu = 0.1
    signam = 0.01
    
    Theta = matrix(rnorm(p*p,mu,sigma),nrow = p)
    Theta = Theta * B
    X = PGMSim(n,p, alpha = rep(0,p), Theta=Theta, maxit=1000)
    X = PGMSim(n,p, alpha = rep(0,p), Theta=Theta, maxit=130)
    
    mydata$X=X
    ##mydata$A=NULL
    ##mydata$sigma=NULL
  }
  
  return(mydata)
}
p <- 15
n <- 100
event<-15
graph.type <- "scale-free"  
sim_result <- XMRF.Sim(n = n, p = event , model = "GGM", graph.type = graph.type)

X <- mvrnorm(n, rep(0, p), diag(1, p))
P <- sim_result$B * 0.5         
diag(P) <- 1.5
Sigma <- cov2cor(solve(P))

Z <- mvrnorm(n, mu = rep(0, event), Sigma=Sigma)
U <- pnorm(Z)
beta_mat <- matrix(0, nrow = ncol(X), ncol = event)  # 18 × 15

for (j in 1:event) {
  beta_mat[sample(1:ncol(X), 3), j] <- 1
}

event_types <- c(rep("COX", 7), rep("AFT", 8)) 
simulate_cox_time <- function(U, linpred, type = "exp", a = 1, b = 1, c = 1) {
  if (type == "exp") {
    # λ₀(t) = 1 → exponential
    return(-log(1 - U) / exp(linpred))
  } else if (type == "weibull") {
    # λ₀(t) = abt^{b-1}
    return((-log(1 - U) / (a * exp(linpred)))^(1 / b))
  } else if (type == "linear") {
    # λ₀(t) = ct 
    return(sqrt(-2 * log(1 - U) / (c * exp(linpred)))
    )
  } else {
    stop("Unsupported baseline hazard type.")
  }
}
cox_types <- c("exp", "exp", "exp", "weibull", "weibull", "linear", "linear")
t_mat <- matrix(NA, nrow = n, ncol = event)
for (j in 1:event) {
  linpred <- X %*% beta_mat[, j]
  if (event_types[j] == "COX") {
    t_mat[, j] <- simulate_cox_time(U[, j], linpred, type = cox_types[j])
  } else {
    # AFT 模型
    t_mat[, j] <- exp(linpred + qnorm(U[, j]))
  }
}#生成t_ij
c_mat <- matrix(rexp(n * event, rate = 0.8), nrow = n)
y_mat <- pmin(t_mat, c_mat)
delta_mat <- (t_mat <= c_mat) * 1
colnames(t_mat) <- paste0("Event", 1:ncol(t_mat))

a<-as.data.frame(colMeans(delta_mat))
colMeans(a)

data_list <- lapply(1:event, function(j) {
  list(
    y = log(y_mat[, j]),
    delta = delta_mat[, j],
    label = paste0("Event", j),
    model = event_types[j]
  )
})

#==========5 independent events
set.seed(168)
n <- 100
p_indep <- 5

X_list <- lapply(1:p_indep, function(j) matrix(rnorm(n * p_indep), n, p_indep))
beta_list <- lapply(1:p_indep, function(j) {
  beta <- rep(0, p_indep)  # 預設全部為 0
  important_idx <- sample(1:p_indep, 5)  
  beta[important_idx] <- 1 
  return(beta)
})


T_indep <- sapply(1:p_indep, function(j) {
  Xj <- X_list[[j]]
  betaj <- beta_list[[j]]
  linpred <- Xj %*% betaj
  base_hazard <- rexp(n, rate = 1) 
  Tj <- base_hazard * exp(linpred)
  return(Tj)
})

colnames(T_indep) <- paste0("T", 16:(15+p_indep))

censor_indep <- rexp(n, rate = 0.2)
y_indep <- pmin(T_indep, censor_indep)
delta_indep <- (T_indep <= censor_indep) * 1



# 合併事件資料
t_all <- cbind(t_mat, T_indep)
delta_all <- cbind(delta_mat, delta_indep)
event_names <- c(colnames(t_mat), colnames(T_indep))
y_all<-cbind(y_mat,y_indep)
colnames(t_all) <- event_names
colnames(delta_all) <- event_names
covariate_map <- c(rep("X_all", 15), rep("Z", p_indep))
names(covariate_map) <- event_names

#========== Boosting ==========
boosting_beta <- function(y, delta, X, max_iter = 400,
                          learning.rate = 0.1, stop_value = 0.001) {
  n <- length(y)
  p <- ncol(X)
  beta <- rep(0, p)
  
  # Kaplan–Meier 估 G(y) = P(C > y)；事件設為 (1 - delta)
  fit_G <- survfit(survival::Surv(y, 1 - delta) ~ 1)
  G_fun <- stepfun(fit_G$time, c(1, fit_G$surv))
  G_y <- G_fun(y)
  G_y[G_y <= 1e-8] <- 1e-8
  w <- delta / G_y  # ICPW 權重
  loss <- function(beta) {
    r <- as.vector(y - X %*% beta)
    sum(w * r^2)/n
  }
  
  loss_vec <- numeric(max_iter)
  for (iter in 1:max_iter) {
    r <- as.vector(y - X %*% beta)
    grad <- as.vector(-2 * t(X) %*% (w * r) / n)
   
    j <- which.max(abs(grad))
    xj <- X[, j]
    num <- sum(w * r * xj)
    den <- sum(w * xj^2)
    if (den <= 1e-12) {
      remaining <- setdiff(1:p, j)
      if (length(remaining) == 0) {
        break  
      } 
      else {
        next
      }
    }
    s_star <- num / den
    step <- learning.rate * s_star
    beta_old <- beta
    beta[j] <- beta[j] + step
    
    loss_vec[iter] <- loss(beta)
    
    if (max(abs(beta - beta_old)) < stop_value) {
      cat("Converged at iteration", iter, "\n")
      loss_vec <- loss_vec[seq_len(iter)]
      break
    }
  }
  
  list(beta = beta, loss = loss_vec)
}

beta_hat_list <- lapply(data_list, function(d) {
  boosting_beta(d$y, d$delta, X)
})
names(beta_hat_list) <- sapply(data_list, function(d) d$label)
extra_event_names <- colnames(T_indep)
extra_beta_list <- lapply(seq_along(extra_event_names), function(j) {
  boosting_beta(
    y = log(y_indep[, j]),
    delta = delta_indep[, j],
    X = X_list[[j]]
  )
})
names(extra_beta_list) <- extra_event_names

beta_hat_list_all <- c(beta_hat_list, extra_beta_list)


#========== print true & estimated beta ==========
true_betas <- vector("list",event)
for (j in 1:event) {
  true_betas[[j]] <- beta_mat[, j]
}
names(true_betas) <- paste0("Event", 1:event)

for (j in 1:event) {
  cat("\n====", names(beta_hat_list)[j], "====\n")
  print(cbind(True = true_betas[[j]], Estimated = round(beta_hat_list[[j]]$beta,3)))
}#print ture & estimated
for (j in 1:p_indep) {
  cat("\n====", paste0("Extra_Event", j), "====\n")
  true_beta <- beta_list[[j]]
  est_beta <- extra_beta_list[[j]]$beta
  names(true_beta) <- colnames(X_list)
  names(est_beta) <- colnames(X_list)
  print(cbind(True = true_beta, Estimated = round(est_beta, 3)))
}

beta_mat_all<-c(true_betas,beta_list)
#========== bias for beta_hat ==========
l1_bias <- sapply(1:(15+p_indep), function(j) {
  est <- beta_hat_list_all[[j]]$beta
  est[abs(est) < 0.2] <- 0  # 將小於 0.2 的估計值歸零
  true <- beta_mat_all[[j]]
  sum(abs(est - true))
})
l2_bias <- sapply(1:15+p_indep, function(j) {
  est <- beta_hat_list_all[[j]]$beta
  est[abs(est) < 0.2] <- 0  
  true <- beta_mat_all[[j]]
  sqrt(sum((est - true)^2))
})
l2_bias <- sum(l2_bias)
l1_bias <- sum(l1_bias)
cat("L1-norm_bias", round(l1_bias, 3), "\n")
cat("L2-norm_bias", round(l2_bias, 3), "\n")
#==========Spe & Sen for beta_hat ==========
true_betas_block2 <- lapply(1:5, function(j) beta_list[[j]])
true_betas_all <- c(true_betas, true_betas_block2)  

#整體的
threshold <- 0.2
beta_hat_all_flat <- unlist(lapply(beta_hat_list_all, function(x) abs(x$beta) >= threshold))
beta_true_all_flat <- unlist(lapply(true_betas_all, function(x) abs(x) > 1e-6))
TP <- sum(beta_hat_all_flat & beta_true_all_flat)
FP <- sum(beta_hat_all_flat & !beta_true_all_flat)
FN <- sum(!beta_hat_all_flat & beta_true_all_flat)
TN <- sum(!beta_hat_all_flat & !beta_true_all_flat)
overall_sen <- TP / (TP + FN)
overall_spe <- TN / (TN + FP)
cat("Overall Sen:", round(overall_sen, 3), "\n")
cat("Overall Spe:", round(overall_spe, 3), "\n")
#

#每個events各自計算後平均
evaluate_beta_threshold <- function(beta_hat, beta_true, threshold = 0.2) {
  selected <- abs(beta_hat) >= threshold
  true <- abs(beta_true) > 1e-6
  TP <- sum(selected & true)
  FP <- sum(selected & !true)
  FN <- sum(!selected & true)
  TN <- sum(!selected & !true)
  sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA
  return(c(Sensitivity = sensitivity, Specificity = specificity))
}
eval_matrix <- sapply(1:20, function(j) {
  beta_hat <- beta_hat_list_all[[j]]$beta
  beta_true <- true_betas_all[[j]]
  evaluate_beta_threshold(beta_hat, beta_true, threshold = 0.2)
})
eval_df <- as.data.frame(t(eval_matrix))
rownames(eval_df) <- c(paste0("Event", 1:15), paste0("Extra_Event", 1:5))
print(eval_df)                     
colMeans(eval_df, na.rm = TRUE)    
#



#==========F_hat==========
K_h <- function(z, h) dnorm(z / h) / h
compute_F_vector <- function(y, delta, Xproj, beta_hat, h = NULL) {
  n <- length(y)
  if (is.null(h)) { 
    h <- n^(-1/5)
    if (h <= 0) h <- 1e-3
  }
  # Kaplan–Meier for censoring
  fit_G <- survfit(Surv(y, 1 - delta) ~ 1)
  G_fun <- stepfun(fit_G$time, c(1, fit_G$surv))
  G_y <- G_fun(y)
  G_y[G_y < 1e-8] <- 1e-8
  
  # weight
  weight <- delta / G_y
  F_hat <- numeric(n)
  
  for (i in 1:n) {
    num <- sum(weight * (y <= y[i]) * K_h(Xproj - Xproj[i], h))
    den <- sum(K_h(Xproj - Xproj[i], h))
    F_hat[i] <- num / den
  }
  return(F_hat)
  
}

estimate_graph_from_F <- function(beta_hat_list, x_list, y_mat, delta_mat, rho = 0.1) {
  m <- length(beta_hat_list)
  n <- nrow(x_list[[1]])
  F_list <- vector("list", m)
  for (j in 1:m) {
    beta <- beta_hat_list_all[[j]]$beta
    Xj <- x_list[[j]]  
    Xproj <- as.vector(Xj %*% beta)
    yj <- y_mat[, j]
    delta <- delta_mat[, j]
    F_list[[j]] <- compute_F_vector(yj, delta, Xproj, beta_hat_list[[j]])
  }
  
  
  
  D <- matrix(0, m, m)
  for (j in 1:m) {
    for (k in j:m) {
      dval <- dcor(F_list[[j]], F_list[[k]])
      D[j, k] <- dval
      D[k, j] <- dval
    }
  }
  threshold<-seq(1e-5,0.1,by=0.0005)
  gla <- glasso(D, rho = rho)
  #
  adj_list <- list()
  metric_list <- list()
  for (thresh in threshold) {
    adj <- (abs(gla$wi) > thresh) * 1
    diag(adj) <- 0
    
    P1<-sim_result$B*1
    P2 <- diag(0, p_indep)
    P <- as.matrix(bdiag(P1, P2))
    adj_matrix_ture<-(abs(P) > 0.1) * 1
    
    sen <- sum(adj & adj_matrix_ture) / sum(adj_matrix_ture)
    spe <- sum(!adj & !adj_matrix_ture) / sum(!adj_matrix_ture)
    bal_acc <- (sen + spe) / 2
    edge_count <- sum(adj) / 2
    
    adj_list[[as.character(thresh)]] <- adj
    metric_list[[as.character(thresh)]] <- c(
      threshold = thresh,
      sensitivity = sen,
      specificity = spe,
      balanced_accuracy = bal_acc,
      edges = edge_count
    )
  }

  return(list(D = D, gla = gla, adj_list = adj_list, metrics = metric_list))
  
}
event_names <- paste0("T", 1:(15+p_indep))
covariate_map <- c(rep("X", 15), rep("X_list", p_indep))
names(covariate_map) <- event_names

x_list_all <- lapply(event_names, function(name) {
  if (covariate_map[name] == "X") {
    return(matrix(X, nrow = nrow(X), ncol = ncol(X)))
  } else {
    j <- as.integer(sub("T", "", name)) - 15
    if (j >= 1 && j <= length(X_list)) {
      return(X_list[[j]])
    } else {
      stop(sprintf("Index j = %d out of bounds for X_list", j))
    }
  }
})
names(x_list_all) <- event_names

logy<-log(y_all)
res <- estimate_graph_from_F(beta_hat_list_all, x_list_all, logy, delta_all, rho = 0.1)
adj_matrix_est <- res$adj_list

metrics_df <- as.data.frame(do.call(rbind, res$metrics))
best_idx <- which.max(metrics_df$balanced_accuracy)
best_thresh <- metrics_df$threshold[best_idx]
adj_matrix_est <- res$adj_list[[as.character(best_thresh)]]

cat(sprintf("Best threshold (by balanced accuracy): %.3f\n", best_thresh))
cat(sprintf("Sensitivity: %.3f\n", metrics_df$sensitivity[best_idx]))
cat(sprintf("Specificity: %.3f\n", metrics_df$specificity[best_idx]))


P1<-sim_result$B*1
P2 <- diag(0, p_indep)
P <- as.matrix(bdiag(P1, P2))
adj_matrix_ture<-(abs(P) > 0.1) * 1

sen_boost <- sum(adj_matrix_est & adj_matrix_ture) / sum(adj_matrix_ture)
spe_boost<- sum(!adj_matrix_est & !adj_matrix_ture) / sum(!adj_matrix_ture)
cat(sprintf("Sen (Boosting): %.3f\n", sen_boost))
cat(sprintf("Spe (Boosting): %.3f\n", spe_boost))


#=========== Plot networld structure ==========
par(mfrow = c(1, 2))
vertex_colors <- ifelse(as.integer(sub("T", "", event_names)) >= 16, "salmon", "skyblue")
# True network
plot(
  graph_from_adjacency_matrix(adj_matrix_ture, mode = "undirected"),
  main = "True network",
  vertex.label = event_names,
  vertex.color = vertex_colors
)

# Estimated network
plot(
  graph_from_adjacency_matrix(adj_matrix_est, mode = "undirected"),
  main = "Estimated network",
  vertex.label = event_names,
  vertex.color = vertex_colors
)









