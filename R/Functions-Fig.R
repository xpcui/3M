require(dplyr)
require(ggplot2)

## POWER
power_comp <- function(RE, title, legend.position){
  df <- RE
  df$method <- factor(df$method, levels = c("3M","P-value Ranking","Lasso","IT"))
  
  m0 <- c("3M","P-value Ranking","Lasso","IT")
  m0 <- factor(m0, levels=m0)
  myColors <- c("red","grey","blue","orange")
  names(myColors) <- levels(m0)
  colScale <- scale_colour_manual(name = "Method", values = myColors)
  colScale1 <- scale_fill_manual(name = "Method", values = myColors)
  
  p <- ggplot(data=df, aes(x=beta, y=mean, col=method, group=method)) +
    geom_point(aes(colour=method), size=5, position=position_dodge(0.1))+
    geom_smooth(data=df,aes(x=beta, y=mean, group=method), method = "lm", formula = y ~ poly(x, 3), se = FALSE, position=position_dodge(0.1)) +
    labs(x = "Predictive effect", y = "Power", title = title) + 
    coord_cartesian(xlim = c(1, 3), ylim = c(0, 1)) +
    scale_x_continuous(breaks = seq(1, 3, 0.5)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    theme_bw() +
    theme(legend.position = legend.position,
          plot.title = element_text(size=22),
          strip.text = element_text(size = 22),
          axis.text = element_text(size = 22),
          axis.title = element_text(size = 22),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 22),
          legend.key.size = unit(1.5,"cm"))
  
  p <- p + facet_grid(cols=vars(n), scales = "fixed", labeller = label_both) + colScale + colScale1
  p
}

power_comp_legend <- function(RE, title, legend.position){
  df <- RE
  df$method <- factor(df$method, levels = c("3M","P-value Ranking","Lasso","IT"), labels = c(" 3M  "," P-value Ranking  "," Lasso  "," IT  "))
  
  m0 <- c(" 3M  "," P-value Ranking  "," Lasso  "," IT  ")
  m0 <- factor(m0, levels=m0)
  myColors <- c("red","grey","blue","orange")
  names(myColors) <- levels(m0)
  colScale <- scale_colour_manual(name = "Method  ", values = myColors)
  colScale1 <- scale_fill_manual(name = "Method  ", values = myColors)
  
  p <- ggplot(data=df, aes(x=beta, y=mean, col=method, group=method)) +
    geom_point(aes(colour=method), size=5, position=position_dodge(0.1))+
    geom_smooth(data=df,aes(x=beta, y=mean, group=method), method = "lm", formula = y ~ poly(x, 3), se = FALSE, position=position_dodge(0.1)) +
    labs(x = "Predictive effect", y = "Power", title = title) + 
    coord_cartesian(xlim = c(1, 3), ylim = c(0, 1)) +
    scale_x_continuous(breaks = seq(1, 3, 0.5)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    theme_bw() +
    theme(legend.position = legend.position,
          plot.title = element_text(size=22),
          strip.text = element_text(size = 22),
          axis.text = element_text(size = 22),
          axis.title = element_text(size = 22),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 22),
          legend.key.size = unit(1.5,"cm"))
  
  p <- p + facet_grid(cols=vars(n), scales = "fixed", labeller = label_both) + colScale + colScale1
  p
}

error_comp <- function(mean, e0){
  type0 <- c("Error A","Error B","Error A or B")
  type <- rep(type0, each = length(e0))
  effect <- rep(e0, 3)
  
  df <- cbind.data.frame(type=type, effect=effect, mean=mean)
  df$type <- factor(df$type, levels = c("Error A","Error B","Error A or B"),
                    labels = c(expression(paste("P(",E[A]^c,")")),
                               expression(paste("P(",E[B]^c,")")),
                               expression(paste("P(",E[A]^c,"U",E[B]^c,")"))))
  
  p <- ggplot(data=df, aes(x=effect, y=mean)) +
    geom_point(aes(colour="red"), size=5) +
    geom_smooth(data=df,aes(x=effect, y=mean, colour = "red"), method = "lm", formula = y ~ poly(x, 3), se = FALSE) +
    geom_hline(yintercept=0.05, linetype="dashed", color = "black", linewidth=1.2) +
    coord_cartesian(xlim = c(1, 1.5),ylim = c(0, 0.05)) +
    scale_x_continuous(breaks = seq(1, 1.5, 0.1)) +
    scale_y_continuous(breaks = seq(0, 0.05, 0.01)) +
    labs(x = "Predictive effect", y = "Error rate", title = "") +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(color="black", size=25),
          strip.text = element_text(size = 25),
          axis.text = element_text(size = 25),
          axis.title = element_text(size = 25)) +
    facet_grid(cols=vars(type), scales = "fixed", label = "label_parsed")
  p
}

plot.est <- function(est, N){
  mean0 <- apply(est, 2, mean)
  se0 <- apply(est, 2, sd)/sqrt(N)
  
  df <- cbind.data.frame(N=N, mean0=mean0, se0=se0)
  
  p <- ggplot(df, aes(x=N, y=mean0)) + 
    geom_line() +
    geom_point(color="red", size = 5) +
    geom_errorbar(aes(ymin=mean0-se0, ymax=mean0+se0), width=.2, position=position_dodge(0.05)) +
    labs(x = "n", y = expression(widehat(beta)^"*")) +
    coord_cartesian(ylim = c(0, 1.2)) +
    scale_y_continuous(breaks=seq(0, 1.2, 0.2)) +
    scale_x_continuous(breaks=seq(0, 1000, 100)) +
    geom_hline(yintercept = 1, color="orange", linetype = "longdash", linewidth = 0.4) +
    theme_bw() +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 25))
  
  p
}

plot.est.ols <- function(est, N){
  mean0 <- apply(est, 2, mean, na.rm=T)
  sd0 <- apply(est, 2, sd, na.rm=T)
  
  df <- cbind.data.frame(N=N, mean0=mean0, sd0=sd0)
  
  p <- ggplot(df, aes(x=log10(N), y=mean0)) + 
    geom_line() +
    geom_point(color="red", size = 5) +
    geom_errorbar(aes(ymin=mean0-sd0, ymax=mean0+sd0), width=.0) +
    labs(x = "log10(N)", y = expression(C[1]==C[1](Z))) +
    coord_cartesian(ylim = c(-0.3, 0.3)) +
    scale_y_continuous(breaks=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3)) +
    geom_hline(yintercept = 0, color="orange", linetype = "longdash", size = 0.4) +
    theme_bw() +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 25))
  
  p
}

plot.bar <- function(mean){
  m0 <- c("3M","P-value Ranking","Lasso","IT")
  size0 <- c("n == 379","n == 500","n == 1000")
  e0 <- c("beta[1] == 3","beta[1] == 4","beta[1] == 5")
  
  m0 <- factor(m0, levels=m0)
  myColors <- c("red","grey","blue","orange")
  names(myColors) <- levels(m0)
  colScale <- scale_colour_manual(name = "Method", values = myColors)
  colScale1 <- scale_fill_manual(name = "Method", values = myColors)
  
  method <- c(rep(m0, each = 3),rep(m0, each = 3),rep(m0, each = 3))
  size <- rep(size0, each = 4*3)
  effect <- rep(e0, 3*4)
  
  df <- cbind.data.frame(effect=effect, mean=mean, method=method, size=size)
  
  df$effect <- factor(df$effect, levels = e0)
  df$size <- factor(df$size, levels = size0)
  
  p <- ggplot(data=df, aes(x=method, y=mean, fill=method)) +
    geom_bar(stat="identity", position=position_dodge()) +
    coord_cartesian(ylim=c(0.4, 1)) +
    scale_y_continuous(breaks = seq(0.4,1,0.2)) +
    theme(legend.position = "none",
          strip.text = element_text(size = 25),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust = 1),
          axis.title = element_text(size = 20)) +
    labs(x="",y="Power",title = "") +
    facet_grid(size ~ effect, scales = "fixed", labeller = label_parsed) + colScale + colScale1
  p
}

plot.QQ <- function(select.seq){
  sample1 <- select.seq[,1]; s1 <- quantile(sample1, probs = seq(0,1,0.01))
  sample2 <- select.seq[,2]; s2 <- quantile(sample2, probs = seq(0,1,0.01))
  sample3 <- select.seq[,3]; s3 <- quantile(sample3, probs = seq(0,1,0.01))
  sample4 <- select.seq[,4]; s4 <- quantile(sample4, probs = seq(0,1,0.01))
  theory0 <- select.seq[,5]; t0 <- quantile(theory0, probs = seq(0,1,0.01))
  
  df <- cbind.data.frame(x=s1, y=t0)
  p1 = ggplot(df, aes(x=x, y=y)) + geom_point() +
    geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=0.8) +
    labs(x="Theoretical Quantiles",y="Sample Quantiles",title = "A") +
    theme(axis.title = element_text(size=20),
          axis.text = element_text(size=20),
          plot.title = element_text(size=20))
  
  df <- cbind.data.frame(x=s2, y=t0)
  p2 = ggplot(df, aes(x=x, y=y)) + geom_point() +
    geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=0.8) +
    labs(x="Theoretical Quantiles",y="Sample Quantiles",title = "B") +
    theme(axis.title = element_text(size=20),
          axis.text = element_text(size=20),
          plot.title = element_text(size=20))
  
  df <- cbind.data.frame(x=s3, y=t0)
  p3 = ggplot(df, aes(x=x, y=y)) + geom_point() +
    geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=0.8) +
    labs(x="Theoretical Quantiles",y="Sample Quantiles",title = "C") +
    theme(axis.title = element_text(size=20),
          axis.text = element_text(size=20),
          plot.title = element_text(size=20))
  
  df <- cbind.data.frame(x=s4, y=t0)
  p4 = ggplot(df, aes(x=x, y=y)) + geom_point() +
    geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=0.8) +
    labs(x="Theoretical Quantiles",y="Sample Quantiles",title = "D") +
    theme(axis.title = element_text(size=20),
          axis.text = element_text(size=20),
          plot.title = element_text(size=20))
  
  grid.arrange(arrangeGrob(p1,p2,p3,p4,ncol=2), heights=c(10, 1))
}

















