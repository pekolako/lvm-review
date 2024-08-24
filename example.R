# Load the packages:
library(gllvm); library(ltm); library(gmf)
# and the Australian ant data (Gibb and Cunningham, 2011) which is included in
# the gmf package.
data("antTraits")
y <- as.matrix(antTraits$abund)
# transform from counts to presence-absence, i.e., 0/1 to facilitate use of ltm
y <- (y>0)*1
set.seed(29052024)

# Fit logistic GLLVMs using EVA with 1 to 10 latent variables, and print the
# resulting information criteria
evas <- list()
criteria <- matrix(NA, ncol=5, nrow=10)
for (i in 1:10) {
  # default link for "binomial" is "probit" so "logit" needs to be specified.
  evai <- gllvm(as.matrix(y), family="binomial", method="EVA", 
                link="logit", num.lv=i)  
  evas[[i]] <- evai
  criteria[i,1] <- AIC(evai)
  criteria[i,2] <- AICc(evai)
  criteria[i,3] <- BIC(evai)
  criteria[i,4] <- evai$logL
  criteria[i,5] <- summary(evai)$df
}
criteria <- cbind(1:10, criteria)
criteria <- as.data.frame(criteria)
colnames(criteria) <- c("p", "AIC", "AICc", "BIC", "logL", "df")
criteria

# For ordination, choose the model with num.lv = 2.
# Next, fit same type of model using the ltm package, but using the same
# starting values for the parameters as in gllvm
eva <- evas[[2]]
u0 <- eva$start$index
b0 <- eva$start$params[,1]  # staring values for fixed column effects
th0 <- t(t(eva$start$params[,2:3]) * eva$start$sigma.lv) # loading matrix
par0 <- cbind(b0, th0)
# Fit using ltm, with gllvm starting values and upper triangle of the loading
# matrix constrained to 0.
ltm <- ltm(y ~ z1 + z2, constraint = matrix(c(1,3,0),byrow=TRUE,ncol=3), 
           start.val=as.matrix(par0), control=list(verbose=TRUE))

# Ordination plots:
setEPS()
postscript("example_ordplots.eps")
par(mfrow=c(2,1),oma=c(1,1,1,1), mar=c(4.1,4.1,2.1,1.1))
fsc0 <- as.matrix(factor.scores(ltm)$score.dat[,c("z1","z2")])
fsc <- fsc0 %*% svd(fsc0)$v  # LVs rotated to their principal direction, which
# is also the default if using gllvm::ordiplot().
colnames(fsc) <- c("z1","z2")
fsc <- as.data.frame(fsc)
z1 <- jitter(fsc$z1, amount=0.1)
z2 <- jitter(fsc$z2, amount=0.1)
# Plot the LVs from the ltm fit:
plot(z1,z2,xlab="Latent variable 1",
     ylab="Latent variable 2", cex=0.6, ylim=c(min(z2)-0.5,max(z2)+0.5),
     xlim=c(min(z1)-0.1,max(z1)+0.1))
mtext("a)", side=3, line=1, adj=-0.1, cex=1)
text(z1,z2, as.character(1:30),cex=0.6,pos=3,col="red")

lvs <- eva$lvs %*% svd(eva$lvs)$v # LVs rotated to their principal direction
colnames(lvs) <- c("z1","z2")
lvs <- as.data.frame(lvs)
z1 <- jitter(lvs$z1, amount=0.1)
z2 <- jitter(lvs$z2, amount=0.1)
# Plot the LVs from the gllvm fit:
plot(z1,z2,xlab="Latent variable 1",
     ylab="Latent variable 2", cex=0.6, ylim=c(min(z2)-0.5,max(z2)+0.5),
     xlim=c(min(z1)-0.1,max(z1)+0.1))
mtext("b)", side=3, line=1, adj=-0.1, cex=1)
text(z1,z2, as.character(1:30),cex=0.6,pos=3,col="red")
dev.off()

# Residual correlation plot:
setEPS()
postscript("example_corplot.eps")
par(mfrow=c(1,1))
# The following packages are needed:
library(corrplot); library(gclus) 
cr <- getResidualCor(evas[[2]]) # eva with 1 lv 
corrplot(cr[order.single(cr), order.single(cr)], diag = FALSE, type = "lower", 
         method = "square", tl.cex = 0.6, tl.srt = 45, tl.col = "black")

dev.off()

# Latex table from the information criterias computed earlier
library(xtable)
print(xtable(criteria),include.rownames = FALSE)

# Residual plots as provided by gllvm
setEPS()
postscript("example_residuals.eps")
#ordiplot(eva, biplot=T)
par(mfrow=c(2,2),oma=c(1,1,1,1), mar=c(4.1,4.1,2.1,1.1))
plot(eva, var.colors = 1, which = 1:4)
dev.off()







