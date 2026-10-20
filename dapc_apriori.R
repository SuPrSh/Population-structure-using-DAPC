library(adegenet)
set.seed(123)

#----load dataset
data <- import2genind("infile_mugger.gtx", missing=0, quiet=FALSE)

#----define group
grp <- pop(data)                 # a priori sampling locations
K   <- nlevels(grp)
table(grp)

#---- PCs by a-score, all K-1 DAs
tmp   <- dapc(data, grp, n.pca = min(50, nInd(data) - 1), n.da = K - 1)
dapc1 <- dapc(data, grp, n.pca = optim.a.score(tmp)$best, n.da = K - 1)
round(cumsum(dapc1$eig) / sum(dapc1$eig) * 100, 1)   # % discrimination on axis

col <- funky(K)                  # or c("orange", "green", ...) length K

#----LD1 vs LD2, labelled by location
scatter(dapc1, xax = 1, yax = 2, col = col, label = levels(grp), clabel = 0.8,
        cstar = 0, cellipse = 1.5, pch = 20, cex = 2, solid = 0.6, bg = "white",
        scree.da = TRUE, posi.da = "bottomleft", ratio.da = 0.18,
        scree.pca = TRUE, posi.pca = "bottomright", ratio.pca = 0.18)

#----LD1 vs LD3
scatter(dapc1, xax = 1, yax = 3, col = col, label = levels(grp), clabel = 0.8,
        cstar = 0, pch = 20, cex = 2, solid = 0.6, bg = "white", scree.da = FALSE)

#----Membership probabilities (individuals ordered by location)
compoplot(dapc1, col.pal = function(k) col, show.lab = TRUE,
          lab = indNames(data), txt.leg = levels(grp), n.col = 2)

#----Reassignment success per location
summary(dapc1)$assign.per.pop

#----Save main figure
jpeg("dapc_apriori.jpeg", width = 8, height = 8, units = "in", res = 500)
scatter(dapc1, col = col, label = levels(grp), clabel = 0.8, cstar = 0,
        cellipse = 1.5, pch = 20, cex = 2, solid = 0.6, bg = "white",
        scree.da = FALSE, scree.pca = TRUE, posi.pca = "bottomright", ratio.pca = 0.18)
dev.off()
#----END----#
