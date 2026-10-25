library(adegenet)
set.seed(123)

#----Load genepop file: skip this if you have gtx file
data <- read.genepop("infile-GENEPOP.gen", ncode = 3)   # ncode = digits per allele (2 or 3)
data
#--sanity check for pop and individuals
nPop(data)              # number of populations
table(pop(data))        # individuals per population
levels(pop(data)) <- c("Site1", "Site2", "Site3", "Site4")

#---load gtx file: DO NOT RUN LINE 12 and 13 if you load from genepop file
#data <- import2genind("infile_mugger.gtx", missing=0, quiet=FALSE)
#data

#----define group
grp <- pop(data)                 # a priori sampling locations
K   <- nlevels(grp)
table(grp)

#---- PCs by a-score, all K-1 DAs
tmp   <- dapc(data, grp, n.pca = min(50, nInd(data) - 1), n.da = K - 1)
dapc1 <- dapc(data, grp, n.pca = optim.a.score(tmp)$best, n.da = K - 1)
round(cumsum(dapc1$eig) / sum(dapc1$eig) * 100, 1)   # % discrimination on axis

col <- funky(K)                  # or c("orange", "green", ...) length K
lab <- c("Pali", "Baroli", "Nadigaon", "Tigri")
ids <- indNames(data)

#----Scatter with cluster labels
scatter(dapc1, col = col, label = lab, clabel = 1, cstar = 1, cellipse = 1.5,
        pch = 20, cex = 2.5, solid = 0.5, bg = "white",
        scree.pca = TRUE, posi.pca = "topright", posi.da = "topleft",
        ratio.pca = 0.18, ratio.da = 0.18)
dev.copy(jpeg, "dapc_ap_1A.jpeg", width = 10, height = 10, units = "in", res = 500); dev.off() #run to save image
#----Plot othe LDs
#----LD1 vs LD2, labelled by location
scatter(dapc1, xax = 1, yax = 2, col = col, label = levels(grp), clabel = 0.8,
        cstar = 0, cellipse = 1.5, pch = 20, cex = 2, solid = 0.6, bg = "white",
        scree.da = TRUE, posi.da = "bottomleft", ratio.da = 0.18,
        scree.pca = TRUE, posi.pca = "bottomright", ratio.pca = 0.18)
dev.copy(jpeg, "dapc_ap_1B.jpeg", width = 10, height = 10, units = "in", res = 500); dev.off() #run to save image
#----LD1 vs LD3
scatter(dapc1, xax = 1, yax = 3, col = col, label = levels(grp), clabel = 0.8,
        cstar = 0, pch = 20, cex = 2, solid = 0.6, bg = "white", scree.da = FALSE)
dev.copy(jpeg, "dapc_ap_1C.jpeg", width = 10, height = 10, units = "in", res = 500); dev.off() #run to save image
#----Membership probabilities (individuals ordered by location)
compoplot(dapc1, col.pal = function(k) col, show.lab = TRUE,
          lab = indNames(data), txt.leg = levels(grp), n.col = 2)
dev.copy(jpeg, "complot_ap_1.jpeg", width = 10, height = 4, units = "in", res = 500); dev.off() #run to save image
#----Reassignment success per location
summary(dapc1)$assign.per.pop

dev.off()
#----END----#
