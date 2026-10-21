library(adegenet)
set.seed(123)
#----Load genepop file: skip this if you have gtx file
data <- read.genepop("infile-GENEPOP.gen", ncode = 3)   # ncode = digits per allele (2 or 3)
data
#--sanity check for pop and individuals
nPop(data1)              # number of populations
table(pop(data1))        # individuals per population
levels(pop(data1)) <- c("Site1", "Site2", "Site3", "Site4")

#---load gtx file: DO NOT RUN LINE 12 and 13 if you load from genepop file
#data <- import2genind("infile_mugger.gtx", missing=0, quiet=FALSE)
#data

#---- User-defined K, PCs and DAs
grp <- find.clusters(data, n.pca = nInd(data) - 1, max.n.clust = 40)   # pick K at BIC elbow
table(pop(data), grp$grp)
K <- nlevels(grp$grp)

dapc2 <- dapc(data, grp$grp)          # interactive: choose PCs, then DAs

summary(dapc2)
scatter(dapc2, posi.da = "topleft", scree.pca = TRUE, posi.pca = "topright",
        bg.inset = "white", ratio.da = .18, ratio.pca = .18)

#----DAPC with number of PCs chosen by a-score
tmp  <- dapc(data, grp$grp, n.pca = 40, n.da = K - 1)       # enough PCs to test 1:40
asc  <- optim.a.score(tmp, n.pca = 1:40, n.sim = 100)
best <- max(asc$best, 2)                                    # avoid the 1-PC artefact

dapc3 <- dapc(data, grp$grp, n.pca = best, n.da = K - 1)

#--assign color
col <- funky(K)                          # or your own colours, length K
lab <- paste("Cluster", 1:K)             # or c("BA-NG-GD", "...", ...), length K
ids <- indNames(data)

#--Scatter with cluster
scatter(dapc3, col=col, posi.da = "topleft", scree.pca = TRUE, posi.pca = "topright",
        bg.inset = "white", ratio.da = .18, ratio.pca = .18)


#--Scatter with individual IDs (no insets, so text() lands correctly)
scatter(dapc3, col = col, label = lab, clabel = 1, cstar = 0,
        pch = 20, cex = 1.5, solid = 0.5, bg = "white", scree.da = FALSE)
text(dapc3$ind.coord[, 1], dapc3$ind.coord[, 2], ids, cex = 0.5, pos = 3)

#----Membership probabilities with individual IDs
compoplot(dapc3, col.pal = function(k) col, show.lab = TRUE, lab = ids,
          txt.leg = lab, n.col = K, posi = "bottomright")

#----Assignment plot
assignplot(dapc3)

#----Save main figure
jpeg("dapc_denovo.jpeg", width = 8, height = 8, units = "in", res = 500)
scatter(dapc2, col = col, label = levels(grp), clabel = 0.8, cstar = 0,
        cellipse = 1.5, pch = 20, cex = 2, solid = 0.6, bg = "white",
        scree.da = FALSE, scree.pca = TRUE, posi.pca = "bottomright", ratio.pca = 0.18)
dev.off()
