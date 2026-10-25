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

#---- User-defined K, PCs and DAs
grp <- find.clusters(data, n.pca = nInd(data) - 1, max.n.clust = 40)   # pick K at BIC elbow
table(pop(data), grp$grp) #individual distribution acrosss cluster
K <- nlevels(grp$grp)     #assigning cluster
dapc2 <- dapc(data, grp$grp)          # interactive: choose PCs, then DAs
summary(dapc2)
col <- funky(K)            # auto choose colors or use col <- c("#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A")
scatter(dapc2, posi.da = "bottomleft", scree.pca = TRUE, posi.pca = "topright",
        bg.inset = "white", ratio.da = .18, ratio.pca = .18, col = col, clabel = 0.7, cstar = 1.0, cellipse = 1.5,
        pch = 20, cex = 2, solid = 0.9, bg = "white") #change eigenvalue inset position by changing "topright" to any position you like e.g. "bottomright" "bottomleft"
dev.copy(jpeg, "dapc2A.jpeg", width = 10, height = 10, units = "in", res = 500); dev.off() #run to save image

#----DAPC with number of PCs chosen by a-score
maxPC <- min(40, nInd(data) - 1)
tmp   <- dapc(data, grp$grp, n.pca = maxPC, n.da = K - 1)
asc   <- optim.a.score(tmp, n.pca = 1:maxPC, n.sim = 100)
best <- max(asc$best, 2)                                    # avoid the 1-PC artefact

dapc3 <- dapc(data, grp$grp, n.pca = best, n.da = K - 1)

#--assign color
col <- funky(K)                          # or your own colours, length K
lab <- paste("Cluster", 1:K)             # or c("BA-NG-GD", "...", ...), length K
ids <- indNames(data)

#--Scatter with cluster
scatter(dapc3, posi.da = "bottomleft", scree.pca = TRUE, posi.pca = "topright",
        bg.inset = "white", ratio.da = .18, ratio.pca = .18, col = col, clabel = 0.5, cstar = 1.0, cellipse = 1.5,
        pch = 20, cex = 2, solid = 0.9, bg = "white")
dev.copy(jpeg, "dapc3A.jpeg", width = 8, height = 8, units = "in", res = 500); dev.off() #run to save image

#----Assignment plot
assignplot(dapc3)
dev.copy(jpeg, "assignment.jpeg", width = 10, height = 10, units = "in", res = 500); dev.off() #run to save image
dev.off()
