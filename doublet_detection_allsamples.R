wt401<- FindNeighbors(wt401, dims = 1:40)
wt401<- FindClusters(wt401, resolution = 1.0)
wt401<-RunUMAP(wt401, dims = 1:40)
DimPlot(wt401)
sweep.wt401<- paramSweep_v3(wt401, PCs = 1:40, sct = FALSE)
annotations<-wt401$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nexp_wt401<-round(0.08*ncol(wt401))
nexp_adj_wt401<-round(nexp_wt401*(1-homotypic.prop))
sweep.stats_wt401<- summarizeSweep(sweep.wt401, GT = FALSE)
bcmvn_wt401 <- find.pK(sweep.stats_wt401)
bcmvn.max_wt401 <- bcmvn_wt401[which.max(bcmvn_wt401$BCmetric),]

wt433<- FindNeighbors(wt433, dims = 1:40)
wt433<- FindClusters(wt433, resolution = 1.0)
wt433<-RunUMAP(wt433, dims = 1:40)
DimPlot(wt433)
sweep.wt433<- paramSweep_v3(wt433, PCs = 1:40, sct = FALSE)
annotations<-wt433$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nexp_wt433<-round(0.08*ncol(wt433))
nexp_adj_wt433<-round(nexp_wt433*(1-homotypic.prop))
sweep.stats_wt433<- summarizeSweep(sweep.wt433, GT = FALSE)
bcmvn_wt433 <- find.pK(sweep.stats_wt433)
bcmvn.max_wt433 <- bcmvn_wt433[which.max(bcmvn_wt433$BCmetric),]

wt438<- FindNeighbors(wt438, dims = 1:40)
wt438<- FindClusters(wt438, resolution = 1.0)
wt438<-RunUMAP(wt438, dims = 1:40)
DimPlot(wt438)
sweep.wt438<- paramSweep_v3(wt438, PCs = 1:40, sct = FALSE)
annotations<-wt438$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nexp_wt438<-round(0.08*ncol(wt438))
nexp_adj_wt438<-round(nexp_wt438*(1-homotypic.prop))
sweep.stats_wt438<- summarizeSweep(sweep.wt438, GT = FALSE)
bcmvn_wt438 <- find.pK(sweep.stats_wt438)
bcmvn.max_wt438 <- bcmvn_wt438[which.max(bcmvn_wt438$BCmetric),]

wt555<- FindNeighbors(wt555, dims = 1:40)
wt555<- FindClusters(wt555, resolution = 1.0)
wt555<-RunUMAP(wt555, dims = 1:40)
DimPlot(wt555)
sweep.wt555<- paramSweep_v3(wt555, PCs = 1:40, sct = FALSE)
annotations<-wt555$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nexp_wt555<-round(0.08*ncol(wt555))
nexp_adj_wt555<-round(nexp_wt555*(1-homotypic.prop))
sweep.stats_wt555<- summarizeSweep(sweep.wt555, GT = FALSE)
bcmvn_wt555 <- find.pK(sweep.stats_wt555)
bcmvn.max_wt555 <- bcmvn_wt555[which.max(bcmvn_wt555$BCmetric),]

ko3<- FindNeighbors(ko3, dims = 1:40)
ko3<- FindClusters(ko3, resolution = 1.0)
ko3<-RunUMAP(ko3, dims = 1:40)
DimPlot(ko3)
sweep.ko3<- paramSweep_v3(ko3, PCs = 1:40, sct = FALSE)
annotations<-ko3$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nexp_ko3<-round(0.08*ncol(ko3))
nexp_adj_ko3<-round(nexp_ko3*(1-homotypic.prop))
sweep.stats_ko3<- summarizeSweep(sweep.ko3, GT = FALSE)
bcmvn_ko3 <- find.pK(sweep.stats_ko3)
bcmvn.max_ko3 <- bcmvn_ko3[which.max(bcmvn_ko3$BCmetric),]

ko12<- FindNeighbors(ko12, dims = 1:40)
ko12<- FindClusters(ko12, resolution = 1.0)
ko12<-RunUMAP(ko12, dims = 1:40)
DimPlot(ko12)
sweep.ko12<- paramSweep_v3(ko12, PCs = 1:40, sct = FALSE)
annotations<-ko12$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nexp_ko12<-round(0.08*ncol(ko12))
nexp_adj_ko12<-round(nexp_ko12*(1-homotypic.prop))
sweep.stats_ko12<- summarizeSweep(sweep.ko12, GT = FALSE)
bcmvn_ko12 <- find.pK(sweep.stats_ko12)
bcmvn.max_ko12 <- bcmvn_ko12[which.max(bcmvn_ko12$BCmetric),]

ko397<- FindNeighbors(ko397, dims = 1:40)
ko397<- FindClusters(ko397, resolution = 1.0)
ko397<-RunUMAP(ko397, dims = 1:40)
DimPlot(ko397)
sweep.ko397<- paramSweep_v3(ko397, PCs = 1:40, sct = FALSE)
annotations<-ko397$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nexp_ko397<-round(0.08*ncol(ko397))
nexp_adj_ko397<-round(nexp_ko397*(1-homotypic.prop))
sweep.stats_ko397<- summarizeSweep(sweep.ko397, GT = FALSE)
bcmvn_ko397 <- find.pK(sweep.stats_ko397)
bcmvn.max_ko397 <- bcmvn_ko397[which.max(bcmvn_ko397$BCmetric),]

ko494<- FindNeighbors(ko494, dims = 1:40)
ko494<- FindClusters(ko494, resolution = 1.0)
ko494<-RunUMAP(ko494, dims = 1:40)
DimPlot(ko494)
sweep.ko494<- paramSweep_v3(ko494, PCs = 1:40, sct = FALSE)
annotations<-ko494$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nexp_ko494<-round(0.08*ncol(ko494))
nexp_adj_ko494<-round(nexp_ko494*(1-homotypic.prop))
sweep.stats_ko494<- summarizeSweep(sweep.ko494, GT = FALSE)
bcmvn_ko494 <- find.pK(sweep.stats_ko494)
bcmvn.max_ko494 <- bcmvn_ko494[which.max(bcmvn_ko494$BCmetric),]

wt401<- doubletFinder_v3(wt401, PCs = 1:40, pN = 0.25, pK = 0.01 , nExp = nexp_adj_wt401, reuse.pANN = FALSE, sct = FALSE)
wt433<- doubletFinder_v3(wt433, PCs = 1:40, pN = 0.25, pK = 0.005 , nExp = nexp_adj_wt433, reuse.pANN = FALSE, sct = FALSE)
wt438<- doubletFinder_v3(wt438, PCs = 1:40, pN = 0.25, pK = 0.01 , nExp = nexp_adj_wt438, reuse.pANN = FALSE, sct = FALSE)
wt555<- doubletFinder_v3(wt555, PCs = 1:40, pN = 0.25, pK = 0.01 , nExp = nexp_adj_wt555, reuse.pANN = FALSE, sct = FALSE)
ko3<- doubletFinder_v3(ko3, PCs = 1:40, pN = 0.25, pK = 0.02 , nExp = nexp_adj_ko3, reuse.pANN = FALSE, sct = FALSE)
ko12<- doubletFinder_v3(ko12, PCs = 1:40, pN = 0.25, pK = 0.01 , nExp = nexp_adj_ko12, reuse.pANN = FALSE, sct = FALSE)
ko397<- doubletFinder_v3(ko397, PCs = 1:40, pN = 0.25, pK = 0.01 , nExp = nexp_adj_ko397, reuse.pANN = FALSE, sct = FALSE)
ko494<- doubletFinder_v3(ko494, PCs = 1:40, pN = 0.25, pK = 0.005 , nExp = nexp_adj_ko494, reuse.pANN = FALSE, sct = FALSE)

##code for doubletFinder run on merged Seurat object
annotations<-p8str.merged$seurat_clusters
annotations<-p8str.merged$seurat_clusters
nexp_p8str.merged<-round(0.08*ncol(p8str.merged))
nexp_adj_p8str.merged<-round(nexp_p8str.merged*(1-homotypic.prop))
sweep.stats_p8str.merged<- summarizeSweep(sweep.p8str.merged, GT = FALSE)
bcmvn.max_p8str.merged <- bcmvn_p8str.merged[which.max(bcmvn_p8str.merged$BCmetric),]
p8str.merged<- doubletFinder_v3(p8str.merged, PCs = 1:40, pN = 0.25, pK = 0.005 , nExp = nexp_adj_p8str.merged, reuse.pANN = FALSE, sct = FALSE)
                         
