results_folder = '~/Desktop/Referenced single cell-RNA sequencing data/E14.5 mouse fetal liver sc-RNA sequencing/3/2'

library(Giotto)

instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  python_path = '/Users/xuhan/Library/r-miniconda-arm64/envs/giotto_env/bin/pythonw')


giotto_SC<-createGiottoObject(expression = readExprMatrix("~/Desktop/Referenced single cell-RNA sequencing data/E14.5 mouse fetal liver sc-RNA sequencing/3/GSM4647231_FetalLiver_1_dge.txt"),
                              instructions = instrs)


giotto_SC<-filterGiotto(gobject = giotto_SC,
                        expression_threshold = 1,
                        feat_det_in_min_cells = 10,
                        min_det_feats_per_cell = 500,
                        expression_values = c('raw'),
                        verbose = T)


giotto_SC <- normalizeGiotto(gobject = giotto_SC, scalefactor = 10000)


library(rtracklayer)
gtf <- import("~/Desktop/Referenced single cell-RNA sequencing data/E14.5 mouse fetal liver sc-RNA sequencing/3/Mus_musculus.GRCm39.108.gtf.gz")
gtf <- gtf[gtf$gene_name!="" & !is.na(gtf$gene_name)]
mito <- gtf$gene_name[as.character(seqnames(gtf)) %in% "Mt"]
mito<-unique(mito)

giotto_SC<-addFeatsPerc(
  giotto_SC,
  feats = mito,
  vector_name = 'perc_mito'
)

giotto_SC<-subsetGiotto(giotto_SC,
                        cell_ids = pDataDT(giotto_SC)[which(pDataDT(giotto_SC)$perc_mito < 20),]$cell_ID)


giotto_SC <- addStatistics(gobject = giotto_SC, expression_values = 'raw')


giotto_SC <- calculateHVF(gobject = giotto_SC)
giotto_SC <- runPCA(gobject = giotto_SC, center = TRUE, scale_unit = TRUE)
screePlot(giotto_SC, ncp = 50, save_param = list(save_name = '3_scree_plot'))


showGiottoDimRed(giotto_SC)
giotto_SC <- createNearestNetwork(gobject = giotto_SC,
                                  dim_reduction_to_use = 'pca', dim_reduction_name = 'pca',
                                  dimensions_to_use = 1:10, k = 15)


giotto_SC = runUMAP(giotto_SC, dimensions_to_use = 1:10)


giotto_SC <- doLeidenCluster(gobject = giotto_SC, resolution = 0.4, n_iterations = 1000)


plotUMAP(gobject = giotto_SC,
         cell_color = 'leiden_clus', show_NN_network = F, point_size = 2.5,
         save_param = list(save_name = "4_Cluster"))

markers_scran = findMarkers_one_vs_all(gobject=giotto_SC, method="scran",
                                       expression_values="normalized", cluster_column='leiden_clus', min_feats=3)
markergenes_scran = unique(markers_scran[, head(.SD, 3), by="cluster"][["feats"]])

plotMetaDataHeatmap(giotto_SC, expression_values = "normalized", metadata_cols = 'leiden_clus',
                    selected_feats = markergenes_scran,
                    y_text_size = 8, show_values = 'zscores_rescaled',
                    save_param = list(save_name = '5_a_metaheatmap'))


topgenes_scran = markers_scran[, head(.SD, 3), by = 'cluster']$feats

violinPlot(giotto_SC, feats = unique(topgenes_scran), cluster_column = 'leiden_clus',
           strip_text = 10, strip_position = 'right',
           save_param = list(save_name = '5_b_violinplot_scran', base_width = 5))



dimFeatPlot2D(giotto_SC, feats = c("Cd34","Lyz2","Ngp","C1qb"), cow_n_col = 2, save_param = list(save_name = "6_featureplot"))


fetalliver_labels<-c("Erythroid progenitor cell",#1
                     "Erythroid progenitor cell",#2
                     "Stem cell",#3
                     "Erythroid progenitor cell",#4
                     "Erythroid progenitor cell",#5
                     "Erythroid progenitor cell",#6
                     "Neutrophil",#7
                     "Other Macrophage",#8
                     "Erythroid progenitor cell",#9
                     "Erythroid cell",#10
                     "Erythroid progenitor cell",#11
                     "Erythroid progenitor cell",#12
                     "C1q+ Macrophage",#13
                     "Hepatocyte",#14
                     "Mast cell",#15
                     "Megakaryocyte",#16
                     "Erythroid progenitor cell")#17
names(fetalliver_labels)<-1:17
giotto_SC<-annotateGiotto(gobject = giotto_SC, annotation_vector = fetalliver_labels ,
                          cluster_column = 'leiden_clus', name = 'fetalliver_labels')
dimPlot2D(gobject = giotto_SC,     dim_reduction_name = 'umap',
          cell_color = "fetalliver_labels", show_NN_network = F, point_size = 2.5,
          save_param = list(save_name = "7_Annotation"))

plotMetaDataHeatmap(giotto_SC, selected_feats = unique(topgenes_scran),
                    metadata_cols = c('fetalliver_labels'),
                    x_text_size = 18,
                    x_text_angle = 45,
                    y_text_size = 18,
                    strip_text_size = 18)

plotMetaDataHeatmap(giotto_SC, expression_values = "normalized", metadata_cols = 'fetalliver_labels',
                    selected_feats = markergenes_scran,
                    x_text_size = 18,
                    x_text_angle = 45,
                    y_text_size = 18,
                    strip_text_size = 18,  show_values = 'zscores_rescaled',
                    save_param = list(save_name = '5_c_metaheatmap'))

mycolorcode = c('plum', 'red', 'lavender','magenta',
                'purple', 'yellowgreen', 'yellow', 'gray', 'orange')
names(mycolorcode) = c('Neutrophil', 'Erythroid cell', 'Mast cell', 'Erythroid progenitor cell',
                       'Stem cell', 'C1q+ Macrophage', 'Megakaryocyte', 'Hepatocyte', 'Other Macrophage')
dimPlot2D(gobject = giotto_SC,     dim_reduction_name = 'umap',
          cell_color = "fetalliver_labels", cell_color_code = mycolorcode, show_NN_network = F, label_size = 8, point_size = 3.5, legend_text = 16,
          legend_symbol_size = 6, axis_text = 16, axis_title = 16,
          save_param = list(save_name = "8_a_Annotation"))

dimPlot2D(gobject = giotto_SC,     dim_reduction_name = 'umap',
          cell_color = "fetalliver_labels", cell_color_code = mycolorcode, show_NN_network = F, label_size = 8, point_size = 3.5, show_legend = F, legend_text = 16,
          legend_symbol_size = 6, axis_text = 16, axis_title = 16,
          save_param = list(save_name = "8_b_Annotation"))

dimPlot2D(gobject = giotto_SC,     dim_reduction_name = 'umap',
          cell_color = "fetalliver_labels", cell_color_code = mycolorcode, show_NN_network = F, show_center_label = F, label_size = 8, point_size = 3.5, show_legend = F, legend_text = 16,
          legend_symbol_size = 6, axis_text = 16, axis_title = 16,
          save_param = list(save_name = "8_c_Annotation"))


GetExpression <- getExpression(
  giotto_SC,
  values = "normalized",
  spat_unit = NULL,
  feat_type = NULL,
  output = c("exprObj", "matrix"),
  set_defaults = TRUE
)

GE <- as.data.frame(as.matrix(GetExpression@exprMat))

LR_data = data.table::fread("~/Desktop/Visium/The second eight mouse samples/E14-5_fetal_liver/mouse_lr_pair.txt")

LR_data[, ligand_det := ifelse(ligand_gene_symbol %in% giotto_SC@feat_ID[['rna']], T, F)]
LR_data[, receptor_det := ifelse(receptor_gene_symbol %in% giotto_SC@feat_ID[['rna']], T, F)]
LR_data_det = LR_data[ligand_det == T & receptor_det == T]
select_ligands = LR_data_det$ligand_gene_symbol
select_receptors = LR_data_det$receptor_gene_symbol

## get statistical significance of gene pair expression changes based on expression ##
expr_only_scores = exprCellCellcom(gobject = giotto_SC,
                                   cluster_column = 'fetalliver_labels',
                                   random_iter = 50,
                                   feat_set_1 = select_ligands,
                                   feat_set_2 = select_receptors)


## select top LR ##
selected_spat = expr_only_scores[p.adj <= 0.05 & abs(log2fc) > 0.1 & lig_nr >= 2 & rec_nr >= 2]
data.table::setorder(selected_spat, -PI)

top_LR_ints = unique(selected_spat[order(-abs(PI))]$LR_comb)[1:50]
top_LR_cell_ints = unique(selected_spat[order(-abs(PI))]$LR_cell_comb)[1:50]

plotCCcomHeatmap(gobject = giotto_SC,
                 comScores = expr_only_scores,
                 selected_LR = top_LR_ints,
                 selected_cell_LR = top_LR_cell_ints,
                 show = 'LR_expr')


plotCCcomDotplot(gobject = giotto_SC,
                 comScores = expr_only_scores,
                 selected_LR = top_LR_ints,
                 selected_cell_LR = top_LR_cell_ints,
                 cluster_on = 'PI')

library(writexl)

write_xlsx(GE, "~/Desktop/Referenced single cell-RNA sequencing data/E14.5 mouse fetal liver sc-RNA sequencing/3/E14.5FL expression matrix.xlsx")


write_xlsx(expr_only_scores, "~/Desktop/Referenced single cell-RNA sequencing data/E14.5 mouse fetal liver sc-RNA sequencing/3/E14.5FL ligand-receptor.xlsx")

