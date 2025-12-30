results_folder = 'C:/Users/ser/Desktop/Visium/WT4_FetalLiver/WT4_FetalLiver/2'

library(Giotto)

python_path <- 'C:/Users/ser/anaconda3/envs/my_python_env/python.exe'

# 3. Create Giotto instructions
# Directly saving plots to the working directory without rendering them in the editor saves time.
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  return_plot = FALSE,
                                  python_path = python_path)

data_path = 'C:/Users/ser/Desktop/Visium/WT4_FetalLiver/WT4_FetalLiver'

visium_E14.5FetalLiver = createGiottoVisiumObject(visium_dir = data_path,
                                                  expr_data = 'raw',
                                                  png_name = 'tissue_lowres_image.png',
                                                  gene_column_index = 2,
                                                  instructions = instrs)

showGiottoImageNames(visium_E14.5FetalLiver)

pDataDT(visium_E14.5FetalLiver)

spatPlot2D(gobject = visium_E14.5FetalLiver, cell_color = 'in_tissue', point_size = 2,
           cell_color_code = c('0' = 'lightgrey', '1' = 'blue'),
           show_image = T, image_name = 'image')

metadata = pDataDT(visium_E14.5FetalLiver)
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
visium_E14.5FetalLiver = subsetGiotto(visium_E14.5FetalLiver, cell_ids = in_tissue_barcodes)


visium_E14.5FetalLiver <- filterGiotto(gobject = visium_E14.5FetalLiver,
                                       expression_threshold = 1,
                                       feat_det_in_min_cells = 50,
                                       min_det_feats_per_cell = 1000,
                                       expression_values = c('raw'),
                                       verbose = T)

visium_E14.5FetalLiver <- normalizeGiotto(gobject = visium_E14.5FetalLiver, scalefactor = 6000, verbose = T)


visium_E14.5FetalLiver <- addStatistics(gobject = visium_E14.5FetalLiver)


spatPlot2D(gobject = visium_E14.5FetalLiver, show_image = T, point_alpha = 1,
           cell_color = 'nr_feats', color_as_factor = F)


visium_E14.5FetalLiver <- calculateHVF(gobject = visium_E14.5FetalLiver, save_plot = TRUE)

gene_metadata = fDataDT(visium_E14.5FetalLiver)
featgenes = gene_metadata[hvf == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$feat_ID

visium_E14.5FetalLiver <- runPCA(gobject = visium_E14.5FetalLiver, feats_to_use = featgenes)

screePlot(visium_E14.5FetalLiver, ncp = 30)

dimPlot2D(gobject = visium_E14.5FetalLiver,
          dim_reduction_to_use = "pca")


visium_E14.5FetalLiver <- runUMAP(visium_E14.5FetalLiver, dimensions_to_use = 1:10)
plotUMAP(gobject = visium_E14.5FetalLiver)

visium_E14.5FetalLiver <- runtSNE(visium_E14.5FetalLiver, dimensions_to_use = 1:10)
plotTSNE(gobject = visium_E14.5FetalLiver)


visium_E14.5FetalLiver <- createNearestNetwork(gobject = visium_E14.5FetalLiver, dimensions_to_use = 1:10, k = 15)


visium_E14.5FetalLiver <- doLeidenCluster(gobject = visium_E14.5FetalLiver, resolution = 0.5, n_iterations = 1000)

dimPlot2D(gobject = visium_E14.5FetalLiver,     dim_reduction_name = 'umap',
          cell_color = "leiden_clus", show_NN_network = F, label_size = 10, point_size = 6.0, legend_text = 20,
          legend_symbol_size = 6, axis_text = 16, axis_title = 16,
          save_param = list(save_name = "7_Annotation"))

dimPlot2D(gobject = visium_E14.5FetalLiver,     dim_reduction_name = 'umap',
          cell_color = "leiden_clus", show_NN_network = F, label_size = 10, point_size = 6.0, show_legend = F, legend_text = 20,
          legend_symbol_size = 6, axis_text = 16, axis_title = 16,
          save_param = list(save_name = "7_a_Annotation"))

dimPlot2D(gobject = visium_E14.5FetalLiver,     dim_reduction_name = 'umap',
          cell_color = "leiden_clus", show_NN_network = F, label_size = 10, show_center_label = F, point_size = 6.0, show_legend = F, legend_text = 20,
          legend_symbol_size = 6, axis_text = 16, axis_title = 16,
          save_param = list(save_name = "7_b_Annotation"))

spatDimPlot(gobject = visium_E14.5FetalLiver, cell_color = 'nr_feats', color_as_factor = F,
            dim_point_size = 2, spat_point_size = 2.5)


spatPlot2D(visium_E14.5FetalLiver, cell_color = 'leiden_clus', label_size = 10,
           coord_fix_ratio = 1, point_size = 6.5, show_legend = F, axis_text = 16,
           axis_title = 16)

spatPlot2D(visium_E14.5FetalLiver, cell_color = 'leiden_clus',
           group_by = 'leiden_clus', coord_fix_ratio = 1,
           cow_n_col = 3, show_legend = F,
           save_param = list(base_width = 14, base_height = 14))


gini_markers_subclusters = findMarkers_one_vs_all(gobject = visium_E14.5FetalLiver,
                                                  method = 'gini',
                                                  expression_values = 'normalized',
                                                  cluster_column = 'leiden_clus',
                                                  min_feats = 20,
                                                  min_expr_gini_score = 0.5,
                                                  min_det_gini_score = 0.5)
topgenes_gini = gini_markers_subclusters[, head(.SD, 2), by = 'cluster']$feats


violinPlot(visium_E14.5FetalLiver, feats = unique(topgenes_gini), cluster_column = 'leiden_clus',
           strip_text = 8, strip_position = 'right',
           save_param = list(base_width = 5, base_height = 10))


plotMetaDataHeatmap(visium_E14.5FetalLiver, selected_feats = unique(topgenes_gini),
                    metadata_cols = c('leiden_clus'),
                    x_text_size = 10, y_text_size = 10)


dimFeatPlot2D(visium_E14.5FetalLiver, expression_values = 'scaled',
              feats = gini_markers_subclusters[, head(.SD, 1), by = 'cluster']$feats,
              cow_n_col = 3, point_size = 2.5,
              save_param = list(base_width = 8, base_height = 8))


scran_markers_subclusters = findMarkers_one_vs_all(gobject = visium_E14.5FetalLiver,
                                                   method = 'scran',
                                                   expression_values = 'normalized',
                                                   cluster_column = 'leiden_clus')
topgenes_scran = scran_markers_subclusters[, head(.SD, 2), by = 'cluster']$feats


violinPlot(visium_E14.5FetalLiver, feats = unique(topgenes_scran), cluster_column = 'leiden_clus',
           strip_text_size = 8, strip_position = 'right',
           save_param = list(base_width = 5))


plotMetaDataHeatmap(visium_E14.5FetalLiver,  x_text_size = 20,
                    x_text_angle = 45,
                    y_text_size = 20,
                    strip_text_size = 20, selected_feats = topgenes_scran,
                    metadata_cols = c('leiden_clus'))


dimFeatPlot2D(visium_E14.5FetalLiver, expression_values = 'scaled',
              feats = scran_markers_subclusters[, head(.SD, 1), by = 'cluster']$feats,
              cow_n_col = 3, point_size = 2.5,
              save_param = list(base_width = 8, base_height = 8))

spatFeatPlot2D(visium_E14.5FetalLiver, expression_values = 'scaled',
               feats = scran_markers_subclusters[, head(.SD, 1), by = 'cluster']$feats,
               cow_n_col = 3, point_size = 2.5,
               save_param = list(base_width = 8, base_height = 8))

Car2highErythroidCell_markers = c("Klf1", "Car2", "Cited4", "Rfesd", "Asns", "Hebp1", "Alad", "Hmgb3", "Blvrb", "Hmbs", "Tmem14c", "Rhd", "Cpox", "Pla2g12a", "Prdx2", "Glrx5", "Icam4", "Cldn13", "Smim1", "Mns1")

TyrobphighMacrophage_markers = c("Lyz2", "Ctsc", "Ifitm3", "Lgals3", "Cd74", "Cd52", "Ms4a6c", "Id2", "F13a1", "Tmsb4x", "Coro1a", "Tyrobp", "Ctss", "Psap", "Fos", "Ly6e", "Anxa2", "Tmsb10", "Cst3", "Emb")

MpohighNeutriphil_markers = c("Elane", "Mpo", "Prtn3", "Ctsg", "Ms4a3", "Plac8", "Serpinb1a", "Fcnb", "Srgn", "Cd63", "Stfa1", "BC100530", "Nkg7", "S100a8", "Gstm1", "Arhgdib", "Tmed3", "Prss57", "Clec12a", "Cmtm7")

HbbbthighErythroidCell_markers = c("Hba-a2", "Hbb-bs", "Hbb-bt", "Hba-a1", "Alas2", "Slc25a37", "Gypa", "Slc4a1", "Pnpo", "Snca", "Fam210b", "Cdc25b", "Ackr1", "Cd59a", "Ube2l6", "Cd24a", "Hbp1b", "Fech", "Gmpr", "Mgst3")

C1qbhighMacrophage_markers = c("Hmox1", "C1qc", "C1qb", "Fcna", "Cd5l", "Cfp", "Lgmn", "C1qa", "Slc40a1", "Grn", "Ctsd", "Ctsb", "Apoe", "Ccl24", "Ctss", "Marco", "Adgre1", "Aif1", "Selenop", "Il18bp")

Mastcell_markers = c("Mcpt8", "Cpa3", "Ifitm1", "Prss34",  "Hdc", "Cd63", "Cd200r3", "Alox5", "Alox5ap", "Ccl3", "Ccl9", "Gata2", "Cd9", "Cyp4f18", "Anxa2", "Fcgr3", "Spp1", "Lilr4b", "Perp", "Tmsb4x")

Hepatocyte_markers = c("Afp", "Alb", "Apoa2", "Ttr", "Rbp4", "Apoa1", "Trf", "Serpina6", "Fgb", "Apom", "Kng1", "Ambp", "Ahsg", "Fgg", "Serpina1b", "H19", "Dlk1", "Serpina1a", "Fabp1", "Serpinf2")

Megakaryocyte_markers = c("Pf4", "Ppbp", "Gp1bb",  "Rap1b", "Tpm4", "Lat", "Bin2", "Sdpr", "Zyx", "Cd9", "Spns1", "Itga2b", "Gp9", "F2rl2", "Treml1", "Rsu1", "Plek", "Tuba8", "Timp3", "Rabgap1l")

StemCell_markers = c("H2afy", "Cd34", "Cmtm7",  "Ptprcap", "Marcksl1", "Psmb8", "Ramp1", "Emb", "Idh2", "BC035044", "Gmfg", "Tmsb10", "Celf2", "Ctla2a", "Ptpn18", "Bmyc", "Ccnd2", "Coro1a", "Mef2c", "Prtn3")


PAGE_matrix_1 = makeSignMatrixPAGE(sign_names = c('Erythroid progenitor cell',
                                                  'Other Macrophage',
                                                  'Neutrophil',
                                                  'Erythroid cell',
                                                  'C1q+ Macrophage',
                                                  'Mast cell',
                                                  'Hepatocyte',
                                                  'Megakaryocyte',
                                                  'Stem cell'),
                                   sign_list = list(Car2highErythroidCell_markers,
                                                    TyrobphighMacrophage_markers,
                                                    MpohighNeutriphil_markers,
                                                    HbbbthighErythroidCell_markers,
                                                    C1qbhighMacrophage_markers,
                                                    Mastcell_markers,
                                                    Hepatocyte_markers,
                                                    Megakaryocyte_markers,
                                                    StemCell_markers))

visium_E14.5FetalLiver = runPAGEEnrich(gobject = visium_E14.5FetalLiver, sign_matrix = PAGE_matrix_1)

cell_types_PAGE = colnames(PAGE_matrix_1)
plotMetaDataCellsHeatmap(gobject = visium_E14.5FetalLiver,
                         metadata_cols = 'leiden_clus',
                         value_cols = cell_types_PAGE,
                         spat_enr_names = 'PAGE',
                         x_text_size = 20,
                         y_text_size = 20)

spatCellPlot2D(gobject = visium_E14.5FetalLiver,
               spat_enr_names = 'PAGE',
               cell_annotation_values = cell_types_PAGE[1:4],
               cow_n_col = 2,coord_fix_ratio = 1, point_size = 3, show_legend = T)


Spatialcorrelation <- getSpatialEnrichment(
  visium_E14.5FetalLiver,
  spat_unit = NULL,
  feat_type = NULL,
  name = "PAGE",
  output = c("spatEnrObj", "data.table"),
  copy_obj = TRUE,
  set_defaults = TRUE
)

write.xlsx(Spatialcorrelation@enrichDT, file = "myworkbookforWT4.xlsx",
           sheetName = "Spatialcorrelation@enrichDT", append = FALSE)





