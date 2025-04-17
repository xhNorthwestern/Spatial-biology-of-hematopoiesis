library(Giotto)

genv_exists = checkGiottoEnvironment()

results_folder = '~/Desktop/Visium/Supplemental assays/Mouse/EDTA5D-3_BM/2'

instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  return_plot = FALSE)

data_path = '~/Desktop/Visium/Supplemental assays/Mouse/EDTA5D-3_BM'

visium_P8BM = createGiottoVisiumObject(visium_dir = data_path,
                                       expr_data = 'raw',
                                       png_name = 'tissue_lowres_image.png',
                                       gene_column_index = 2,
                                       instructions = instrs)

showGiottoImageNames(visium_P8BM)

pDataDT(visium_P8BM)

spatPlot2D(gobject = visium_P8BM, cell_color = 'in_tissue', point_size = 2,
           cell_color_code = c('0' = 'lightgrey', '1' = 'blue'),
           show_image = T, image_name = 'image')

metadata = pDataDT(visium_P8BM)
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
visium_P8BM = subsetGiotto(visium_P8BM, cell_ids = in_tissue_barcodes)


visium_P8BM <- filterGiotto(gobject = visium_P8BM,
                            expression_threshold = 1,
                            feat_det_in_min_cells = 50,
                            expression_values = c('raw'),
                            verbose = T)

visium_P8BM <- normalizeGiotto(gobject = visium_P8BM, scalefactor = 6000, verbose = T)


visium_P8BM <- addStatistics(gobject = visium_P8BM)


spatPlot2D(gobject = visium_P8BM, show_image = T, point_alpha = 0.7,
           cell_color = 'nr_feats', color_as_factor = F)

spatPlot2D(gobject = visium_P8BM, show_image = F, point_alpha = 0.7,
           cell_color = 'nr_feats', color_as_factor = F)


visium_P8BM <- calculateHVF(gobject = visium_P8BM, save_plot = TRUE)

gene_metadata = fDataDT(visium_P8BM)
featgenes = gene_metadata[hvf == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$feat_ID


visium_P8BM <- runPCA(gobject = visium_P8BM,
                      feats_to_use = featgenes)

screePlot(visium_P8BM, ncp = 30)

plotPCA(gobject = visium_P8BM)


visium_P8BM <- runUMAP(visium_P8BM, dimensions_to_use = 1:10)
plotUMAP(gobject = visium_P8BM)

visium_P8BM <- runtSNE(visium_P8BM, dimensions_to_use = 1:10)
plotTSNE(gobject = visium_P8BM)


visium_P8BM <- createNearestNetwork(gobject = visium_P8BM, dimensions_to_use = 1:10, k = 15)


visium_P8BM <- doLeidenCluster(gobject = visium_P8BM, resolution = 0.4, n_iterations = 1000)


dimPlot2D(gobject = visium_P8BM,     dim_reduction_name = 'umap',
          cell_color = "leiden_clus", show_NN_network = F, label_size = 10, point_size = 6.0, legend_text = 20,
          legend_symbol_size = 6, axis_text = 16, axis_title = 16,
          save_param = list(save_name = "7_Annotation"))

dimPlot2D(gobject = visium_P8BM,     dim_reduction_name = 'umap',
          cell_color = "leiden_clus", show_NN_network = F, label_size = 10, point_size = 6.0, show_legend = F, legend_text = 20,
          legend_symbol_size = 6, axis_text = 16, axis_title = 16,
          save_param = list(save_name = "7_a_Annotation"))

dimPlot2D(gobject = visium_P8BM,     dim_reduction_name = 'umap',
          cell_color = "leiden_clus", show_NN_network = F, label_size = 10, show_center_label = F, point_size = 6.0, show_legend = F, legend_text = 20,
          legend_symbol_size = 6, axis_text = 16, axis_title = 16,
          save_param = list(save_name = "7_b_Annotation"))


spatDimPlot(gobject = visium_P8BM, cell_color = 'leiden_clus',
            dim_point_size = 2, spat_point_size = 2.5)

spatDimPlot(gobject = visium_P8BM, cell_color = 'nr_feats', color_as_factor = F,
            dim_point_size = 2, spat_point_size = 2.5)


spatPlot2D(visium_P8BM, cell_color = 'leiden_clus', label_size = 10,
           coord_fix_ratio = 1, point_size = 5.4, show_legend = F, axis_text = 16,
           axis_title = 16)

spatPlot2D(visium_P8BM, cell_color = 'leiden_clus',
           group_by = 'leiden_clus', coord_fix_ratio = 1,
           cow_n_col = 4, show_legend = F,
           save_param = list(base_width = 14, base_height = 14))


gini_markers_subclusters = findMarkers_one_vs_all(gobject = visium_P8BM,
                                                  method = 'gini',
                                                  expression_values = 'normalized',
                                                  cluster_column = 'leiden_clus',
                                                  min_feats = 20,
                                                  min_expr_gini_score = 0.5,
                                                  min_det_gini_score = 0.5)
topgenes_gini = gini_markers_subclusters[, head(.SD, 2), by = 'cluster']$feats


violinPlot(visium_P8BM, feats = unique(topgenes_gini), cluster_column = 'leiden_clus',
           strip_text = 8, strip_position = 'right',
           save_param = list(base_width = 5, base_height = 10))


plotMetaDataHeatmap(visium_P8BM, selected_feats = unique(topgenes_gini),
                    metadata_cols = c('leiden_clus'),
                    x_text_size = 10, y_text_size = 10)


dimFeatPlot2D(visium_P8BM, expression_values = 'scaled',
              feats = gini_markers_subclusters[, head(.SD, 1), by = 'cluster']$feats,
              cow_n_col = 4, point_size = 0.75,
              save_param = list(base_width = 8, base_height = 8))


scran_markers_subclusters = findMarkers_one_vs_all(gobject = visium_P8BM,
                                                   method = 'scran',
                                                   expression_values = 'normalized',
                                                   cluster_column = 'leiden_clus')
topgenes_scran = scran_markers_subclusters[, head(.SD, 2), by = 'cluster']$feats


violinPlot(visium_P8BM, feats = unique(topgenes_scran), cluster_column = 'leiden_clus',
           strip_text = 10, strip_position = 'right',
           save_param = list(base_width = 5))


plotMetaDataHeatmap(visium_P8BM,  x_text_size = 20,
                    x_text_angle = 45,
                    y_text_size = 20,
                    strip_text_size = 20, selected_feats = topgenes_scran [c(1:12, 14:16, 18:22)],
                    metadata_cols = c('leiden_clus'))


dimFeatPlot2D(visium_P8BM, expression_values = 'scaled',
              feats = scran_markers_subclusters[, head(.SD, 1), by = 'cluster']$feats,
              cow_n_col = 3, point_size = 1.75,
              save_param = list(base_width = 8, base_height = 8))


Neutrophil1_markers = c("Chil3", "Camp", "Ngp", "S100a8", "Fcnb", "Hmgn2", "S100a9", "Ube2c", "Cebpe", "Lcn2", "Orm1", "Ltf", "Wfdc21", "Nusap1", "Spc25", "Hmgb2", "Tuba1c", "Asf1b", "Rflnb", "Serpinb1a")

Neutrophil2_markers = c("Mmp8", "Ifitm6", "Ltf", "Mmp9", "Lcn2", "Ly6g", "AA467197", "Wfdc21", "Ceacam10", "Slfn4", "Fpr2", "Anxa1", "Plbd1", "S100a9", "I830127L07Rik", "Cd177", "S100a6", "Retnlg", "Adpgk", "Chil1")

Myeloidcell_markers = c("Retnlg", "Il1b", "S100a6", "Ccl6", "Clec4d", "Cxcl2", "Junb", "Mmp8", "Cxcr2", "Fpr1", "H2-Q10", "Grina", "Dusp1", "Csf3r", "Mmp9", "Slpi", "Lilr4b", "Ccr1", "Slfn4", "S100a11")

Neutrophil3_markers = c("Elane", "Mpo", "Prtn3", "Ctsg", "Srgn", "Nkg7", "Ms4a3", "Gstm1", "Cst7", "Prss57", "Plac8", "Rgcc", "Mt1", "Calr", "Cd63", "Dmkn", "P4hb", "Mgl2", "Pdia6", "Vat1")

Macrophage1_markers = c("Ms4a6c", "F13a1", "Crip", "S100a4", "S100a10", "Ifi30", "Lgals1", "Fn1", "Ly6e", "Ly6c2", "Psap", "Ccr2", "Ass1", "Tmsb10", "Cst3", "C1galt1c1", "Prdx4", "Emb", "Irf8", "Ctsc")

Neutrophil4_markers = c("Ltf", "Lcn2", "St3gal5", "Ngp", "AA467197", "Ifitm6", "Camp", "Lyz2", "Zmpste24", "Chil3", "Wfdc21", "Orm1", "Anxa1", "Cybb", "Plbd1", "S100a8", "Ffar2", "Ltb4r1", "I830127L07Rik", "Itgb2l")

PreBcell_markers = c("Vpreb3", "Chchd10", "Igll1", "Ebf1", "Myl4", "Ighm", "Vpreb1", "Cd79a", "Mzb1", "Ptprcap", "Pafah1b3", "Blnk", "Fcrla", "Cd79b", "Cnp", "Ptma", "Smarca4", "Pou2af1", "Tifa", "Spib")

HSPC_markers = c("Cd34", "H2afy", "Mif", "Adgrg1", "Npm1", "Tmem176b", "Cdk6", "Rpl13", "Rpl3", "Rpl18a", "Gas5", "Plac8", "Rps5", "Rps3a1", "Rpl14", "Rps24", "Rpl11", "Rpl32", "Rps26", "Rpl18")

Erythroidcell_markers = c("Hba-a1", "Hbb-bs", "Hbb-bt", "Car2", "Hba-a2", "Car1", "Blvrb", "Prdx2", "Ctse", "Rhd", "Hmbs", "Aqp1", "Alad", "Klf1", "Cldn13", "Hebp1", "Gypa", "Glrx5", "Tmem14c", "Isg20")

Neutrophil5_markers = c("Fcnb", "Hmgn2", "Camp", "Serpinb1a", "Elane", "Chil3", "Cd63", "Rgcc", "Tmsb4x", "S100a8", "Ffar2", "Mogat2", "Cebpe", "Hsd11b1", "Lta4h", "Tuba8", "S100a9", "Tmem216", "Ms4a3", "Clec12a")

Bcell_markers = c("Igkc", "Igha", "Jchain", "Iglv1", "Iglc1", "Ighm", "Cd74", "Iglc2", "Ms4a1", "H2-Eb1", "Ly6d", "Iglc3", "H2-Aa", "Cd79a", "H2-Ab1", "Fcmr", "Cd79b", "Ighd", "Cd19", "Mzb1")

Macrophage2_markers = c("Ctss", "S100a4", "Ifitm3", "Wfdc17", "Apoe", "Clec4a3", "Ccr2", "Ms4a4c", "Fn1", "Clec4a1", "Lyz2", "Gngt2", "Ccl9", "Ms4a6c", "Ifi27l2a", "Psap", "Ms4a6b", "Cd68", "Atf3", "Crip1")

Tcell1_markers = c("Ccl5", "Ms4a4b", "Trbc2", "AW112010", "Trbc1", "Thy1", "Lck", "Gimap4", "Ctsw", "Gimap3", "Cd3g", "Cd3d", "Trac", "H2-Q7", "Klrd1", "Il2rb", "Lat", "Shisa5", "Rps27", "Gimap6")

Dendriticcell_markers = c("Siglech", "Cox6a2", "Bst2", "Ctsl", "Ly6d", "Irf8", "Cd7", "Klk1", "Tcf4", "H2-Aa", "Rnase6", "Cd74", "Ccr9", "Mpeg1", "Tsc22d1", "Ly6a", "Pltp", "Pld4", "Spib", "Upb1")

Basophil_markers = c("Prss34", "Mcpt8", "Cpa3", "Fcer1a", "Ccl3", "Ccl4", "Ms4a2", "Ier3", "Cd200r3", "Ifitm1", "Gata2", "Csrp3", "Rgs1", "Alox15", "Itga2b", "Cd63", "Cyp11a1", "Rnase12", "Stx3", "Ccl9")

Macrophage3_markers = c("Fcna", "C1qc", "C1qb", "C1qa", "Hmox1", "Vcam1", "Apoe", "Selenop", "Cd5l", "Pld3", "Axl", "Ftl1", "Igf1", "Ctsb", "Slc40a1", "Trf", "H2-Eb1", "H2-Aa", "Cd68", "Sdc3")

Neutrophil6_markers = c("Mpo", "Ctsg", "Prtn3", "Plac8", "Elane", "Cst7", "Srgn", "Nkg7", "Ms4a3", "Calr", "Gstm1", "Rps2", "Rgcc", "Slc25a5", "Ssr2", "Rps12", "Rps8", "Rps23", "Rpl23", "Rpl32")


PAGE_matrix_1 = makeSignMatrixPAGE(sign_names = c('Neutrophil_Chil3 high',
                                                  'Neutrophil_Mmp8 high',
                                                  'Myeloid Cell',
                                                  'Neutrophil_Elane high',
                                                  'Macrophage_Ms4a6c high',
                                                  'Neutrophil_Ltf high',
                                                  'Pre B Cell',
                                                  'Hematopoietic Stem and Progenitor Cell',
                                                  'Erythroid Cell',
                                                  'Neutrophil_Fcnb high',
                                                  'B Cell',
                                                  'Macrophage_Ctss high',
                                                  'T Cell_Ccl5 high',
                                                  'Plasmacytoid Dendritic Cell',
                                                  'Basophil',
                                                  'C1q+ Macrophage',
                                                  'Neutrophil_Mpo high'),
                                   sign_list = list(Neutrophil1_markers,
                                                    Neutrophil2_markers,
                                                    Myeloidcell_markers,
                                                    Neutrophil3_markers,
                                                    Macrophage1_markers,
                                                    Neutrophil4_markers,
                                                    PreBcell_markers,
                                                    HSPC_markers,
                                                    Erythroidcell_markers,
                                                    Neutrophil5_markers,
                                                    Bcell_markers,
                                                    Macrophage2_markers,
                                                    Tcell1_markers,
                                                    Dendriticcell_markers,
                                                    Basophil_markers,
                                                    Macrophage3_markers,
                                                    Neutrophil6_markers))

visium_P8BM = runPAGEEnrich(gobject = visium_P8BM, sign_matrix = PAGE_matrix_1)

cell_types_PAGE = colnames(PAGE_matrix_1)
plotMetaDataCellsHeatmap(gobject = visium_P8BM,
                         metadata_cols = 'leiden_clus',
                         value_cols = cell_types_PAGE,
                         spat_enr_names = 'PAGE',
                         x_text_size = 18,
                         y_text_size = 20)


spatCellPlot2D(gobject = visium_P8BM,
               spat_enr_names = 'PAGE',
               cell_annotation_values = cell_types_PAGE[13:16],
               cow_n_col = 2,coord_fix_ratio = 1, point_size = 3, show_legend = T)


spatDimCellPlot2D(gobject = visium_P8BM,
                  spat_enr_names = 'PAGE',
                  cell_annotation_values = cell_types_PAGE[13:16],
                  cow_n_col = 1, spat_point_size = 1,
                  plot_alignment = 'horizontal',
                  save_param = list(base_width=7, base_height=10))


Spatialcorrelation <- getSpatialEnrichment(
  visium_P8BM,
  spat_unit = NULL,
  feat_type = NULL,
  name = "PAGE",
  output = c("spatEnrObj", "data.table"),
  copy_obj = TRUE,
  set_defaults = TRUE
)

library(writexl)

write_xlsx(Spatialcorrelation@enrichDT, "~/Desktop/Visium/Supplemental assays/Mouse/EDTA5D-3_BM/2/spatial enrichment.xlsx")



















