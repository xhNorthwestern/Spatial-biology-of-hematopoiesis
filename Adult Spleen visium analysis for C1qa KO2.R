results_folder = 'C:/Users/ser/Desktop/Visium/Spleen/KO2_spleen/KO2_spleen/2'

library(Giotto)

python_path <- 'C:/Users/ser/anaconda3/envs/my_python_env/python.exe'

# 3. Create Giotto instructions
# Directly saving plots to the working directory without rendering them in the editor saves time.
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  return_plot = FALSE,
                                  python_path = python_path)

data_path = 'C:/Users/ser/Desktop/Visium/Spleen/KO2_spleen/KO2_spleen'

visium_Spleen = createGiottoVisiumObject(visium_dir = data_path,
                                        expr_data = 'raw',
                                        png_name = 'tissue_lowres_image.png',
                                        gene_column_index = 2,
                                        instructions = instrs)

showGiottoImageNames(visium_Spleen)

pDataDT(visium_Spleen)

spatPlot2D(gobject = visium_Spleen, cell_color = 'in_tissue', point_size = 2,
           cell_color_code = c('0' = 'lightgrey', '1' = 'blue'),
           show_image = T, image_name = 'image')

metadata = pDataDT(visium_Spleen)
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
visium_Spleen = subsetGiotto(visium_Spleen, cell_ids = in_tissue_barcodes)


visium_Spleen <- filterGiotto(gobject = visium_Spleen,
                             expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000,
                             expression_values = c('raw'),
                             verbose = T)

visium_Spleen <- normalizeGiotto(gobject = visium_Spleen, scalefactor = 6000, verbose = T)


visium_Spleen <- addStatistics(gobject = visium_Spleen)


spatPlot2D(gobject = visium_Spleen, show_image = T, point_alpha = 1.5,
           cell_color = 'nr_feats', color_as_factor = F)


visium_Spleen <- calculateHVF(gobject = visium_Spleen, save_plot = TRUE)

gene_metadata = fDataDT(visium_Spleen)
featgenes = gene_metadata[hvf == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$feat_ID


visium_Spleen <- runPCA(gobject = visium_Spleen,
                       feats_to_use = featgenes)

screePlot(visium_Spleen, ncp = 30)

plotPCA(gobject = visium_Spleen)


visium_Spleen <- runUMAP(visium_Spleen, dimensions_to_use = 1:10)
plotUMAP(gobject = visium_Spleen)

visium_Spleen <- runtSNE(visium_Spleen, dimensions_to_use = 1:10)
plotTSNE(gobject = visium_Spleen)


visium_Spleen <- createNearestNetwork(gobject = visium_Spleen, dimensions_to_use = 1:10, k = 15)


visium_Spleen <- doLeidenCluster(gobject = visium_Spleen, resolution = 0.4, n_iterations = 1000)

dimPlot2D(gobject = visium_Spleen,     dim_reduction_name = 'umap',
          cell_color = "leiden_clus", show_NN_network = F, label_size = 10, point_size = 6.0, legend_text = 20,
          legend_symbol_size = 6, axis_text = 16, axis_title = 16,
          save_param = list(save_name = "7_Annotation"))

dimPlot2D(gobject = visium_Spleen,     dim_reduction_name = 'umap',
          cell_color = "leiden_clus", show_NN_network = F, label_size = 10, point_size = 6.0, show_legend = F, legend_text = 20,
          legend_symbol_size = 6, axis_text = 16, axis_title = 16,
          save_param = list(save_name = "7_a_Annotation"))

dimPlot2D(gobject = visium_Spleen,     dim_reduction_name = 'umap',
          cell_color = "leiden_clus", show_NN_network = F, label_size = 10, show_center_label = F, point_size = 6.0, show_legend = F, legend_text = 20,
          legend_symbol_size = 6, axis_text = 16, axis_title = 16,
          save_param = list(save_name = "7_b_Annotation"))

spatDimPlot(gobject = visium_Spleen, cell_color = 'leiden_clus',
            dim_point_size = 2, spat_point_size = 2.5)

spatDimPlot(gobject = visium_Spleen, cell_color = 'nr_feats', color_as_factor = F,
            dim_point_size = 2, spat_point_size = 2.5)


spatPlot2D(visium_Spleen, cell_color = 'leiden_clus', label_size = 10,
           coord_fix_ratio = 1, point_size = 5.5, show_legend = F, axis_text = 16,
           axis_title = 16)

spatPlot2D(visium_Spleen, cell_color = 'leiden_clus',
           group_by = 'leiden_clus', coord_fix_ratio = 1,
           cow_n_col = 4, show_legend = F,
           save_param = list(base_width = 14, base_height = 14))


gini_markers_subclusters = findMarkers_one_vs_all(gobject = visium_Spleen,
                                                  method = 'gini',
                                                  expression_values = 'normalized',
                                                  cluster_column = 'leiden_clus',
                                                  min_feats = 20,
                                                  min_expr_gini_score = 0.5,
                                                  min_det_gini_score = 0.5)
topgenes_gini = gini_markers_subclusters[, head(.SD, 2), by = 'cluster']$feats


violinPlot(visium_Spleen, feats = unique(topgenes_gini), cluster_column = 'leiden_clus',
           strip_text = 8, strip_position = 'right',
           save_param = list(base_width = 5, base_height = 10))


plotMetaDataHeatmap(visium_Spleen, selected_feats = unique(topgenes_gini),
                    metadata_cols = c('leiden_clus'),
                    x_text_size = 10, y_text_size = 10)


dimFeatPlot2D(visium_Spleen, expression_values = 'scaled',
              feats = gini_markers_subclusters[, head(.SD, 1), by = 'cluster']$feats,
              cow_n_col = 3, point_size = 1.75,
              save_param = list(base_width = 8, base_height = 8))


scran_markers_subclusters = findMarkers_one_vs_all(gobject = visium_Spleen,
                                                   method = 'scran',
                                                   expression_values = 'normalized',
                                                   cluster_column = 'leiden_clus')
topgenes_scran = scran_markers_subclusters[, head(.SD, 2), by = 'cluster']$feats


violinPlot(visium_Spleen, feats = unique(topgenes_scran), cluster_column = 'leiden_clus',
           strip_text = 10, strip_position = 'right',
           save_param = list(base_width = 5))


plotMetaDataHeatmap(visium_Spleen,  x_text_size = 20,
                    x_text_angle = 45,
                    y_text_size = 20,
                    strip_text_size = 20, selected_feats = topgenes_scran,
                    metadata_cols = c('leiden_clus'))


dimFeatPlot2D(visium_Spleen, expression_values = 'scaled',
              feats = scran_markers_subclusters[, head(.SD, 1), by = 'cluster']$feats,
              cow_n_col = 3, point_size = 1.75,
              save_param = list(base_width = 8, base_height = 8))

Bcell1_markers = c("Ighd", "Ly6d", "Cd79a", "Iglc3", "Ms4a1", "RP24-338A5.4", "Cd19", "Spib", "Fcmr", "Fcer2a", "Cd79b", "Cd74", "Igkc", "Siglecg", "Ebf1", "Serpinb1a", "Snn", "Hvcn1", "Iglc1", "Blnk")

Tcell_markers = c("Trac", "Trbc2", "Cd8b1", "Cd3g", "Ms4a4b", "Cd8a", "Lat", "Saraf", "Thy1", "Trbc1", "Tcf7", "Ms4a6b", "Dapl1", "Cd247", "Cd27", "Atp1b3", "Prkcq", "Lck", "Cd3d", "Actn1")

NKcell_markers = c("Gzma", "Klra8", "Klra13-ps", "Klrk1", "Klre1", "Ccl5", "Ncr1", "Klrb1c", "Klra9", "Klra4", "Nkg7", "Klra7", "Klrd1", "Klrc2", "AW112010", "Ctsw", "Xcl1", "Klri2", "Klrc1", "Il2rb")

DendriticCell1_markers = c("S100a4", "Ffar2", "Plbd1", "Tmem176b", "Cst3", "Rgs2", "Avpi1", "Bcl2a1b", "Cyp4f16", "Tmem176a", "Asb2", "Rogdi", "Itgax", "Bcl2a1d", "H2-Eb1", "Flt3", "I830077J02Rik", "Dnase1l3", "Gm2a", "Ccnd1")

Monocyte_markers = c("Lyz2", "Ifitm3", "S100a4", "Ms4a6c", "Lgals3", "Ccr2", "Ccl9", "Gm9733", "Ccl6", "Alox5ap", "Fcer1g", "Wfdc17", "Cst3", "Clec4a3", "Psap", "F13a1", "Clec4a1", "Apoe", "Anxa2", "Sat1")

Neutrophil1_markers = c("S100a8", "S100a9", "Retnlg", "Il1b", "Ccl6", "Msrb1", "Mmp9", "Hp", "Slpi", "Csf3r", "S100a11", "Clec4d", "Cd300ld", "Junb", "Gsr", "Cxcr2", "Anxa2", "Fos", "Ifitm1", "Pla2g7")

Macrophage_markers = c("C1qc", "Hmox1", "Fcna",  "C1qa", "Vcam1", "C1qb", "Trf", "Slc40a1", "Ctsb", "Cfp", "Slpi", "Creg1", "Axl", "Apoe", "Ftl1-ps1", "Ftl1", "Spic", "Csf1r", "Itgad", "Maf")

Bcell2_markers = c("Jchain", "Iglv2", "Igkv5-48", "Ighm", "Igha", "Iglc2", "Iglv1", "Igkc", "Igkv1-117", "Igkv3-2", "Igkv16-104", "Ighg2c", "Iglc1", "Mzb1", "Hsp90b1", "Edem1", "Txndc5", "Fkbp11", "Creld2", "Xbp1")

DendriticCell2_markers = c("Siglech", "Bst2", "Ctsl", "Klk1", "Cox6a2", "Cd209d", "Irf8", "Cd209a", "Mpeg1", "Klra17", "Pltp", "Atp1b1", "Klk1b27", "Grn", "Upb1", "Gm21762", "Tyrobp", "Ccr9", "Mvb12a", "Smim5")

ErythroidCell_markers = c("Hba-a1", "Car2", "Hbb-bs", "Hbb-bt", "Car1", "Hba-a2", "Klf1", "Hmbs", "Prdx2", "Mt1", "Fam132a", "Atpif1", "Ube2c", "Rhd", "Blvrb", "Birc5", "Ctse", "Tuba1b", "Gypa", "Asns")

Neutrophil2_markers = c("Ngp", "Camp", "S100a9", "S100a8", "Lcn2", "Ltf", "Chil3", "Wfdc21", "Anxa1", "Fcnb", "I830127L07Rik", "Lyz2", "Cd177", "Hmgn2", "Mmp8", "Ly6g", "Elane", "Lrg1", "Pglyrp1", "Hp")

PAGE_matrix_1 = makeSignMatrixPAGE(sign_names = c('Marginal zone B cell',
                                                  'T cell',
                                                  'NK cell',
                                                  'Dendritic cell_S100a4 high',
                                                  'Monocyte',
                                                  'Neutrophil_S100a8 high',
                                                  'C1q+ Macrophage',
                                                  'B cell',
                                                  'Plasmacytoid dendritic cell_Siglech high',
                                                  'Erythroid cell',
                                                  'Neutrophil_Ngp high'),
                                   sign_list = list(Bcell1_markers,
                                                    Tcell_markers,
                                                    NKcell_markers,
                                                    DendriticCell1_markers,
                                                    Monocyte_markers,
                                                    Neutrophil1_markers,
                                                    Macrophage_markers,
                                                    Bcell2_markers,
                                                    DendriticCell2_markers,
                                                    ErythroidCell_markers,
                                                    Neutrophil2_markers))

visium_Spleen = runPAGEEnrich(gobject = visium_Spleen, sign_matrix = PAGE_matrix_1)

cell_types_PAGE = colnames(PAGE_matrix_1)
plotMetaDataCellsHeatmap(gobject = visium_Spleen,
                         metadata_cols = 'leiden_clus',
                         value_cols = cell_types_PAGE,
                         spat_enr_names = 'PAGE',
                         x_text_size = 20,
                         y_text_size = 20)


spatCellPlot2D(gobject = visium_Spleen,
               spat_enr_names = 'PAGE',
               cell_annotation_values = cell_types_PAGE[8:11],
               cow_n_col = 2,coord_fix_ratio = 1, point_size = 3, show_legend = T)


spatDimCellPlot2D(gobject = visium_Spleen,
                  spat_enr_names = 'PAGE',
                  cell_annotation_values = cell_types_PAGE[1:4],
                  cow_n_col = 1, spat_point_size = 1.5,
                  plot_alignment = 'horizontal',
                  save_param = list(base_width=7, base_height=10))


Spatialcorrelation <- getSpatialEnrichment(
  visium_Spleen,
  spat_unit = NULL,
  feat_type = NULL,
  name = "PAGE",
  output = c("spatEnrObj", "data.table"),
  copy_obj = TRUE,
  set_defaults = TRUE
)


write.xlsx(Spatialcorrelation@enrichDT, file = "myworkbookforKO2Spleen.xlsx",
           sheetName = "Spatialcorrelation@enrichDT", append = FALSE)
