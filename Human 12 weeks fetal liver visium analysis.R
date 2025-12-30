
library(Giotto)

genv_exists = checkGiottoEnvironment()

results_folder = '~/Desktop/Visium/The third twelve human fetus samples/Matrix/The first eight samples/NMS21-13027_A4/2'

instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  python_path = '/Users/xuhan/Library/r-miniconda-arm64/envs/giotto_env/bin/python')

data_path = '~/Desktop/Visium/The third twelve human fetus samples/Matrix/The first eight samples/NMS21-13027_A4'

visium_human12weeksFetalLiver = createGiottoVisiumObject(visium_dir = data_path,
                                                         expr_data = 'raw',
                                                         png_name = 'tissue_lowres_image.png',
                                                         gene_column_index = 2,
                                                         instructions = instrs)

showGiottoImageNames(visium_human12weeksFetalLiver)

pDataDT(visium_human12weeksFetalLiver)

spatPlot2D(gobject = visium_human12weeksFetalLiver, cell_color = 'in_tissue', point_size = 2,
           cell_color_code = c('0' = 'lightgrey', '1' = 'blue'),
           show_image = T, image_name = 'image')

metadata = pDataDT(visium_human12weeksFetalLiver)
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
visium_human12weeksFetalLiver = subsetGiotto(visium_human12weeksFetalLiver, cell_ids = in_tissue_barcodes)


visium_human12weeksFetalLiver <- filterGiotto(gobject = visium_human12weeksFetalLiver,
                                              expression_threshold = 1,
                                              feat_det_in_min_cells = 50,
                                              min_det_feats_per_cell = 1000,
                                              expression_values = c('raw'),
                                              verbose = T)

visium_human12weeksFetalLiver <- normalizeGiotto(gobject = visium_human12weeksFetalLiver, scalefactor = 6000, verbose = T)


visium_human12weeksFetalLiver <- addStatistics(gobject = visium_human12weeksFetalLiver)


spatPlot2D(gobject = visium_human12weeksFetalLiver, show_image = T, point_alpha = 0.7,
           cell_color = 'nr_feats', color_as_factor = F)


visium_human12weeksFetalLiver <- calculateHVF(gobject = visium_human12weeksFetalLiver, save_plot = TRUE)

gene_metadata = fDataDT(visium_human12weeksFetalLiver)
featgenes = gene_metadata[hvf == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$feat_ID


visium_human12weeksFetalLiver <- runPCA(gobject = visium_human12weeksFetalLiver,
                                        genes_to_use = featgenes)

screePlot(visium_human12weeksFetalLiver, ncp = 30)

plotPCA(gobject = visium_human12weeksFetalLiver)


visium_human12weeksFetalLiver <- runUMAP(visium_human12weeksFetalLiver, dimensions_to_use = 1:10)
plotUMAP(gobject = visium_human12weeksFetalLiver)


visium_human12weeksFetalLiver <- createNearestNetwork(gobject = visium_human12weeksFetalLiver, dimensions_to_use = 1:10, k = 15)


visium_human12weeksFetalLiver <- doLeidenCluster(gobject = visium_human12weeksFetalLiver, resolution = 0.6, n_iterations = 1000)


dimPlot2D(gobject = visium_human12weeksFetalLiver,     dim_reduction_name = 'umap',
          cell_color = "leiden_clus", show_NN_network = F, label_size = 10, point_size = 6.0, legend_text = 20,
          legend_symbol_size = 6, axis_text = 16, axis_title = 16,
          save_param = list(save_name = "7_Annotation"))

dimPlot2D(gobject = visium_human12weeksFetalLiver,     dim_reduction_name = 'umap',
          cell_color = "leiden_clus", show_NN_network = F, label_size = 10, point_size = 6.0, show_legend = F, legend_text = 20,
          legend_symbol_size = 6, axis_text = 16, axis_title = 16,
          save_param = list(save_name = "7_a_Annotation"))

dimPlot2D(gobject = visium_human12weeksFetalLiver,     dim_reduction_name = 'umap',
          cell_color = "leiden_clus", show_NN_network = F, label_size = 10, show_center_label = F, point_size = 6.0, show_legend = F, legend_text = 20,
          legend_symbol_size = 6, axis_text = 16, axis_title = 16,
          save_param = list(save_name = "7_b_Annotation"))

dimPlot2D(gobject = visium_human12weeksFetalLiver,     dim_reduction_name = 'umap',
          cell_color = "leiden_clus", show_NN_network = F, label_size = 10, show_center_label = F, point_size = 5.0, show_legend = F, legend_text = 20,
          legend_symbol_size = 6, axis_text = 16, axis_title = 16,
          save_param = list(save_name = "7_c_Annotation"))

spatDimPlot(gobject = visium_human12weeksFetalLiver, cell_color = 'leiden_clus',
            dim_point_size = 2, spat_point_size = 2.5)

spatDimPlot(gobject = visium_human12weeksFetalLiver, cell_color = 'nr_feats', color_as_factor = F,
            dim_point_size = 2, spat_point_size = 2.5)

spatPlot2D(visium_human12weeksFetalLiver, cell_color = 'leiden_clus', label_size = 10,
           coord_fix_ratio = 1, point_size = 3.7, show_legend = F, axis_text = 16,
           axis_title = 16)

spatPlot2D(visium_human12weeksFetalLiver, cell_color = 'leiden_clus',
           group_by = 'leiden_clus', coord_fix_ratio = 1,
           cow_n_col = 4, show_legend = F,
           save_param = list(base_width = 14, base_height = 14))


gini_markers_subclusters = findMarkers_one_vs_all(gobject = visium_human12weeksFetalLiver,
                                                  method = 'gini',
                                                  expression_values = 'normalized',
                                                  cluster_column = 'leiden_clus',
                                                  min_feats = 20,
                                                  min_expr_gini_score = 0.5,
                                                  min_det_gini_score = 0.5)
topgenes_gini = gini_markers_subclusters[, head(.SD, 2), by = 'cluster']$feats


violinPlot(visium_human12weeksFetalLiver, feats = unique(topgenes_gini), cluster_column = 'leiden_clus',
           strip_text = 8, strip_position = 'right',
           save_param = list(base_width = 5, base_height = 10))


plotMetaDataHeatmap(visium_human12weeksFetalLiver, selected_feats = unique(topgenes_gini),
                    metadata_cols = c('leiden_clus'),
                    x_text_size = 10, y_text_size = 10)


dimFeatPlot2D(visium_human12weeksFetalLiver, expression_values = 'scaled',
              feats = gini_markers_subclusters[, head(.SD, 1), by = 'cluster']$feats,
              cow_n_col = 3, point_size = 1.75,
              save_param = list(base_width = 8, base_height = 8))


scran_markers_subclusters = findMarkers_one_vs_all(gobject = visium_human12weeksFetalLiver,
                                                   method = 'scran',
                                                   expression_values = 'normalized',
                                                   cluster_column = 'leiden_clus')
topgenes_scran = scran_markers_subclusters[, head(.SD, 2), by = 'cluster']$feats


violinPlot(visium_human12weeksFetalLiver, feats = unique(topgenes_scran), cluster_column = 'leiden_clus',
           strip_text = 10, strip_position = 'right',
           save_param = list(base_width = 5))


plotMetaDataHeatmap(visium_human12weeksFetalLiver,  x_text_size = 20,
                    x_text_angle = 45,
                    y_text_size = 20,
                    strip_text_size = 20, selected_feats = topgenes_scran,
                    metadata_cols = c('leiden_clus'))


dimFeatPlot2D(visium_human12weeksFetalLiver, expression_values = 'scaled',
              feats = scran_markers_subclusters[, head(.SD, 1), by = 'cluster']$feats,
              cow_n_col = 3, point_size = 1.75,
              save_param = list(base_width = 8, base_height = 8))

Bcell_markers = c("CD74", "IGHM", "TCL1A", "CD79B", "HLA-DRA", "CD79A", "CD52", "VPREB3", "IGLL1", "CD24", "IGKC", "ACTG1", "CD37", "CXCR4", "HLA-B", "LTB", "TMSB10", "DNAJB1", "MS4A1", "TMSB4X")

Tcell_markers = c("NKG7", "KLRB1", "TMSB4X", "B2M", "CD7", "HLA-B", "BTG1", "CCL4", "IL32", "DNAJB1", "GZMA", "PTPRC", "IFITM1", "IFITM2", "HLA-A", "HLA-C", "CD69", "CORO1A", "CCL5", "CD52")

Endothelialcell_markers = c("FCN3", "DNASE1L3", "IFITM3", "AKAP12", "CRHBP", "CLEC1B", "JUN", "OIT3", "FOS", "FGF23", "SOCS3", "ZFP36", "JUNB", "MARCKSL1", "RAMP2", "SLC2A3", "S100A16", "ACP5", "C8orf4", "HSPA5")

HSPC_markers = c("PRSS57", "ACTG1", "HNRNPA1", "ZFAS1", "SPINK2", "RPS4X", "NPM1", "LDHB", "RPS18", "RPL10A", "RPS3", "RPS14", "RPL3", "FXYD5", "EEF1B2", "SOX4", "RPL7", "RP11-620J15.3", "ENO1", "GNB2L1")

Hepatocyte_markers = c("ALB", "APOA2", "MT1G", "SERPINA1", "APOA1", "MT1H", "FABP1", "APOC3", "MT2A", "TTR", "AHSG", "MT1E", "AMBP", "AFP", "RBP4", "APOH", "MT1X", "MT1F", "APOC1", "IGF2")

Megakaryocyte_markers = c("PF4", "TMSB4X", "PPBP", "PLEK", "GP9", "PDLIM1", "LIMS1", "FERMT3", "ACTG1", "CMTM5", "ACTB", "RAP1B", "ITGA2B", "SDPR", "TAGLN2", "LAT", "ILK", "TPM4", "RSU1", "PFN1")

Macrophage1_markers = c("C1QC", "C1QB", "C1QA",  "CST3", "FTL", "CD74", "MS4A7", "LGMN", "CTSB", "TYROBP", "SEPP1", "SAT1", "LIPA", "PSAP", "HMOX1", "FCGRT", "HLA-DRA", "SLC40A1", "CD68", "HLA-DRB1")

Macrophage2_markers = c("CD74", "HLA-DRA", "CST3", "HLA-DPA1", "HLA-DRB1", "HLA-DPB1", "LGALS1", "TMSB4X", "TMSB10", "HLA-DQB1", "S100A11", "TYROBP", "LYZ", "HLA-DQA1", "VIM", "AIF1", "SRGN", "LSP1", "CLEC10A", "SAT1")

StromalCell_markers = c("SPARC", "COL3A1", "IGFBP7", "COL1A1", "DCN", "IGF2", "IGFBP3", "VIM", "BGN", "CALD1", "COL1A2", "COLEC11", "PTN", "DLK1", "IFITM3", "MDK", "RBP1", "LGALS1", "ASPN", "MEST")

ErythroidCell_markers = c("HBG2", "AHSP", "PRDX2", "HBA2", "HBA1", "HBG1", "HBM", "HBB", "BLVRB", "GYPA", "ALAS2", "SLC25A37", "GYPB", "MYL4", "HMBS", "UROD", "S100A6", "GYPC", "MGST3", "TUBA1B")

ErythroidProgenitorCell_markers = c("FAM178B", "HMGA1", "APOC1", "MYC", "NPM1", "RPL22L1", "LDHB", "EEF1B2", "HNRNPA1", "HMGB1", "PRSS57", "S100A6", "RPS2", "NCL", "H2AFZ", "RP11-620J15.3", "RPL7A", "HIST1H4C", "RPS3", "SNRPD1")

PAGE_matrix_1 = makeSignMatrixPAGE(sign_names = c('B cell',
                                                  'T cell',
                                                  'Endothelial cell',
                                                  'Hematopoietic stem and progenitor cell',
                                                  'Hepatocyte',
                                                  'Megakaryocyte',
                                                  'C1Q+ macrophage',
                                                  'Other macrophage',
                                                  'Stromal cell',
                                                  'Erythroid cell',
                                                  'Erythroid progenitor cell'),
                                   sign_list = list(Bcell_markers,
                                                    Tcell_markers,
                                                    Endothelialcell_markers,
                                                    HSPC_markers,
                                                    Hepatocyte_markers,
                                                    Megakaryocyte_markers,
                                                    Macrophage1_markers,
                                                    Macrophage2_markers,
                                                    StromalCell_markers,
                                                    ErythroidCell_markers,
                                                    ErythroidProgenitorCell_markers))

visium_human12weeksFetalLiver = runPAGEEnrich(gobject = visium_human12weeksFetalLiver, sign_matrix = PAGE_matrix_1)

cell_types_PAGE = colnames(PAGE_matrix_1)
plotMetaDataCellsHeatmap(gobject = visium_human12weeksFetalLiver,
                         metadata_cols = 'leiden_clus',
                         value_cols = cell_types_PAGE,
                         spat_enr_names = 'PAGE',
                         x_text_size = 10,
                         y_text_size = 10)


spatCellPlot2D(gobject = visium_human12weeksFetalLiver,
               spat_enr_names = 'PAGE',
               cell_annotation_values = cell_types_PAGE[2],
               cow_n_col = 2,coord_fix_ratio = 1, point_size = 5, show_legend = T, show_plot = T)


spatDimCellPlot2D(gobject = visium_human12weeksFetalLiver,
                  spat_enr_names = 'PAGE',
                  cell_annotation_values = cell_types_PAGE[1:4],
                  cow_n_col = 1, spat_point_size = 1.5,
                  plot_alignment = 'horizontal',
                  save_param = list(base_width=7, base_height=10))


Spatialcorrelation <- runSpatialEnrich(
  visium_human12weeksFetalLiver,
  sign_matrix = PAGE_matrix_1,
  spat_unit = NULL,
  feat_type = NULL,
  enrich_method = "PAGE",
  output_enrichment = c("original", "zscore")
)


library(xlsx)


write.xlsx(Spatialcorrelation@spatial_enrichment[["cell"]][["rna"]][["PAGE"]]@enrichDT, file = "human12weeksFLenrichment.xlsx",
           sheetName = "Spatialcorrelation@enrichDT", append = FALSE)
