results_folder = '~/Desktop/Visium/The fourth sixteen human bone marrow samples/NMB22-955_D1/2'

library(Giotto)

instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  python_path = '/Users/xuhan/Library/r-miniconda-arm64/envs/giotto_env/bin/pythonw')

data_path = '~/Desktop/Visium/The fourth sixteen human bone marrow samples/NMB22-955_D1'

visium_955D1 = createGiottoVisiumObject(visium_dir = data_path,
                                        expr_data = 'raw',
                                        png_name = 'tissue_lowres_image.png',
                                        gene_column_index = 2,
                                        instructions = instrs)

showGiottoImageNames(visium_955D1)

pDataDT(visium_955D1)

spatPlot2D(gobject = visium_955D1, cell_color = 'in_tissue', point_size = 2,
           cell_color_code = c('0' = 'lightgrey', '1' = 'blue'),
           show_image = T, image_name = 'image')

metadata = pDataDT(visium_955D1)
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
visium_955D1 = subsetGiotto(visium_955D1, cell_ids = in_tissue_barcodes)


visium_955D1 <- filterGiotto(gobject = visium_955D1,
                             min_det_feats_per_cell = 50,
                             expression_values = c('raw'),
                             verbose = T)

visium_955D1 <- normalizeGiotto(gobject = visium_955D1, scalefactor = 6000, verbose = T)


visium_955D1 <- addStatistics(gobject = visium_955D1)


spatPlot2D(gobject = visium_955D1, show_image = T, point_alpha = 0.7,
           cell_color = 'nr_feats', color_as_factor = F)


visium_955D1 <- calculateHVF(gobject = visium_955D1, save_plot = TRUE)

gene_metadata = fDataDT(visium_955D1)
featgenes = gene_metadata[hvf == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$feat_ID


visium_955D1 <- runPCA(gobject = visium_955D1,
                       feats_to_use = featgenes)

screePlot(visium_955D1, ncp = 30)

plotPCA(gobject = visium_955D1)


visium_955D1 <- runUMAP(visium_955D1, dimensions_to_use = 1:10)
plotUMAP(gobject = visium_955D1)


visium_955D1 <- createNearestNetwork(gobject = visium_955D1, dimensions_to_use = 1:10, k = 15)


visium_955D1 <- doLeidenCluster(gobject = visium_955D1, resolution = 0.05, n_iterations = 1000)

dimPlot2D(gobject = visium_955D1,     dim_reduction_name = 'umap',
          cell_color = "leiden_clus", show_NN_network = F, label_size = 10, point_size = 6.0, legend_text = 20,
          legend_symbol_size = 6, axis_text = 16, axis_title = 16,
          save_param = list(save_name = "7_Annotation"))

dimPlot2D(gobject = visium_955D1,     dim_reduction_name = 'umap',
          cell_color = "leiden_clus", show_NN_network = F, label_size = 10, point_size = 6.0, show_legend = F, legend_text = 20,
          legend_symbol_size = 6, axis_text = 16, axis_title = 16,
          save_param = list(save_name = "7_a_Annotation"))


spatDimPlot(gobject = visium_955D1, cell_color = 'leiden_clus',
            dim_point_size = 2, spat_point_size = 2.5)

spatDimPlot(gobject = visium_955D1, cell_color = 'nr_feats', color_as_factor = F,
            dim_point_size = 2, spat_point_size = 2.5)

spatPlot2D(visium_955D1, cell_color = 'leiden_clus', label_size = 10,
           coord_fix_ratio = 1, point_size = 4.0, show_legend = F, axis_text = 16,
           axis_title = 16)

spatPlot2D(visium_955D1, cell_color = 'leiden_clus',
           group_by = 'leiden_clus', coord_fix_ratio = 1,
           cow_n_col = 4, show_legend = F,
           save_param = list(base_width = 14, base_height = 14))


gini_markers_subclusters = findMarkers_one_vs_all(gobject = visium_955D1,
                                                  method = 'gini',
                                                  expression_values = 'normalized',
                                                  cluster_column = 'leiden_clus',
                                                  min_feats = 20,
                                                  min_expr_gini_score = 0.5,
                                                  min_det_gini_score = 0.5)
topgenes_gini = gini_markers_subclusters[, head(.SD, 2), by = 'cluster']$feats


violinPlot(visium_955D1, feats = unique(topgenes_gini), cluster_column = 'leiden_clus',
           strip_text = 8, strip_position = 'right',
           save_param = list(base_width = 5, base_height = 10))


plotMetaDataHeatmap(visium_955D1, selected_feats = unique(topgenes_gini),
                    metadata_cols = c('leiden_clus'),
                    x_text_size = 10, y_text_size = 10)


dimFeatPlot2D(visium_955D1, expression_values = 'scaled',
              feats = gini_markers_subclusters[, head(.SD, 1), by = 'cluster']$feats,
              cow_n_col = 3, point_size = 0.75,
              save_param = list(base_width = 8, base_height = 8))


scran_markers_subclusters = findMarkers_one_vs_all(gobject = visium_955D1,
                                                   method = 'scran',
                                                   expression_values = 'normalized',
                                                   cluster_column = 'leiden_clus')
topgenes_scran = scran_markers_subclusters[, head(.SD, 2), by = 'cluster']$feats


violinPlot(visium_955D1, feats = unique(topgenes_scran), cluster_column = 'leiden_clus',
           strip_text = 10, strip_position = 'right',
           save_param = list(base_width = 5))

plotMetaDataHeatmap(visium_955D1,  x_text_size = 20,
                    x_text_angle = 45,
                    y_text_size = 20,
                    strip_text_size = 20, selected_feats = topgenes_scran[c(1:9, 11:14, 16:18)],
                    metadata_cols = c('leiden_clus'))


dimFeatPlot2D(visium_955D1, expression_values = 'scaled',
              feats = scran_markers_subclusters[, head(.SD, 1), by = 'cluster']$feats,
              cow_n_col = 3, point_size = 1,
              save_param = list(base_width = 8, base_height = 8))


Neutrophil1_markers = c("DEFA3", "SLPI", "BPI", "MS4A3", "DEFA4", "CD24", "SRGN", "CEACAM8", "NSUN5", "NKG7", "CD63", "CEACAM6", "HMGN2", "CD59", "RNASE2", "OLFM4", "RNASE3", "DUT", "RETN", "CYBA")

Neutrophil2_markers = c("CTSG", "PRTN3", "ELANE", "CST7", "AZU1", "PRSS57", "MPO", "CLEC11A", "DEFA4", "SRGN", "CALR", "LAIR2", "RNASE3", "P4HB", "SMIM24", "NUCB2", "SERPINB1", "RPS14", "HSPA5", "RNASE2")

Neutrophil3_markers = c("S100A12", "PGLYRP1", "LTF", "LCN2", "CAMP", "CRISP3", "S100A9", "S100A8", "MMP8", "PADI4", "PLBD1", "MMP9", "CD177", "TFF3", "HP", "ORM1", "GCA", "CLC", "IFITM2", "FPR1")

Monocyte_markers = c("LGALS1", "CST3", "ANXA2", "HLA-DRB1", "S100A10", "HLA-DRA", "CD74", "HLA-DPA1", "MS4A6A", "HLA-DPB1", "HLA-DRB6", "FOS", "NFKBIA", "ZFP36", "S100A4", "TMSB10", "AP1S2", "SELL", "GPX1", "TNFSF13B")

Erythroidcell_markers = c("HBA2", "HBA1", "HBD", "HBB", "CA1", "HBM", "GYPB", "AHSP", "SLC4A1", "PRDX2", "CA2", "BLVRB", "GYPA", "ALAS2", "HMBS", "SLC25A37", "EPB42", "TFRC", "UROD", "HEMGN")

Erythroidprogenitorcell_markers = c("APOC1", "FAM178B", "CNRIP1", "REXO2", "BLVRB", "PVT1", "PRDX2", "AKR1C3", "FTH1", "ATPIF1", "AHSP", "ALDH1A1", "NPM1", "HEBP1", "HSPE1", "SOD1", "TXNL1", "SYNGR1", "TMEM14C", "KLF1")

Bcell1_markers = c("IGKC", "IGHG1", "IGHG3",  "IGLC2", "IGHG4", "IGHA1", "IGLC3", "IGHA2", "IGHG2", "MZB1", "JCHAIN", "IGHGP", "DERL3", "HERPUD1", "SEC11C", "PRDX4", "SSR4", "XBP1", "TNFRSF17", "ITM2C")

Bcell2_markers = c("IGLL1", "VPREB3", "VPREB1",  "TCL1A", "DNTT", "IGHM", "CD79A", "CD79B", "CMTM8", "SOX4", "CD52", "STMN1", "PRDX1", "HMGB1", "HLA-DRA", "VDAC1", "CD74", "RRM2", "PTMA", "C12orf57")

Tcell_markers = c("CCL5", "GZMA", "IL32",  "TRAC", "TRBC2", "CD3D", "TRBC1", "GZMH", "CD2", "CD52", "MMP9", "IL2RG", "PSMB9", "B2M", "CD53", "DDX5", "BIN2", "NKG7", "HLA-E", "LITAF")

Macrophage_markers = c("SEPP1", "C1QB", "C1QA",  "VCAM1", "CCL3L3", "CCL4", "MS4A7", "C1QC", "CXCL12", "CCL3", "CTSB", "LGMN", "HMOX1", "CD163", "MS4A4A", "HLA-DPA1", "CD5L", "CCL4L2", "RGS1", "RNASE1")

PAGE_matrix_1 = makeSignMatrixPAGE(sign_names = c('Neutrophil_S100A9 high',
                                                  'Neutrophil_PRTN3 high',
                                                  'Neutrophil_LTF high',
                                                  'Monocyte',
                                                  'Erythroid Cell',
                                                  'Erythroid Progenitor Cell',
                                                  'B Cell (Plasmocyte)',
                                                  'B Cell (Centrocyte)',
                                                  'T Cell',
                                                  'C1Q+ Macrophage'),
                                   sign_list = list(Neutrophil1_markers,
                                                    Neutrophil2_markers,
                                                    Neutrophil3_markers,
                                                    Monocyte_markers,
                                                    Erythroidcell_markers,
                                                    Erythroidprogenitorcell_markers,
                                                    Bcell1_markers,
                                                    Bcell2_markers,
                                                    Tcell_markers,
                                                    Macrophage_markers))

visium_955D1 = runPAGEEnrich(gobject = visium_955D1, sign_matrix = PAGE_matrix_1)

cell_types_PAGE = colnames(PAGE_matrix_1)
plotMetaDataCellsHeatmap(gobject = visium_955D1,
                         metadata_cols = 'leiden_clus',
                         value_cols = cell_types_PAGE,
                         spat_enr_names = 'PAGE',
                         x_text_size = 20,
                         y_text_size = 20)


spatCellPlot2D(gobject = visium_955D1,
               spat_enr_names = 'PAGE',
               cell_annotation_values = cell_types_PAGE[1:4],
               cow_n_col = 2,coord_fix_ratio = 1, point_size = 1.25, show_legend = T)


spatDimCellPlot2D(gobject = visium_955D1,
                  spat_enr_names = 'PAGE',
                  cell_annotation_values = cell_types_PAGE[1:4],
                  cow_n_col = 1, spat_point_size = 1,
                  plot_alignment = 'horizontal',
                  save_param = list(base_width=7, base_height=10))


Spatialcorrelation <- get_spatial_enrichment(
  visium_955D1,
  spat_unit = NULL,
  feat_type = NULL,
  enrichm_name = "PAGE",
  output = c("spatEnrObj", "data.table"),
  copy_obj = TRUE,
  set_defaults = TRUE
)

library(writexl)

write_xlsx(Spatialcorrelation@enrichDT, "~/Desktop/Visium/The fourth sixteen human bone marrow samples/NMB22-955_D1/2/955D1 correlation analysis.xlsx")


