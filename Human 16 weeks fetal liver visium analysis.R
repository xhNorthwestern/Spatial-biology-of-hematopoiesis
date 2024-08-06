results_folder = '~/Desktop/Visium/The third twelve human fetus samples/Matrix/The first eight samples/NMA21-00016_A4'

instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  python_path = '/Users/xuhan/Library/r-miniconda-arm64/envs/Giotto/bin/python')

data_path = '~/Desktop/Visium/The third twelve human fetus samples/Matrix/The first eight samples/NMA21-00016_A4'

visium_human16weeksFetalLiver = createGiottoVisiumObject(visium_dir = data_path,
                                                         expr_data = 'raw',
                                                         png_name = 'tissue_lowres_image.png',
                                                         gene_column_index = 2,
                                                         instructions = instrs)

showGiottoImageNames(visium_human16weeksFetalLiver)

pDataDT(visium_human16weeksFetalLiver)

spatPlot2D(gobject = visium_human16weeksFetalLiver, cell_color = 'in_tissue', point_size = 2,
           cell_color_code = c('0' = 'lightgrey', '1' = 'blue'),
           show_image = T, image_name = 'image')

metadata = pDataDT(visium_human16weeksFetalLiver)
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
visium_human16weeksFetalLiver = subsetGiotto(visium_human16weeksFetalLiver, cell_ids = in_tissue_barcodes)


visium_human16weeksFetalLiver <- filterGiotto(gobject = visium_human16weeksFetalLiver,
                                              expression_threshold = 1,
                                              feat_det_in_min_cells = 50,
                                              min_det_feats_per_cell = 1000,
                                              expression_values = c('raw'),
                                              verbose = T)

visium_human16weeksFetalLiver <- normalizeGiotto(gobject = visium_human16weeksFetalLiver, scalefactor = 6000, verbose = T)


visium_human16weeksFetalLiver <- addStatistics(gobject = visium_human16weeksFetalLiver)


spatPlot2D(gobject = visium_human16weeksFetalLiver, show_image = T, point_alpha = 0.7,
           cell_color = 'nr_feats', color_as_factor = F)


visium_human16weeksFetalLiver <- calculateHVF(gobject = visium_human16weeksFetalLiver, save_plot = TRUE)

gene_metadata = fDataDT(visium_human16weeksFetalLiver)
featgenes = gene_metadata[hvf == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$feat_ID


visium_human16weeksFetalLiver <- runPCA(gobject = visium_human16weeksFetalLiver,
                                        genes_to_use = featgenes)

screePlot(visium_human16weeksFetalLiver, ncp = 30)

plotPCA(gobject = visium_human16weeksFetalLiver)


visium_human16weeksFetalLiver <- runUMAP(visium_human16weeksFetalLiver, dimensions_to_use = 1:10)
plotUMAP(gobject = visium_human16weeksFetalLiver)

visium_human16weeksFetalLiver <- runtSNE(visium_human16weeksFetalLiver, dimensions_to_use = 1:10)
plotTSNE(gobject = visium_human16weeksFetalLiver)


visium_human16weeksFetalLiver <- createNearestNetwork(gobject = visium_human16weeksFetalLiver, dimensions_to_use = 1:10, k = 15)


visium_human16weeksFetalLiver <- doLeidenCluster(gobject = visium_human16weeksFetalLiver, resolution = 0.6, n_iterations = 1000)

plotUMAP(gobject = visium_human16weeksFetalLiver,
         cell_color = 'leiden_clus', show_NN_network = F, point_size = 2.5)


spatDimPlot(gobject = visium_human16weeksFetalLiver, cell_color = 'leiden_clus',
            dim_point_size = 2, spat_point_size = 2.5)

spatDimPlot(gobject = visium_human16weeksFetalLiver, cell_color = 'nr_feats', color_as_factor = F,
            dim_point_size = 2, spat_point_size = 2.5)


spatPlot2D(visium_human16weeksFetalLiver, cell_color = 'leiden_clus',
           coord_fix_ratio = 1)

spatPlot2D(visium_human16weeksFetalLiver, cell_color = 'leiden_clus',
           group_by = 'leiden_clus', coord_fix_ratio = 1,
           cow_n_col = 4, show_legend = F,
           save_param = list(base_width = 14, base_height = 14))


spatPlot2D(visium_human16weeksFetalLiver, cell_color = 'leiden_clus',
           select_cell_groups = '2', coord_fix_ratio = 1, show_other_cells = TRUE,
           cell_color_code = c('2' = 'red'), other_cell_color = "grey", other_point_size = 1.5,
           save_param = list(base_width = 7, base_height = 7))



DG_subset = subsetGiottoLocs(visium_human16weeksFetalLiver,
                             x_max = 14000, x_min = 8000,
                             y_max = -4000, y_min = -8500,
                             return_gobject = TRUE)

spatDimPlot(gobject = DG_subset,
            cell_color = 'leiden_clus', spat_point_size = 5)


gini_markers_subclusters = findMarkers_one_vs_all(gobject = visium_human16weeksFetalLiver,
                                                  method = 'gini',
                                                  expression_values = 'normalized',
                                                  cluster_column = 'leiden_clus',
                                                  min_feats = 20,
                                                  min_expr_gini_score = 0.5,
                                                  min_det_gini_score = 0.5)
topgenes_gini = gini_markers_subclusters[, head(.SD, 2), by = 'cluster']$feats


violinPlot(visium_human16weeksFetalLiver, feats = unique(topgenes_gini), cluster_column = 'leiden_clus',
           strip_text = 8, strip_position = 'right',
           save_param = list(base_width = 5, base_height = 10))


plotMetaDataHeatmap(visium_human16weeksFetalLiver, selected_feats = unique(topgenes_gini),
                    metadata_cols = c('leiden_clus'),
                    x_text_size = 10, y_text_size = 10)


dimFeatPlot2D(visium_human16weeksFetalLiver, expression_values = 'scaled',
              feats = gini_markers_subclusters[, head(.SD, 1), by = 'cluster']$feats,
              cow_n_col = 3, point_size = 0.75,
              save_param = list(base_width = 8, base_height = 8))


scran_markers_subclusters = findMarkers_one_vs_all(gobject = visium_human16weeksFetalLiver,
                                                   method = 'scran',
                                                   expression_values = 'normalized',
                                                   cluster_column = 'leiden_clus')
topgenes_scran = scran_markers_subclusters[, head(.SD, 2), by = 'cluster']$feats


violinPlot(visium_human16weeksFetalLiver, feats = unique(topgenes_scran), cluster_column = 'leiden_clus',
           strip_text = 10, strip_position = 'right',
           save_param = list(base_width = 5))


plotMetaDataHeatmap(visium_human16weeksFetalLiver, selected_feats = topgenes_scran,
                    metadata_cols = c('leiden_clus'))


dimFeatPlot2D(visium_human16weeksFetalLiver, expression_values = 'scaled',
              feats = scran_markers_subclusters[, head(.SD, 1), by = 'cluster']$feats,
              cow_n_col = 3, point_size = 2,
              save_param = list(base_width = 8, base_height = 8))

ErythroidCell1_markers = c("TCL1A", "POU5F1", "PRDX1", "CFLAR", "NDUFAB1", "DPPA5", "LAMTOR5", "ATP6AP2", "PDPN", "SLIRP", "SYCP3", "HENMT1", "DNAJC15", "UQCRH", "APOC1", "SAT1", "METTL5", "MRPL14", "PMAIP1", "GTSF1")

ErythroidProgenitorCell1_markers = c("MT1L", "FAM210B", "PRKAR2B", "VAT1L", "KITLG", "GJA1", "RELN", "TOX3", "CNTNAP4", "GATM", "CRYM", "BEX1", "ENC1", "RTN4", "DDX17", "SERINC5", "MT3", "RNASE1", "FREM2", "CAV2")

ErythroidProgenitorCell2_markers = c("NR4A2", "ARL5B", "BTG2", "KLF10", "NR4A1", "CXCL2", "BHLHE40", "NR4A3", "EGR3", "DIRAS3", "EGR1", "CSRNP1", "SLC2A3", "DUSP2", "THBS1", "MAT2A", "NAMPT", "ZFP36", "FOSL2", "VEGFA")

ErythroidCell2_markers = c("GTSF1", "DPPA5", "TCL1A", "APOC1", "LUZP4", "DNAJC15", "SCGB3A2", "SYCP3", "POU5F1", "PRDX4", "PDPN", "NDUFAB1", "KHDC3L", "S100A10", "HENMT1", "CENPH", "SLC25A5", "WBSCR22", "SNRPC", "LAMTOR5")

ConventionalDendriticCell_markers = c("CYSTM1", "SPRR2F", "IFITM3", "IFITM2", "PODXL", "KRT19", "SPARC", "ANXA5", "KIAA0101", "MDK", "EMX2", "MGST3", "TPD52L1", "PRDX2", "ARL4A", "MT1E", "FABP5", "ID2", "KRT18", "PPDPF")

Neutrophil1_markers = c("CCNB1", "KPNA2", "CKS2", "CKS1B", "UBE2C", "CCNB2", "PTTG1", "HMGB3", "HMGB2", "HIST1H1C", "NUSAP1", "NDUFAB1", "TUBA1C", "DNAJC15", "CALM2", "HN1", "PMAIP1", "RAN", "CDC20", "MND1")

MastProgenitorCell_markers = c("NPY", "PCP4", "TAC1",  "PRSS35", "CPE", "NPNT", "ATP1A2", "GSTA1", "IFI6", "RASL11B", "MGP", "RP5-1050E16.2", "PLAC1", "CTSV", "C3orf58", "NR2F2", "RBP1", "BST2", "JUN", "TSPAN7")

HSPC_markers = c("DCN", "COL1A1", "C7", "DLK1", "COL1A2", "GPC3", "COL6A3", "POSTN", "LGALS1", "OLFML3", "PDGFRA", "COL3A1", "MFAP4", "ISLR", "SEPP1", "ADAMTS1", "DNAJB1", "RGS2", "TIMP1", "IGFBP5")

EffectorTCell_markers = c("TOP2A", "UBE2C", "PTTG1",  "CDK1", "NUSAP1", "HMGB2", "TPX2", "ARL6IP1", "KRT19", "CDKN3", "CENPF", "CDC20", "CCNA2", "PRC1", "KIAA0101", "SMC4", "CCNB2", "NUF2", "FBXO5", "LGALS1")

BCell_markers = c("PAGE2B", "PAGE5", "PAGE2",  "MAGEA4", "CENPH", "WBSCR22", "PRAME", "SYCP3", "SOD1", "PRDX4", "HMGB3", "GTSF1", "PSMD10", "TEX30", "MAEL", "ERP29", "CCT2", "HENMT1", "HPRT1", "SMS")

Neutrophil2_markers = c("CCL21", "A2M", "TM4SF1", "FN1", "ESAM", "GNG11", "PRCP", "ANXA2", "IGFBP7", "CDH5", "RAMP2", "CD74", "CALCRL", "CD34", "EMCN", "CCL2", "C8orf4", "TGFBR2", "SPP1", "PCAT19")

Neutrophil3_markers = c("CCL4", "CCL3", "C1QB",  "CD74", "HLA-DRA", "CXCL10", "CCL4L2", "AIF1", "SRGN", "C1QA", "HLA-DPA1", "B2M", "TYROBP", "IFNG", "FCER1G", "CCL3L3", "RGS1", "HLA-B", "C1QC", "HLA-DPB1")


PAGE_matrix_1 = makeSignMatrixPAGE(sign_names = c('HBM high_Erythroid Cell',
                                                  'REXO2 high_Erythroid Progenitor Cell',
                                                  'NPM1 high_Erythroid Progenitor Cell',
                                                  'HBB high_Erythroid Cell',
                                                  'Conventional Dendritic Cell',
                                                  'ELANE high_Neutrophil',
                                                  'Mast Progenitor Cell',
                                                  'HSPC',
                                                  'Effector T Cell',
                                                  'B Cell',
                                                  'PRTN3 high_Neutrophil',
                                                  'C1Q+ Macrophage'),
                                   sign_list = list(ErythroidCell1_markers,
                                                    ErythroidProgenitorCell1_markers,
                                                    ErythroidProgenitorCell2_markers,
                                                    ErythroidCell2_markers,
                                                    ConventionalDendriticCell_markers,
                                                    Neutrophil1_markers,
                                                    MastProgenitorCell_markers,
                                                    HSPC_markers,
                                                    EffectorTCell_markers,
                                                    BCell_markers,
                                                    Neutrophil2_markers,
                                                    Neutrophil3_markers))

visium_human16weeksFetalLiver = runPAGEEnrich(gobject = visium_human16weeksFetalLiver, sign_matrix = PAGE_matrix_1)

cell_types_PAGE = colnames(PAGE_matrix_1)
plotMetaDataCellsHeatmap(gobject = visium_human16weeksFetalLiver,
                         metadata_cols = 'leiden_clus',
                         value_cols = cell_types_PAGE,
                         spat_enr_names = 'PAGE',
                         x_text_size = 8,
                         y_text_size = 8)


spatCellPlot2D(gobject = visium_human16weeksFetalLiver,
               spat_enr_names = 'PAGE',
               cell_annotation_values = cell_types_PAGE[1:4],
               cow_n_col = 2,coord_fix_ratio = 1, point_size = 1.25, show_legend = T)


spatDimCellPlot2D(gobject = visium_human16weeksFetalLiver,
                  spat_enr_names = 'PAGE',
                  cell_annotation_values = cell_types_PAGE[1:4],
                  cow_n_col = 1, spat_point_size = 1,
                  plot_alignment = 'horizontal',
                  save_param = list(base_width=7, base_height=10))


Spatialcorrelation <- getSpatialEnrichment(
  visium_human16weeksFetalLiver,
  spat_unit = NULL,
  feat_type = NULL,
  enrichm_name = "PAGE",
  output = c("spatEnrObj", "data.table"),
  copy_obj = TRUE,
  set_defaults = TRUE
)

library(xlsx)


write.xlsx(Spatialcorrelation@enrichDT, file = "myworkbookhumanfetalliver16weeks.xlsx",
           sheetName = "Spatialcorrelation@enrichDT", append = FALSE)


