library(Giotto)


genv_exists = checkGiottoEnvironment()

# 1. ** SET WORKING DIRECTORY WHERE PROJECT OUPUTS WILL SAVE TO **
results_folder = '~/Desktop/Xenium 5K assay/Mouse/Giotto Analysis/E135 FL/2'

# 3. Create Giotto instructions
# Directly saving plots to the working directory without rendering them in the editor saves time.
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  return_plot = FALSE)

# ** SET PATH TO FOLDER CONTAINING XENIUM DATA **
xenium_folder = '~/Desktop/Xenium 5K assay/Mouse/Giotto Analysis/E135 FL'

# general files (some are supplemental files)
settings_path = paste0(xenium_folder, 'experiment.xenium')
he_img_path = paste0(xenium_folder, 'morphology_focus')
if_img_path = paste0(xenium_folder, 'morphology.ome.tif')

# files (SUBCELLULAR): (tutorial focuses on working with these files)
cell_bound_path = paste0('~/Desktop/Xenium 5K assay/Mouse/Giotto Analysis/E135 FL/cell_boundaries.csv.gz')
nuc_bound_path = paste0('~/Desktop/Xenium 5K assay/Mouse/Giotto Analysis/E135 FL/nucleus_boundaries.csv.gz')
tx_path = paste0('~/Desktop/Xenium 5K assay/Mouse/Giotto Analysis/E135 FL/transcripts.csv.gz')
feat_meta_path = paste0('~/Desktop/Xenium 5K assay/Mouse/Giotto Analysis/E135 FL/cell_feature_matrix/features.tsv.gz') # (also used in aggregate)

# files (AGGREGATE):
expr_mat_path = paste0(xenium_folder, 'cell_feature_matrix')
cell_meta_path = paste0(xenium_folder, 'cells.csv.gz') # contains spatlocs


# load features metadata
# (make sure cell_feature_matrix folder is unpacked)
feature_dt = data.table::fread(feat_meta_path, header = FALSE)
colnames(feature_dt) = c('ensembl_ID','feat_name','feat_type')

# find the feature IDs that belong to each feature type
feature_dt[, table(feat_type)]
feat_types = names(feature_dt[, table(feat_type)])

feat_types_IDs = lapply(
  feat_types, function(type) feature_dt[feat_type == type, unique(feat_name)]
)
names(feat_types_IDs) = feat_types

tx_dt = data.table::fread(tx_path)
data.table::setnames(x = tx_dt,
                     old = c('feature_name', 'x_location', 'y_location'),
                     new = c('feat_ID', 'x', 'y'))
cat('Transcripts info available:\n ', paste0('"', colnames(tx_dt), '"'), '\n',
    'with', tx_dt[,.N], 'unfiltered detections\n')

# filter by qv (Phred score)
tx_dt_filtered = tx_dt[qv >= 20]
cat('and', tx_dt_filtered[,.N], 'filtered detections\n\n')

# separate detections by feature type
tx_dt_types = lapply(
  feat_types_IDs, function(types) tx_dt_filtered[feat_ID %in% types]
)

invisible(lapply(seq_along(tx_dt_types), function(x) {
  cat(names(tx_dt_types)[[x]], 'detections: ', tx_dt_types[[x]][,.N], '\n')
}))

gpoints_list = lapply(
  tx_dt_types, function(x) createGiottoPoints(x = x)
) # 208.499 sec elapsed


library(raster)


# preview QC probe detections
plot(gpoints_list$`Unassigned Codeword`,
     point_size = 0.3,
     main = 'Unassigned Codeword')
plot(gpoints_list$`Negative Control Codeword`,
     point_size = 0.3,
     main = 'Negative Control Codeword')
plot(gpoints_list$`Negative Control Probe`,
     point_size = 0.3,
     main = 'Negative Control Probe')

# preview two genes (slower)
plot(gpoints_list$`Gene Expression`,  # 77.843 sec elapsed
     feats = c('Alas2', 'C1qb'))
tx_dt_types$`Gene Expression`[feat_ID %in% c('Alas2', 'C1qb'), table(feat_ID)]

cellPoly_dt = data.table::fread(cell_bound_path)
nucPoly_dt = data.table::fread(nuc_bound_path)

data.table::setnames(cellPoly_dt,
                     old = c('cell_id', 'vertex_x', 'vertex_y'),
                     new = c('poly_ID', 'x', 'y'))
data.table::setnames(nucPoly_dt,
                     old = c('cell_id', 'vertex_x', 'vertex_y'),
                     new = c('poly_ID', 'x', 'y'))

gpoly_cells = createGiottoPolygonsFromDfr(segmdfr = cellPoly_dt,
                                          name = 'cell',
                                          calc_centroids = TRUE)
gpoly_nucs = createGiottoPolygonsFromDfr(segmdfr = nucPoly_dt,
                                         name = 'nucleus',
                                         calc_centroids = TRUE)


plot(x = gpoly_nucs, point_size = 0.1, type = 'centroid')

xenium_gobj = createGiottoObjectSubcellular(
  gpoints = list(rna = gpoints_list$`Gene Expression`,
                 blank_code = gpoints_list$`Unassigned Codeword`,
                 neg_code = gpoints_list$`Negative Control Codeword`,
                 neg_probe = gpoints_list$`Negative Control Probe`),
  gpolygons = list(cell = gpoly_cells,
                   nucleus = gpoly_nucs),
  instructions = instrs
)

showGiottoSpatialInfo(xenium_gobj)

showGiottoSpatLocs(xenium_gobj)

spatPlot2D(xenium_gobj,
           spat_unit = 'cell',
           point_shape = 'no_border',
           point_size = 0.5,
           point_alpha = 0.4,
           save_param = list(
             base_width = 7,
             base_height = 7,
             save_name = '1_spatplot'))

xenium_gobj = calculateOverlapRaster(xenium_gobj,
                                     spatial_info = 'cell',
                                     feat_info = 'rna')

showGiottoSpatialInfo(xenium_gobj)

xenium_gobj = overlapToMatrix(xenium_gobj,
                              poly_info = 'cell',
                              feat_info = 'rna',
                              name = 'raw')

showGiottoExpression(xenium_gobj)

# Print a preview of all available features metadata
showGiottoFeatMetadata(xenium_gobj)

xenium_gobj = filterGiotto(gobject = xenium_gobj,
                           spat_unit = 'cell',
                           poly_info = 'cell',
                           expression_threshold = 1,
                           feat_det_in_min_cells = 3,
                           min_det_feats_per_cell = 5)


xenium_gobj = addStatistics(xenium_gobj, expression_values = 'raw')

showGiottoCellMetadata(xenium_gobj)
showGiottoFeatMetadata(xenium_gobj)

xenium_gobj = normalizeGiotto(gobject = xenium_gobj,
                              spat_unit = 'cell',
                              scalefactor = 5000,
                              verbose = T)

xenium_gobj = calculateHVF(gobject = xenium_gobj,
                           spat_unit = 'cell',
                           save_param = list(
                             save_name = '2_HVF'))

cat(fDataDT(xenium_gobj)[, sum(hvf == 'yes')], 'hvf found')

xenium_gobj = runPCA(gobject = xenium_gobj,
                     spat_unit = 'cell',
                     expression_values = 'scaled',
                     feats_to_use = NULL,
                     scale_unit = F,
                     center = F)

# Visualize Screeplot and PCA
screePlot(xenium_gobj,
          ncp = 20,
          save_param = list(
            save_name = '3a_screePlot'))
showGiottoDimRed(xenium_gobj)
plotPCA(xenium_gobj,
        spat_unit = 'cell',
        dim_reduction_name = 'pca',
        dim1_to_use = 1,
        dim2_to_use = 2)

xenium_gobj = runtSNE(xenium_gobj,
                      dimensions_to_use = 1:10,
                      spat_unit = 'cell')
xenium_gobj = runUMAP(xenium_gobj,
                      dimensions_to_use = 1:10,
                      spat_unit = 'cell')

plotTSNE(xenium_gobj,
         point_size = 0.01,
         save_param = list(
           save_name = '4a_tSNE'))
plotUMAP(xenium_gobj,
         point_size = 0.01,
         save_param = list(
           save_name = '4b_UMAP'))

xenium_gobj = createNearestNetwork(xenium_gobj,
                                   dimensions_to_use = 1:10,
                                   k = 10,
                                   spat_unit = 'cell')

xenium_gobj = doLeidenCluster(xenium_gobj,
                              resolution = 0.6,
                              n_iterations = 100,
                              spat_unit = 'cell')

# visualize UMAP cluster results
plotUMAP(gobject = xenium_gobj,
         spat_unit = 'cell',
         cell_color = 'leiden_clus',
         show_legend = FALSE,
         point_size = 0.01,
         point_shape = 'no_border',
         save_param = list(save_name = '5a_umap_leiden'))

plotTSNE(gobject = xenium_gobj,
         spat_unit = 'cell',
         cell_color = 'leiden_clus',
         show_legend = FALSE,
         point_size = 0.01,
         point_shape = 'no_border',
         save_param = list(save_name = '5_umap_leiden'))

spatPlot2D(gobject = xenium_gobj,
           spat_unit = 'cell',
           cell_color = 'leiden_clus',
           point_size = 0.1,
           point_shape = 'no_border',
           background_color = 'black',
           show_legend = TRUE,
           save_param = list(
             save_name = '6_spat_leiden',
             base_width = 15,
             base_height = 15))

spatInSituPlotPoints(xenium_gobj,
                     show_image = FALSE,
                     feats = NULL,
                     point_size = 0.05,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_alpha = 1,
                     polygon_color = 'black',
                     polygon_line_size = 0.01,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = TRUE,
                     coord_fix_ratio = TRUE,
                     save_para = list(
                       save_name = '7_polys'))

markers = findMarkers_one_vs_all(gobject = xenium_gobj,
                                 method = 'gini',
                                 expression_values = 'normalized',
                                 cluster_column = 'leiden_clus',
                                 min_feats = 1, rank_score = 2)
# Display details about the marker genes in-console
markers[, head(.SD, 2), by = 'cluster']

# violinplot
topgini_genes = unique(markers[, head(.SD, 3), by = 'cluster']$feats)
violinPlot(xenium_gobj, feats = topgini_genes, cluster_column = 'leiden_clus', strip_position = 'right')

plotMetaDataHeatmap(xenium_gobj, expression_values = 'scaled',
                    metadata_cols = c('leiden_clus'),
                    selected_feats = topgini_genes)

markers_scran = findMarkers_one_vs_all(gobject=xenium_gobj, method="scran",
                                       expression_values="normalized", cluster_column='leiden_clus', min_feats=3)
markergenes_scran = unique(markers_scran[, head(.SD, 3), by="cluster"][["feats"]])

plotMetaDataHeatmap(xenium_gobj, expression_values = "normalized", metadata_cols = 'leiden_clus',
                    selected_feats = markergenes_scran,
                    y_text_size = 8, show_values = 'zscores_rescaled',
                    save_param = list(save_name = '9_a_metaheatmap'))

topgenes_scran = markers_scran[, head(.SD, 1), by = 'cluster']$feats
# violinplot
violinPlot(xenium_gobj, feats = unique(topgenes_scran), cluster_column = 'leiden_clus',
           strip_text = 10, strip_position = 'right',
           save_param = list(save_name = '9_b_violinplot_scran', base_width = 5))

Xenium_labels<-c("Stromal cell",#1
                 "Erythroid progenitor cell",#2
                 "Hepatocyte",#3
                 "C1q+ macrophage",#4
                 "Hepatocyte",#5
                 "Erythroid cell",#6
                 "Erythroid cell",#7
                 "Erythroid cell",#8
                 "Erythroid cell",#9
                 "Endothelial cell",#10
                 "Hepatocyte",#11
                 "Stromal cell",#12
                 "Erythroid cell",#13
                 "Endothelial cell",#14
                 "Erythroid progenitor cell",#15
                 "Erythroid progenitor cell",#16
                 "C1q+ macrophage",#17
                 "Hepatocyte",#18
                 "Hepatocyte",#19
                 "Neutrophil",#20
                 "Erythroid progenitor cell",#21
                 "Erythroid cell",#22
                 "Stromal cell",#23
                 "Stromal cell",#24
                 "Megakaryocyte",#25
                 "Erythroid cell",#26
                 "Erythroid cell",#27
                 "Hepatocyte",#28
                 "Erythroid cell",#29
                 "Stromal cell",#30
                 "Hepatocyte",#31
                 "Erythroid cell")#32
names(Xenium_labels)<-1:32
xenium_gobj<-annotateGiotto(gobject = xenium_gobj, annotation_vector = Xenium_labels ,
                            cluster_column = 'leiden_clus', name = 'Xenium_labels')
dimPlot2D(gobject = xenium_gobj,     dim_reduction_name = 'umap',
          cell_color = "Xenium_labels", show_NN_network = F, point_size = 1.5,
          save_param = list(save_name = "10_Annotation"))

plotMetaDataHeatmap(xenium_gobj, expression_values = 'scaled',
                    metadata_cols = c('Xenium_labels'),
                    x_text_size = 14,
                    x_text_angle = 45,
                    y_text_size = 14,
                    strip_text_size = 14,
                    selected_feats = topgini_genes)

plotMetaDataHeatmap(xenium_gobj, expression_values = "scaled", metadata_cols = c('Xenium_labels'),
                    selected_feats = markergenes_scran,
                    x_text_size = 14,
                    x_text_angle = 45,
                    y_text_size = 14,
                    strip_text_size = 14, show_values = 'zscores_rescaled',
                    save_param = list(save_name = '10_a_metaheatmap'))

mycolorcode = c('plum', 'red', 'magenta', 
                'purple', 'yellowgreen', 'yellow', 'gray', 'orange')
names(mycolorcode) = c('Neutrophil', 'Erythroid cell', 'Erythroid progenitor cell', 
                       'Stromal cell', 'C1q+ macrophage', 'Megakaryocyte', 'Hepatocyte', 'Endothelial cell')
dimPlot2D(gobject = xenium_gobj,     dim_reduction_name = 'umap',
          cell_color = "Xenium_labels", cell_color_code = mycolorcode, show_NN_network = F, label_size = 6, point_size = 1.5, legend_text = 16,
          legend_symbol_size = 6, axis_text = 16, axis_title = 16,
          save_param = list(save_name = "10_a_Annotation"))

dimPlot2D(gobject = xenium_gobj,     dim_reduction_name = 'umap',
          cell_color = "Xenium_labels", cell_color_code = mycolorcode, show_NN_network = F, label_size = 6, point_size = 1.5, show_legend = F, legend_text = 16,
          legend_symbol_size = 6, axis_text = 16, axis_title = 16,
          save_param = list(save_name = "10_d_Annotation"))

dimPlot2D(gobject = xenium_gobj,     dim_reduction_name = 'umap',
          cell_color = "Xenium_labels", cell_color_code = mycolorcode, show_NN_network = F, show_center_label = F, label_size = 6, point_size = 1.5, show_legend = F, legend_text = 16,
          legend_symbol_size = 6, axis_text = 16, axis_title = 16,
          save_param = list(save_name = "10_e_Annotation"))

spatInSituPlotPoints(xenium_gobj,
                     show_image = FALSE,
                     feats = NULL,
                     point_size = 0.05,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_alpha = 1,
                     polygon_color = 'black',
                     polygon_line_size = 0.01,
                     polygon_fill = 'Xenium_labels',
                     polygon_fill_code = mycolorcode,
                     polygon_fill_as_factor = TRUE,
                     coord_fix_ratio = TRUE,
                     save_para = list(
                       save_name = '8_polys'))

plotStatDelaunayNetwork(gobject = xenium_gobj, maximum_distance = 400)

xenium_gobj = createSpatialNetwork(gobject = xenium_gobj, minimum_k = 2,
                                   maximum_distance_delaunay = 400)
xenium_gobj = createSpatialNetwork(gobject = xenium_gobj, minimum_k = 2,
                                   method = 'kNN', k = 10)
showGiottoSpatNetworks(xenium_gobj)
# visualize the two different spatial networks
spatPlot(gobject = xenium_gobj, show_network = T,
         network_color = 'blue', spatial_network_name = 'Delaunay_network',
         point_size = 2.5, cell_color = 'Xenium_labels')
spatPlot(gobject = xenium_gobj, show_network = T,
         network_color = 'blue', spatial_network_name = 'kNN_network',
         point_size = 2.5, cell_color = 'Xenium_labels')

set.seed(seed = 2841)
cell_proximities = cellProximityEnrichment(gobject = xenium_gobj,
                                           cluster_column = 'Xenium_labels',
                                           spatial_network_name = 'Delaunay_network',
                                           adjust_method = 'fdr',
                                           number_of_simulations = 1000)
# barplot
cellProximityBarplot(gobject = xenium_gobj,
                     CPscore = cell_proximities,
                     min_orig_ints = 5, min_sim_ints = 5, p_val = 0.05)
## heatmap
cellProximityHeatmap(gobject = xenium_gobj, CPscore = cell_proximities,
                     order_cell_types = T, scale = T,
                     color_breaks = c(-1.5, 0, 1.5),
                     color_names = c('blue', 'white', 'red'))

# network
cellProximityNetwork(gobject = xenium_gobj, CPscore = cell_proximities,
                     remove_self_edges = T, only_show_enrichment_edges = T)

# network with self-edges
cellProximityNetwork(gobject = xenium_gobj, CPscore = cell_proximities,
                     remove_self_edges = F, self_loop_strength = 0.3,
                     only_show_enrichment_edges = F,
                     rescale_edge_weights = T,
                     node_size = 8,
                     edge_weight_range_depletion = c(1, 2),
                     edge_weight_range_enrichment = c(2,5))


subloc1 = subsetGiottoLocs(xenium_gobj,
                           x_min = 550, x_max = 750,
                           y_min = 350, y_max = 500,
                           poly_info = 'cell')

# show subset of genes

spatInSituPlotPoints(subloc1,
                     show_image = FALSE,
                     feats = list('rna' = c(
                       "Car2", "C1qb", "Kit")),
                     feats_color_code = c(
                       "Car2" = 'green',
                       'C1qb' = 'blue',
                       'Kit' = 'red'),
                     point_size = 0.5,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'white',
                     polygon_line_size = 0.1,
                     polygon_fill = 'Xenium_labels',
                     polygon_fill_code = mycolorcode,
                     polygon_fill_as_factor = TRUE,
                     coord_fix_ratio = TRUE,
                     save_param = list(
                       save_name = '14_subset_in_situ'))

spatInSituPlotPoints(subloc1,
                     show_image = FALSE,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'white',
                     polygon_line_size = 0.1,
                     polygon_fill = 'Xenium_labels',
                     polygon_fill_code = mycolorcode,
                     polygon_fill_as_factor = TRUE,
                     coord_fix_ratio = TRUE,
                     axis_text = 12,
                     axis_title = 12,
                     legend_text = 16,
                     show_legend = TRUE,
                     save_param = list(
                       save_name = '14_b_subset_in_situ'))


spatPlot2D(gobject = subloc1, point_size = 6.0,
           legend_text = 16,
           legend_symbol_size = 4,
           axis_text = 12,
           axis_title = 12,
           cell_color = 'Xenium_labels', cell_color_code = mycolorcode)

spatPlot2D(gobject = subloc1, point_size = 6.0,
           legend_text = 16,
           legend_symbol_size = 4,
           axis_text = 12,
           axis_title = 12,
           cell_color = 'Xenium_labels', cell_color_code = mycolorcode,
           select_cell_groups = c('C1q+ macrophage','Erythroid cell'), show_other_cells = F)


library(writexl)

write_xlsx(xenium_gobj@cell_metadata[["cell"]][["rna"]]@metaDT, "~/Desktop/Xenium 5K assay/Mouse/Giotto Analysis/E135 FL/2/E135FLcellmeta.xlsx")


write_xlsx(cell_proximities[["enrichm_res"]], "~/Desktop/Xenium 5K assay/Mouse/Giotto Analysis/E135 FL/2/E135FLcellproximity.xlsx")

