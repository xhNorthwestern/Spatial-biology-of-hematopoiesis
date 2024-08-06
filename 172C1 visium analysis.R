results_folder = '~/Desktop/Visium/The fourth sixteen human bone marrow samples/NMB21-172_C1'

instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  python_path = '/Users/xuhan/Library/r-miniconda-arm64/envs/Giotto/bin/python')

data_path = '~/Desktop/Visium/The fourth sixteen human bone marrow samples/NMB21-172_C1'

visium_172C1 = createGiottoVisiumObject(visium_dir = data_path,
                                        expr_data = 'raw',
                                        png_name = 'tissue_lowres_image.png',
                                        gene_column_index = 2,
                                        instructions = instrs)

showGiottoImageNames(visium_172C1)

pDataDT(visium_172C1)

spatPlot2D(gobject = visium_172C1, cell_color = 'in_tissue', point_size = 2,
           cell_color_code = c('0' = 'lightgrey', '1' = 'blue'),
           show_image = T, image_name = 'image')

metadata = pDataDT(visium_172C1)
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
visium_172C1 = subsetGiotto(visium_172C1, cell_ids = in_tissue_barcodes)


visium_172C1 <- filterGiotto(gobject = visium_172C1,
                             expression_threshold = 1,
                             feat_det_in_min_cells = 10,
                             min_det_feats_per_cell = 50,
                             expression_values = c('raw'),
                             verbose = T)

visium_172C1 <- normalizeGiotto(gobject = visium_172C1, scalefactor = 6000, verbose = T)


visium_172C1 <- addStatistics(gobject = visium_172C1)


spatPlot2D(gobject = visium_164C1, show_image = T, point_alpha = 0.7,
           cell_color = 'nr_feats', color_as_factor = F)


visium_164C1 <- calculateHVF(gobject = visium_164C1, save_plot = TRUE)

gene_metadata = fDataDT(visium_164C1)
featgenes = gene_metadata[hvf == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$feat_ID


visium_164C1 <- runPCA(gobject = visium_164C1,
                       genes_to_use = featgenes)

screePlot(visium_164C1, ncp = 30)

plotPCA(gobject = visium_164C1)


visium_164C1 <- runUMAP(visium_164C1, dimensions_to_use = 1:10)
plotUMAP(gobject = visium_164C1)


visium_164C1 <- createNearestNetwork(gobject = visium_164C1, dimensions_to_use = 1:10, k = 15)


visium_164C1 <- doLeidenCluster(gobject = visium_164C1, resolution = 0.4, n_iterations = 1000)

plotUMAP(gobject = visium_164C1,
         cell_color = 'leiden_clus', show_NN_network = F, point_size = 2.5)


spatDimPlot(gobject = visium_164C1, cell_color = 'leiden_clus',
            dim_point_size = 2, spat_point_size = 2.5)

spatDimPlot(gobject = visium_164C1, cell_color = 'nr_feats', color_as_factor = F,
            dim_point_size = 2, spat_point_size = 2.5)


spatPlot2D(visium_164C1, cell_color = 'leiden_clus',
           coord_fix_ratio = 1)

spatPlot2D(visium_164C1, cell_color = 'leiden_clus',
           group_by = 'leiden_clus', coord_fix_ratio = 1,
           cow_n_col = 4, show_legend = F,
           save_param = list(base_width = 14, base_height = 14))


spatPlot2D(visium_164C1, cell_color = 'leiden_clus',
           select_cell_groups = '2', coord_fix_ratio = 1, show_other_cells = TRUE,
           cell_color_code = c('2' = 'red'), other_cell_color = "grey", other_point_size = 1.5,
           save_param = list(base_width = 7, base_height = 7))


gini_markers_subclusters = findMarkers_one_vs_all(gobject = visium_164C1,
                                                  method = 'gini',
                                                  expression_values = 'normalized',
                                                  cluster_column = 'leiden_clus',
                                                  min_feats = 20,
                                                  min_expr_gini_score = 0.5,
                                                  min_det_gini_score = 0.5)
topgenes_gini = gini_markers_subclusters[, head(.SD, 2), by = 'cluster']$feats


violinPlot(visium_164C1, feats = unique(topgenes_gini), cluster_column = 'leiden_clus',
           strip_text = 8, strip_position = 'right',
           save_param = list(base_width = 5, base_height = 10))


plotMetaDataHeatmap(visium_164C1, selected_feats = unique(topgenes_gini),
                    metadata_cols = c('leiden_clus'),
                    x_text_size = 10, y_text_size = 10)


dimFeatPlot2D(visium_164C1, expression_values = 'scaled',
              feats = gini_markers_subclusters[, head(.SD, 1), by = 'cluster']$feats,
              cow_n_col = 3, point_size = 0.75,
              save_param = list(base_width = 8, base_height = 8))


scran_markers_subclusters = findMarkers_one_vs_all(gobject = visium_164C1,
                                                   method = 'scran',
                                                   expression_values = 'normalized',
                                                   cluster_column = 'leiden_clus')
topgenes_scran = scran_markers_subclusters[, head(.SD, 2), by = 'cluster']$feats


violinPlot(visium_164C1, feats = unique(topgenes_scran), cluster_column = 'leiden_clus',
           strip_text = 10, strip_position = 'right',
           save_param = list(base_width = 5))


plotMetaDataHeatmap(visium_164C1, selected_feats = topgenes_scran,
                    metadata_cols = c('leiden_clus'))


dimFeatPlot2D(visium_164C1, expression_values = 'scaled',
              feats = scran_markers_subclusters[, head(.SD, 1), by = 'cluster']$feats,
              cow_n_col = 3, point_size = 1,
              save_param = list(base_width = 8, base_height = 8))









hemato_markers = c("Gata1", "Tal1", "Runx1", "Klf1", "Nfe2l2", "Elf1", "Lmo2", "Nfe2", "Creg1", "Ikzf1")

neuron_markers = c("Sox2", "Sox21", "Pax6", "Pou3f3", "En2", "Lhx5", "Lhx2", "Pou3f2", "Fez1")

epithe_markers = c("Epcam", "Cldn6", "Cldn7", "Cdh1", "Tff3", "Grhl2", "Krt8", "Krt18", "Spint2")

mesoderm_markers = c("Col3a1", "Pcolce", "Cdh11",  "Cnn2", "Serpinh1", "Rhoc", "Maged2", "Fbn2", "Col1a2", "Sparc")

PAGE_matrix_1 = makeSignMatrixPAGE(sign_names = c('Hematopoietic_cells',
                                                  'Neuronal_cells',
                                                  'epithelial_cells',
                                                  'mesoderm_derived cells'),
                                   sign_list = list(hemato_markers,
                                                    neuron_markers,
                                                    epithe_markers,
                                                    mesoderm_markers))

visium_Spleen = runPAGEEnrich(gobject = visium_Spleen, sign_matrix = PAGE_matrix_1)

cell_types_PAGE = colnames(PAGE_matrix_1)
plotMetaDataCellsHeatmap(gobject = visium_Spleen,
                         metadata_cols = 'leiden_clus',
                         value_cols = cell_types_PAGE,
                         spat_enr_names = 'PAGE',
                         x_text_size = 8,
                         y_text_size = 8)


spatCellPlot2D(gobject = visium_Spleen,
               spat_enr_names = 'PAGE',
               cell_annotation_values = cell_types_PAGE[1:4],
               cow_n_col = 2,coord_fix_ratio = 1, point_size = 1.25, show_legend = T)


spatDimCellPlot2D(gobject = visium_Spleen,
                  spat_enr_names = 'PAGE',
                  cell_annotation_values = cell_types_PAGE[1:4],
                  cow_n_col = 1, spat_point_size = 1,
                  plot_alignment = 'horizontal',
                  save_param = list(base_width=7, base_height=10))

visium_Spleen <- createSpatialGrid(gobject = visium_Spleen,
                                   sdimx_stepsize = 400,
                                   sdimy_stepsize = 400,
                                   minimum_padding = 0)

showGiottoSpatGrids(visium_Spleen)

spatPlot2D(visium_Spleen, cell_color = 'leiden_clus', show_grid = T,
           grid_color = 'red', spatial_grid_name = 'spatial_grid')

visium_Spleen <- createSpatialNetwork(gobject = visium_Spleen,
                                      method = 'kNN', k = 5,
                                      maximum_distance_knn = 400,
                                      name = 'spatial_network')

showGiottoSpatNetworks(visium_Spleen)

spatPlot2D(gobject = visium_Spleen,  show_network= T,
           network_color = 'blue', spatial_network_name = 'spatial_network')


ranktest = binSpect(visium_Spleen, bin_method = 'rank',
                    calc_hub = T, hub_min_int = 5,
                    spatial_network_name = 'spatial_network')

spatFeatPlot2D(visium_Spleen, expression_values = 'scaled',
               feats = ranktest$feats[1:6], cow_n_col = 2, point_size = 1.5)



ext_spatial_genes = ranktest[1:1500,]$feats


spat_cor_netw_DT = detectSpatialCorFeats(visium_Spleen,
                                         method = 'network',
                                         spatial_network_name = 'spatial_network',
                                         subset_feats = ext_spatial_genes)


top10_genes = showSpatialCorFeats(spat_cor_netw_DT, feats = 'Car2', show_top_feats = 10)

spatFeatPlot2D(visium_Spleen, expression_values = 'scaled',
               feats = top10_genes$variable[1:4], point_size = 3)


spat_cor_netw_DT = clusterSpatialCorFeats(spat_cor_netw_DT, name = 'spat_netw_clus', k = 20)


heatmSpatialCorFeats(visium_Spleen,
                     spatCorObject = spat_cor_netw_DT,
                     use_clus_name = 'spat_netw_clus',
                     heatmap_legend_param = list(title = NULL),
                     save_param = list(base_height = 6, base_width = 8, units = 'cm'))


netw_ranks = rankSpatialCorGroups(visium_Spleen,
                                  spatCorObject = spat_cor_netw_DT, use_clus_name = 'spat_netw_clus',
                                  save_param = list(  base_height = 3, base_width = 5))



top_netw_spat_cluster = showSpatialCorFeats(spat_cor_netw_DT, use_clus_name = 'spat_netw_clus',
                                            selected_clusters = 6, show_top_feats = 1)



cluster_genes_DT = showSpatialCorFeats(spat_cor_netw_DT, use_clus_name = 'spat_netw_clus', show_top_feats = 1)
cluster_genes = cluster_genes_DT$clus; names(cluster_genes) = cluster_genes_DT$feat_ID

visium_Spleen = createMetafeats(visium_Spleen, feat_clusters = cluster_genes, name = 'cluster_metagene')

showGiottoSpatEnrichments(visium_Spleen)

spatCellPlot(visium_Spleen,
             spat_enr_names = 'cluster_metagene',
             cell_annotation_values = netw_ranks$clusters,
             point_size = 1, cow_n_col = 5, save_param = list(base_width = 15))


table(spat_cor_netw_DT$cor_clusters$spat_netw_clus)
coexpr_dt = data.table::data.table(genes = names(spat_cor_netw_DT$cor_clusters$spat_netw_clus),
                                   cluster = spat_cor_netw_DT$cor_clusters$spat_netw_clus)
data.table::setorder(coexpr_dt, cluster)
top30_coexpr_dt = coexpr_dt[, head(.SD, 30) , by = cluster]
my_spatial_genes <- top30_coexpr_dt$genes



visium_Spleen <- runPCA(gobject = visium_Spleen,
                        feats_to_use = my_spatial_genes,
                        name = 'custom_pca')
visium_Spleen <- runUMAP(visium_Spleen, dim_reduction_name = 'custom_pca', dimensions_to_use = 1:20,
                         name = 'custom_umap')
visium_Spleen <- createNearestNetwork(gobject = visium_Spleen,
                                      dim_reduction_name = 'custom_pca',
                                      dimensions_to_use = 1:20, k = 5,
                                      name = 'custom_NN')
visium_Spleen <- doLeidenCluster(gobject = visium_Spleen, network_name = 'custom_NN',
                                 resolution = 0.15, n_iterations = 1000,
                                 name = 'custom_leiden')


cell_meta = pDataDT(visium_Spleen)
cell_clusters = unique(cell_meta$custom_leiden)

selected_colors = getDistinctColors(length(cell_clusters))
names(selected_colors) = cell_clusters

spatPlot2D(visium_Spleen, cell_color = 'custom_leiden', cell_color_code = selected_colors, coord_fix_ratio = 1)


plotUMAP(gobject = visium_Spleen, cell_color = 'custom_leiden', cell_color_code = selected_colors, point_size = 1.5)



hmrf_folder = paste0(results_folder,'/','11_HMRF/')
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)

HMRF_spatial_genes = doHMRF(gobject = visium_Spleen,
                            expression_values = 'scaled',
                            spatial_genes = my_spatial_genes, k = 20,
                            spatial_network_name="spatial_network",
                            betas = c(0, 10, 5),
                            output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_topgenes_k20_scaled'))

visium_Spleen = addHMRF(gobject = visium_Spleen, HMRFoutput = HMRF_spatial_genes,
                        k = 20, betas_to_add = c(0, 10, 20, 30, 40),
                        hmrf_name = 'HMRF')

spatPlot2D(gobject = visium_Spleen, cell_color = 'HMRF_k20_b.40')



set.seed(seed = 2841)
cell_proximities = cellProximityEnrichment(gobject = visium_Spleen,
                                           cluster_column = 'leiden_clus',
                                           spatial_network_name = 'spatial_network',
                                           adjust_method = 'fdr',
                                           number_of_simulations = 1000)
# barplot
cellProximityBarplot(gobject = visium_Spleen,
                     CPscore = cell_proximities,
                     min_orig_ints = 5, min_sim_ints = 5, p_val = 0.5)

## heatmap
cellProximityHeatmap(gobject = visium_Spleen, CPscore = cell_proximities,
                     order_cell_types = T, scale = T,
                     color_breaks = c(-1.5, 0, 1.5),
                     color_names = c('blue', 'white', 'red'))

# network
cellProximityNetwork(gobject = visium_Spleen, CPscore = cell_proximities,
                     remove_self_edges = F, only_show_enrichment_edges = T)


# network with self-edges
cellProximityNetwork(gobject = visium_Spleen, CPscore = cell_proximities,
                     remove_self_edges = F, self_loop_strength = 0.3,
                     only_show_enrichment_edges = F,
                     rescale_edge_weights = T,
                     node_size = 8,
                     edge_weight_range_depletion = c(1, 2),
                     edge_weight_range_enrichment = c(2,5))


# Option 1
spec_interaction = "cluster 4--cluster 6"
cellProximitySpatPlot2D(gobject = visium_Spleen,
                        interaction_name = spec_interaction,
                        show_network = T,
                        cluster_column = 'leiden_clus',
                        cell_color = 'leiden_clus',
                        cell_color_code = c('cluster 4' = 'lightblue', 'cluster 6' = 'red'),
                        point_size_select = 4, point_size_other = 2)


## select top 25th highest expressing genes
gene_metadata = fDataDT(visium_Spleen)
plot(gene_metadata$nr_cells, gene_metadata$mean_expr)
plot(gene_metadata$nr_cells, gene_metadata$mean_expr_det)

quantile(gene_metadata$mean_expr_det)
high_expressed_genes = gene_metadata[mean_expr_det > 4]$feat_ID

## identify genes that are associated with proximity to other cell types
ICFscoresHighGenes =  findICF(gobject = visium_Spleen,
                              selected_feats = high_expressed_genes,
                              spatial_network_name = 'spatial_network',
                              cluster_column = 'leiden_clus',
                              diff_test = 'permutation',
                              adjust_method = 'fdr',
                              nr_permutations = 500,
                              do_parallel = T)

## visualize all genes
plotCellProximityGenes(visium_Spleen, cpgObject = ICFscoresHighGenes, method = 'dotplot')



## filter genes
ICFscoresFilt = filterICF(ICFscoresHighGenes, min_cells = 2, min_int_cells = 2, min_fdr = 0.1,
                          min_spat_diff = 0.1, min_log2_fc = 0.1, min_zscore = 1)

## visualize subset of interaction changed genes (ICGs)
ICF_genes = c('Hbb-bs', 'Cd74', 'Eef1a1', 'H2-Aa', 'Tpt1')
ICF_genes_types = c('cluster 6', 'cluster 5', 'cluster 1', 'cluster 3', 'cluster 6')
names(ICF_genes) = ICF_genes_types

plotICF(gobject = visium_Spleen,
        cpgObject = ICFscoresHighGenes,
        source_type = 'cluster 1',
        source_markers = 'Cd74',
        ICF_feats = ICF_genes)


LR_data = data.table::fread("~/Desktop/Visium/The second eight mouse samples/Adult_spleen/mouse_lr_pair.txt")

LR_data[, ligand_det := ifelse(ligand_gene_symbol %in% visium_Spleen@feat_ID[['rna']], T, F)]
LR_data[, receptor_det := ifelse(receptor_gene_symbol %in% visium_Spleen@feat_ID[['rna']], T, F)]
LR_data_det = LR_data[ligand_det == T & receptor_det == T]
select_ligands = LR_data_det$ligand_gene_symbol
select_receptors = LR_data_det$receptor_gene_symbol

## get statistical significance of gene pair expression changes based on expression ##
expr_only_scores = exprCellCellcom(gobject = visium_Spleen,
                                   cluster_column = 'leiden_clus',
                                   random_iter = 50,
                                   feat_set_1 = select_ligands,
                                   feat_set_2 = select_receptors)

## get statistical significance of gene pair expression changes upon cell-cell interaction
spatial_all_scores = spatCellCellcom(visium_Spleen,
                                     spat_unit = 'cluster',
                                     feat_type = 'rna',
                                     spatial_network_name = 'cellProximityNetwork',
                                     cluster_column = 'leiden_clus',
                                     random_iter = 50,
                                     feat_set_1 = select_ligands,
                                     feat_set_2 = select_receptors,
                                     adjust_method = 'fdr',
                                     do_parallel = T,
                                     cores = 4,
                                     verbose = 'none')

## * plot communication scores ####

## select top LR ##
selected_spat = spatial_all_scores[p.adj <= 0.5 & abs(log2fc) > 0.1 & lig_nr >= 2 & rec_nr >= 2]
data.table::setorder(selected_spat, -PI)

top_LR_ints = unique(selected_spat[order(-abs(PI))]$LR_comb)[1:33]
top_LR_cell_ints = unique(selected_spat[order(-abs(PI))]$LR_cell_comb)[1:33]

plotCCcomHeatmap(gobject = visium_Spleen,
                 comScores = spatial_all_scores,
                 selected_LR = top_LR_ints,
                 selected_cell_LR = top_LR_cell_ints,
                 show = 'LR_expr')
