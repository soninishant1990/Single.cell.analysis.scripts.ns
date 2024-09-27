
### Cellchat main function
## It generates the cellchat object for all objects seperately and generate the basic figure, like chord plot, bubble plot etc
Cellchat.fun <- function(srt = NULL, Genotype = NULL, out_file1=NULL, force = TRUE, 
organism = NULL, pathways.show.list = NULL){
    for (Genotype.name in unique(srt@meta.data[, Genotype])){
      #Genotype.name = 'CTRL'
    message(Genotype.name)
    #Genotype.name = 'Nf1'
    out_file <- paste0(out_file1,Genotype.name, '/')
    dir.create(out_file)
    file.name = paste0(out_file, Genotype.name, '.cellchat.object.file.rds')
    print(paste0('file name = ', file.name))
    if (!file.exists (file.name) | force) {
      ## check with this range(srt@assays$RNA@data)
      data.input = srt@assays$RNA@data
      meta = srt@meta.data

      cell.use = rownames(meta)[meta[,Genotype] == Genotype.name] # extract the cell names from disease data


      # Prepare input data for CelChat analysis
      data.input = data.input[, cell.use]
      meta = meta[cell.use, ]

      #unique(meta[, Annotation_level_col]) # check the cell labels

      ### Create a CellChat object
      cellchat <- createCellChat(object = data.input, meta = meta, group.by = Annotation_level_col)


      ## Add cell information into meta slot of the object (Optional)
      cellchat <- addMeta(cellchat, meta = meta)
      cellchat <- setIdent(cellchat, ident.use = Annotation_level_col) # set "labels" as default cell identity
      #levels(cellchat@idents) # show factor levels of the cell labels
      groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group


      ### Set the ligand-receptor interaction database
      ## for mouse = CellChatDB.mouse
      ## for human = CellChatDB.human
      if (organism == 'Mouse'){
          CellChatDB <- CellChatDB.mouse
      } else if (organism == 'Human'){
          CellChatDB <- CellChatDB.human
      } else{
          print('Provide organism type, Human or Mouse')
      }

      ## Save cellchat database db Category
      cellchat.database.db.Category <- showDatabaseCategory(CellChatDB)
      png (paste0(out_file,'Cellchat.database.db.Category.png'), width = 1500, height = 1500, pointsize=10, res = 300, type="cairo")
      print (wrap_plots (cellchat.database.db.Category))
      dev.off()

      # Show the structure of the database
      #dplyr::glimpse(CellChatDB$interaction)

      # use a subset of CellChatDB for cell-cell communication analysis
      # CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
      # use all CellChatDB for cell-cell communication analysis
       CellChatDB.use <- CellChatDB # simply use the default CellChatDB

      # set the used database in the object
      cellchat@DB <- CellChatDB.use


      ### Preprocessing the expression data for cell-cell communication analysis
      # subset the expression data of signaling genes to save computation cost
      cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
      future::plan("multisession", workers = 4) # do parallel
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat)

      ### Part II: Inference of cell-cell communication network
      ### Compute the communication probability and infer cellular communication network
      cellchat <- computeCommunProb(cellchat, type = "triMean")

      ## Users can filter out cell-cell communication if there are only a few cells in certain cell groups. By default, the minimum number of cells required in each cell group for cell-cell communication is 10.
      cellchat <- filterCommunication(cellchat, min.cells = 10)

      ## Infer the cell-cell communication at a signaling pathway level
      cellchat <- computeCommunProbPathway(cellchat)

      ## Calculate the aggregated cell-cell communication network
      cellchat <- aggregateNet(cellchat)

      ## save cellchat object
      saveRDS(cellchat, file = file.name)
    } #else{
      cellchat <- readRDS(file.name)
      print(cellchat)
      Visualize.the.aggregated.cell.cell.communication.network(cellchat = cellchat, out_file = out_file, 
      Genotype.name = Genotype.name)
      Visualization.of.cell.cell.communication.network(cellchat = cellchat, out_file = out_file, 
      Genotype.name = Genotype.name, pathways.show.list = pathways.show.list)
      Visualize.cell.cell.communication.mediated.by.multiple.ligand.receptors.or.signaling.pathways(cellchat = cellchat, out_file = out_file, 
      Genotype.name = Genotype.name)
    #}
    
  }
}


Visualize.the.aggregated.cell.cell.communication.network <- function(cellchat = NULL, out_file = NULL, 
      Genotype.name = NULL){
    out_fig <- paste0(out_file, 'Visualize.the.aggregated.cell-cell.communication.network/')
    dir.create(out_fig)
    # Visualize the aggregated cell-cell communication network
    groupSize <- as.numeric(table(cellchat@idents))
    par(mfrow = c(1, 2), xpd = TRUE)

    # Plot based on Number.of.interactions
    pdf(paste0(out_fig, Genotype.name, '.Number.of.interactions.pdf'), width = 12, height = 12)
    p1 <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
    grid.text("Number of interactions", x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
    print(p1)
    dev.off()

    # Plot based on Interaction.weights.strength
    pdf(paste0(out_fig, Genotype.name, '.Interaction.weights.strength.pdf'), width = 12, height = 12)
    p1 <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")
    grid.text("Interaction weights/strength", x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
    print(p1)
    dev.off()

    # Network Plots
    mat <- cellchat@net$weight
    ggplot_list <- list()
    mat.row <- rownames(mat)
    pdf(paste0(out_fig, Genotype.name, '.Interaction.weights.strength.at.celltype.pdf'), width = 12, height = 12)
    for (i in 1:nrow(mat)) {
      mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      mat2[i, ] <- mat[i, ]
      # Add title directly in netVisual_circle
      p1 <- netVisual_circle(
        mat2,
        vertex.weight = groupSize,
        weight.scale = TRUE,
        edge.weight.max = max(mat),
        title.name = paste0(i, ". Interaction.weights.strength")
      )
      grid.text(paste0(mat.row[[i]], ". Interaction.weights.strength"), x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
      ggplot_list[[mat.row[[i]]]] <- p1
    }
    # Arrange and print the grobs in a grid
    print(cowplot::plot_grid(plotlist = ggplot_list, ncol = 4))
    dev.off()

    # Network Plots
    mat <- cellchat@net$count
    ggplot_list <- list()
    mat.row <- rownames(mat)
    pdf(paste0(out_fig, Genotype.name, '.Number.of.interactions.at.celltype.pdf'), width = 12, height = 12)
    for (i in 1:nrow(mat)) {
      mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      mat2[i, ] <- mat[i, ]
      # Add title directly in netVisual_circle
      p1 <- netVisual_circle(
        mat2,
        vertex.weight = groupSize,
        weight.scale = TRUE,
        edge.weight.max = max(mat),
        title.name = paste0(i, ". Number of interactions")
      )
      grid.text(paste0(mat.row[[i]], ". Number of interactions"), x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
      ggplot_list[[mat.row[[i]]]] <- p1
    }
    # Arrange and print the grobs in a grid
    print(cowplot::plot_grid(plotlist = ggplot_list, ncol = 4))
    dev.off()
}


#Visualization.of.cell.cell.communication.network <- function(out_file = NULL, srt = NULL, high_level_celltype_col=NULL){
Visualization.of.cell.cell.communication.network <- function(cellchat = NULL, out_file = NULL, 
      Genotype.name = NULL, pathways.show.list = NULL){
  out_file.fig <- paste0(out_file, 'Visualization.of.cell.cell.communication.network/')
  dir.create(out_file.fig)
  for (pathways.show in pathways.show.list){

    tryCatch({
    out_fig <- paste0(out_file.fig, pathways.show,'.gene.communication.network/')
    dir.create(out_fig)
    # Hierarchy plot
    pdf(paste0(out_fig, 'Hierarchy.plot.Visualization.of.cell-cell.communication.network.only.',pathways.show,'.gene.LR.pdf'), width = 8, height = 8)
    # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
    vertex.receiver = seq(1,4) # a numeric vector.
    p1 <- netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
    grid.text(paste0('Hierarchy plot Visualization of cell-cell communication network only ',pathways.show, ' Gene LR'), x = 0.5, y = 0.95, gp = gpar(fontsize = 8, fontface = "bold"))
    print(p1)
    dev.off()
    
    
    # Circle plot
    pdf(paste0(out_fig, 'Circle.plot.Visualization.of.cell-cell.communication.network.only.',pathways.show,'.gene.LR.pdf'), width = 8, height = 8)
    par(mfrow=c(1,1))
    p1 <- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
    grid.text(paste0('Circle plot Visualization of cell-cell communication network only ',pathways.show, ' Gene LR'), x = 0.5, y = 0.95, gp = gpar(fontsize = 8, fontface = "bold"))
    print(p1)
    dev.off()
    
    
    # Chord diagram
    pdf(paste0(out_fig, 'Chord.diagram.Visualization.of.cell-cell.communication.network.only.',pathways.show,'.gene.LR.pdf'), width = 8, height = 8)
    par(mfrow=c(1,1))
    p1 <- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
    #grid.text(paste0('Chord diagram Visualization of cell-cell communication network only ',pathways.show, ' Gene LR'), x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
    print(p1)
    dev.off()
    
    
    # Heatmap
    pdf(paste0(out_fig, 'Heatmap.Visualization.of.cell-cell.communication.network.only.',pathways.show,'.gene.LR.pdf'), width = 8, height = 8)
    par(mfrow=c(1,1))
    p1 <- netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
    #grid.text(paste0('Heatmap Visualization of cell-cell communication network only ',pathways.show, ' Gene LR'), x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
    print(p1)
    dev.off()
    #> Do heatmap based on a single object
    
    
    if(!is.null(high_level_celltype_col)){
      # Assuming srt is your data frame
      unique_celltype <- unique(srt@meta.data[,Annotation_level_col])
      # Create an empty list to store key-value pairs
      result_list <- list()
      # Iterate over unique high_level_celltypes
      for (key in unique_celltype) {
        # Get corresponding values from celltype
        values <- srt@meta.data[,high_level_celltype_col][srt@meta.data[,Annotation_level_col] == key]
        # Remove duplicates if needed
        values <- unique(values)
        # Assign to the result list
        result_list[[key]] <- values
      }
      pdf(paste0(out_fig, 'Chord.diagram.Visualization.of.cell-cell.communication.network.only.',pathways.show,'.gene.LR_Group_by_higher_celltype.pdf'), width = 8, height = 8)
      p1 <- netVisual_chord_cell(cellchat, signaling = pathways.show, group = unlist(result_list), title.name = paste0(pathways.show, " signaling network"))
      print(p1)
      dev.off()
    }
    
    
    ### Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
    pdf(paste0(out_fig, 'Compute.the.contribution.of.each.ligand-receptor.pair.to.the.overall.signaling.pathway.and.visualize.cell-cell.communication.mediated.by.a.single.ligand-receptor.pair.only.',pathways.show,'.pdf'), width = 8, height = 8)
    p1 <- netAnalysis_contribution(cellchat, signaling = pathways.show)
    print(p1)
    dev.off()
    
    # We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
    pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
    LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
    
    # Hierarchy plot
    vertex.receiver = seq(1,4) # a numeric vector
    #netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
    
    
    # Circle plot
    pdf(paste0(out_fig, 'Circle.plot.Visualization.of.cell-cell.communication.network.only.highest.interaction.pair.',LR.show,'.gene.LR.pdf'), width = 8, height = 8)
    p1 <- netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
    print(p1)
    dev.off()
    
    
    # Chord diagram
    pdf(paste0(out_fig, 'Chord.diagram.Visualization.of.cell-cell.communication.network.only.highest.interaction.pair.',LR.show,'.gene.LR.pdf'), width = 8, height = 8)
    p1 <- netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
    print(p1)
    dev.off()

  }, error = function(e) {
        # Handle the error gracefully
        cat("Error occurred for pathway ", pathways.show, ":", conditionMessage(e), "\n")
        # You can choose to skip this iteration or take other appropriate actions
  })
  }
}



### Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
Visualize.cell.cell.communication.mediated.by.multiple.ligand.receptors.or.signaling.pathways <- function(cellchat = NULL, out_file = NULL, 
      Genotype.name = NULL){
    out_fig <- paste0(out_file, 'Visualize.cell-cell.communication.mediated.by.multiple.ligand-receptors.or.signaling.pathways/')
    dir.create(out_fig)
    
    ## Bubble plot Visualize cell-cell communication mediated by multiple ligand-receptors
    mat <- levels(cellchat@idents)
    ggplot_list <- list()
    #mat.row <- rownames(mat)
    pdf(paste0(out_fig, 'Bubble.plot.Visualize.cell-cell.communication.mediated.by.multiple.ligand-receptors.pdf'), width = 10, height = 10)
    for (i in 1:length(mat)) {
        tryCatch({
            p1 <- netVisual_bubble(cellchat, sources.use = i, remove.isolate = FALSE)
            ggplot_list[[i]] <- p1
        }, error = function(e) {
            # Handle the error gracefully
            cat("Error occurred for source", i, ":", conditionMessage(e), "\n")
            # You can choose to skip this iteration or take other appropriate actions
        })
    }
    # Arrange and print the grobs in a grid
    #print(cowplot::plot_grid(plotlist = ggplot_list, ncol = 4))
    print(ggplot_list)
    dev.off()
    
    ## Chord diagram Visualize all the significant interactions ligand-receptors sending from celltype
    mat <- levels(cellchat@idents)
    ggplot_list <- list()
    mat.row <- rownames(mat)
    pdf(paste0(out_fig, 'Chord.diagram.Visualize.all.the.significant.interactions.ligand-receptors_sending.from.celltype.pdf'), width = 10, height = 10)
    for (i in 1:length(mat)) {
        tryCatch({
            p1 <- netVisual_chord_gene(cellchat, sources.use = i, lab.cex = 0.5,legend.pos.y = 30)
            #netVisual_chord_gene(cellchat, sources.use = i, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)
            ### Define target.use for getting desire interections
            #grid.text(i, x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
            ggplot_list[[i]] <- p1
        }, error = function(e) {
            # Handle the error gracefully
            cat("Error occurred for source", i, ":", conditionMessage(e), "\n")
            # You can choose to skip this iteration or take other appropriate actions
        })
    }
    # Arrange and print the grobs in a grid
    #print(cowplot::plot_grid(plotlist = ggplot_list, ncol = 4))
    print(ggplot_list)
    dev.off()
    
    
    ## Chord diagram Visualize all the significant interactions ligand-receptors received from celltype
    mat <- levels(cellchat@idents)
    ggplot_list <- list()
    mat.row <- rownames(mat)
    pdf(paste0(out_fig, 'Chord.diagram.Visualize.all.the.significant.interactions.ligand-receptors_received.from.celltype.pdf'), width = 10, height = 10)
    for (i in 1:length(mat)) {
        tryCatch({
            p1 <- netVisual_chord_gene(cellchat, sources.use = 1:length(mat), targets.use = i, lab.cex = 0.5,legend.pos.y = 30)
            #netVisual_chord_gene(cellchat, sources.use = i, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)
            ### Define target.use for getting desire interections
            #grid.text(i, x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
            ggplot_list[[i]] <- p1
        }, error = function(e) {
            # Handle the error gracefully
            cat("Error occurred for source", i, ":", conditionMessage(e), "\n")
            # You can choose to skip this iteration or take other appropriate actions
        })
    }
    # Arrange and print the grobs in a grid
    #print(cowplot::plot_grid(plotlist = ggplot_list, ncol = 4))
    print(ggplot_list)
    dev.off()
    
}




genotype_com_fun <- function(out_file1 = NULL, out_file2 = NULL, Genotype = NULL, 
force = TRUE, pathways.show.list = NULL){
cellchat.file.name = paste0(out_file1, "cellchat_object.list.rds")
if (!file.exists (cellchat.file.name) | force) {
  ### Pairwise comparision between genotype
  ### Get celltype list across the sample for liftup cellchat
  level.names <- list()
  for (Genotype.name in unique(srt@meta.data[, Genotype])){
    message(Genotype.name)
    out_file.path <- paste0(out_file1,Genotype.name, '/')
    file.name = readRDS(paste0(out_file.path, Genotype.name, '.cellchat.object.file.rds'))
    group.new = levels(file.name@idents)
    level.names <- c(level.names, group.new)
  }
  level.name.list <- unique(level.names)

  ### Merge all genotype objects into one
  cellchatobj.list <- list()
  for (Genotype.name in unique(srt@meta.data[, Genotype])){
    message(Genotype.name)
    #Genotype.name = 'Nf1'
    out_file.path <- paste0(out_file1,Genotype.name, '/')
    file.name = readRDS(paste0(out_file.path, Genotype.name, '.cellchat.object.file.rds'))
    file.name <- liftCellChat(file.name, level.name.list)
    cellchatobj.list[[Genotype.name]] <- file.name
  }
  print(cellchatobj.list)
  cellchat <- mergeCellChat(cellchatobj.list, add.names = names(cellchatobj.list))


  ## Save cellchat both object list object and merge object
  # Users can now export the merged CellChat object and the list of the two separate objects for later use
  save(cellchatobj.list, file = paste0(out_file1, "cellchat_object.list.rds"))
  save(cellchat, file = paste0(out_file1, "cellchat_merged.rds"))
} else{
  message('Reading cellchat object')
  load(paste0(out_file1, "cellchat_object.list.rds")) ## load by cellchatobj.list
  load(paste0(out_file1, "cellchat_merged.rds")) ## load by cellchat
}


# Define the genotypes for comparison
genotypes <- unique(srt@meta.data[,Genotype])
# Create a pairwise comparison of genotypes
pair_list <- combn(genotypes, 2, simplify = FALSE)

#pair_list[[1]]
for (pair.name in pair_list){
  #pair.name <- pair_list[[1]]
  message(pair.name)
  pair1 <- pair.name[[1]]
  cellchatobj.list[[pair1]]
  pair2 <- pair.name[[2]]
  cellchatobj.list[[pair2]]

  #creat pair folder
  ## Create new directory for output file
  out_file3 <- paste0(out_file2,pair1, '.vs.', pair2, '/')
  dir.create(out_file3)

  file.name = paste0(out_file3, "cellchat_merged.rds")
  if (!file.exists (file.name) | force) {
    n.cellchatobj.list <- cellchatobj.list[c(pair1, pair2)]
    cellchat <- mergeCellChat(n.cellchatobj.list, add.names = names(n.cellchatobj.list))
    save(cellchat, file = paste0(out_file3, "cellchat_merged.rds"))
  } else{
    n.cellchatobj.list <- cellchatobj.list[c(pair1, pair2)]
    load(paste0(out_file3, "cellchat_merged.rds"))
  }


  ## Part I: Predict general principles of cell-cell communication
  ## CellChat starts with the big picture to predict general principles of cell-cell communication. When comparing cell-cell communication among multiple biological conditions, it can answer the following biological questions:
  ## Whether the cell-cell communication is enhanced or not
  ## The interaction between which cell types is significantly changed
  ## How the major sources and targets change from one condition to another
  ## Compare the total number of interactions and interaction strength

  out_file4 <- paste0(out_file3,'Compare_number_of_interactions_and_interaction_strength/')
  dir.create(out_file4)
  gg1 <- compareInteractions(cellchat, show.legend = F)
  gg2 <- compareInteractions(cellchat, show.legend = F, measure = "weight")
  png (paste0(out_file4,'Compare.the.total.number.of.interactions.and.interaction.strength.png'), width = 1500, height = 1200, pointsize=10, res = 300, type="cairo")
  print (wrap_plots (gg1 + gg2))
  dev.off()


  ## Differential number of interactions or interaction strength among different cell populations
  # Plot based on Number.of.interactions
  pdf(paste0(out_file4, 'Number.of.interactions.and.weights.strength.pdf'), width = 13, height = 12)
  p1 <- netVisual_diffInteraction(cellchat, weight.scale = T, title.name = 'Number of interactions')
  grid.text("Number of interactions", x = 0.2, y = 0.90, gp = gpar(fontsize = 16, fontface = "bold"))
  print(p1)
  p2 <- netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", title.name = 'Interaction weights/strength')
  grid.text("Interaction weights/strength", x = 0.2, y = 0.90, gp = gpar(fontsize = 16, fontface = "bold"))
  print(p2)
  dev.off()

  pdf(paste0(out_file4, 'Do.heatmap.based.on.a.merged.object.pdf'), width = 8, height = 8)
  #> Do heatmap based on a merged object
  p1 <- netVisual_heatmap(cellchat, title.name = 'Number of interactions')
  print(p1)
  #> Do heatmap based on a merged object
  p2 <- netVisual_heatmap(cellchat, measure = "weight", title.name = 'Interaction weights/strength')
  print(p2)
  dev.off()


  ## Seperate circle plot
  pdf(paste0(out_file4, 'Number.of.interactions.at_genotype_level.pdf'), width = 8, height = 8)
  weight.max <- getMaxWeight(n.cellchatobj.list, attribute = c("idents","count"))
  plot.list <- list()
  for (i in 1:length(n.cellchatobj.list)) {
    p1 <- netVisual_circle(n.cellchatobj.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(n.cellchatobj.list)[i]))
    grid.text(paste0("Number of interactions - ", names(n.cellchatobj.list)[i]), x = 0.2, y = 0.88, gp = gpar(fontsize = 7, fontface = "bold"))
    print(p1)
  }
  dev.off()



  # Compute the network centrality scores
  n.cellchatobj.list[[pair1]] <- netAnalysis_computeCentrality(n.cellchatobj.list[[pair1]], slot.name = "netP") # the slot 'netP' means the inferred intercellular c
  n.cellchatobj.list[[pair2]] <- netAnalysis_computeCentrality(n.cellchatobj.list[[pair2]], slot.name = "netP") # the slot 'netP' means the inferred intercellular c


  ### Compare the major sources and targets in 2D space
  num.link <- sapply(n.cellchatobj.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
  weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
  pdf(paste0(out_file4, 'Compare.the.major.sources.and.targets.in.2D.space.pdf'), width = 12, height = 6)
  gg <- list()
  for (i in 1:length(n.cellchatobj.list)) {
    gg[[i]] <- netAnalysis_signalingRole_scatter(n.cellchatobj.list[[i]], title = names(n.cellchatobj.list)[i], weight.MinMax = weight.MinMax)
  }
  #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  print(patchwork::wrap_plots(plots = gg))
  dev.off()

  ## Identify and visualize the conserved and context-specific signaling pathways
  ## By comparing the information flow/interaction strength of each signaling pathway, we can identify signaling pathways (i) turn off, (ii) decrease, (iii) turn on, or (iv) increase by changing their information flow at one condition as compared to another condition.
  out_file4 <- paste0(out_file3,'Compare_the_overall_information_flow_of_each_signaling_pathway/')
  dir.create(out_file4)
  ### Identify and visualize the conserved and context-specific signaling pathways
  ### Compare the overall information flow of each signaling pathway
  gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
  gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
  pdf(paste0(out_file4, 'Compare.the.overall.information.flow.of.each.signaling.pathway.pdf'), width = 12, height = 6)
  print(gg1 + gg2)
  dev.off()


  ### Compare outgoing (or incoming) signaling associated with each cell population
  i = 1
  # Combining all the identified signaling pathways from different datasets 
  pathway.union <- union(n.cellchatobj.list[[i]]@netP$pathways, n.cellchatobj.list[[i+1]]@netP$pathways)
  ht1 = netAnalysis_signalingRole_heatmap(n.cellchatobj.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(n.cellchatobj.list)[i], width = 5, height = 16)
  ht2 = netAnalysis_signalingRole_heatmap(n.cellchatobj.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(n.cellchatobj.list)[i+1], width = 5, height = 16)
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  pdf(paste0(out_file4, 'Compare.outgoing.incoming.signaling.associated.with.each.cell.population.pdf'), width = 8, height = 10)
  print(ht1 + ht2)
  dev.off()

  
  
  ## Part III: Identify the up-regulated and down-regulated signaling ligand-receptor pairs
  ## Identify dysfunctional signaling by comparing the communication probabilities
  ### Identify the up-regulated and down-regulated signaling ligand-receptor pairs
  out_file4 <- paste0(out_file3,'Identify_the_upgulated_and_down.regulated_signaling_ligand.receptor_pairs/')
  dir.create(out_file4)
  cellchat <- setIdent(cellchat, ident.use = Annotation_level_col) # set "labels" as default cell identity
  mat <- levels(cellchat@idents)
  mat.row <- rownames(mat)
  ggplot_list <- list()
  pdf(paste0(out_file4, 'Bubble.plot.Visualize.cell-cell.communication.mediated.by.multiple.ligand-receptors.pdf'), width = 10, height = 10)
  for (i in 1:length(mat)) {
      tryCatch({
          p1 <- netVisual_bubble(cellchat, sources.use = i,  comparison = c(1, 2), angle.x = 45)#, remove.isolate = FALSE)
          #p1 <- netVisual_bubble(cellchat, sources.use = i, remove.isolate = FALSE)
          ggplot_list[[i]] <- p1
      }, error = function(e) {
          # Handle the error gracefully
          cat("Error occurred for source", i, ":", conditionMessage(e), "\n")
          # You can choose to skip this iteration or take other appropriate actions
      })
  }
  print(ggplot_list)
  dev.off()


  ### Part IV: Visually compare cell-cell communication using a Hierarchy plot, Circle plot, or Chord diagram
  out_file4 <- paste0(out_file3, 'Visualization.of.cell.cell.communication.network/')
  dir.create(out_file4)
  for (pathways.show in pathways.show.list){
    tryCatch({

        out_file5 <- paste0(out_file4, pathways.show,'.gene.communication.network/')
        dir.create(out_file5)

        pdf(paste0(out_file5, 'Circle.plot.Visualization.of.cell-cell.communication.network.only.',pathways.show,'.gene.LR.pdf'), width = 8, height = 8)
        weight.max <- getMaxWeight(n.cellchatobj.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
        #par(mfrow = c(1,2), xpd=TRUE)
        plot_list <- list()
        for (i in 1:length(n.cellchatobj.list)) {
          p1 <- netVisual_aggregate(n.cellchatobj.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(n.cellchatobj.list)[i]))
          grid.text(paste0(paste(pathways.show, names(n.cellchatobj.list)[i])), x = 0.2, y = 0.85, gp = gpar(fontsize = 8, fontface = "bold"))
          plot_list[[i]] <- p1
        }
        print(plot_list)
        dev.off()


        pdf(paste0(out_file5, 'heatmap.plot.Visualization.of.cell-cell.communication.network.only.',pathways.show,'.gene.LR.pdf'), width = 8, height = 8)
        ht <- list()
        for (i in 1:length(n.cellchatobj.list)) {
          ht[[i]] <- netVisual_heatmap(n.cellchatobj.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(n.cellchatobj.list)[i]))
          #grid.text(paste0(paste(pathways.show, names(n.cellchatobj.list)[i])), x = 0.2, y = 0.85, gp = gpar(fontsize = 8, fontface = "bold"))
        }
        print(ht)
        dev.off()

        # Chord diagram
        pdf(paste0(out_file5, 'Chord.plot.Visualization.of.cell-cell.communication.network.only.',pathways.show,'.gene.LR.pdf'), width = 8, height = 8)
        cp <- list()
        for (i in 1:length(n.cellchatobj.list)) {
          cp[[i]] <- netVisual_aggregate(n.cellchatobj.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(n.cellchatobj.list)[i]))
        }
        print(cp)
        dev.off()


        if(!is.null(high_level_celltype_col)){
          # Assuming srt is your data frame
          unique_celltype <- unique(srt@meta.data[,Annotation_level_col])
          # Create an empty list to store key-value pairs
          result_list <- list()
          # Iterate over unique high_level_celltypes
          for (key in unique_celltype) {
            # Get corresponding values from celltype
            values <- srt@meta.data[,high_level_celltype_col][srt@meta.data[,Annotation_level_col] == key]
            # Remove duplicates if needed
            values <- unique(values)
            # Assign to the result list
            result_list[[key]] <- values
          }
          pdf(paste0(out_file5, 'Chord.diagram.Visualization.of.cell-cell.communication.network.only.',pathways.show,'.gene.LR_Group_by_higher_celltype.pdf'), width = 8, height = 8)
          #p1 <- netVisual_chord_cell(cellchat, signaling = pathways.show, group = unlist(result_list), title.name = paste0(pathways.show, " signaling network"))
          cp <- list()
          for (i in 1:length(n.cellchatobj.list)) {
            cp[[i]] <- netVisual_chord_cell(n.cellchatobj.list[[i]], signaling = pathways.show, group = unlist(result_list), title.name = paste0(pathways.show, " signaling network - ", names(n.cellchatobj.list)[i]))
          }          
          print(cp)
          dev.off()
        }




      }, error = function(e) {
        # Handle the error gracefully
        cat("Error occurred for source", i, ":", conditionMessage(e), "\n")
        # You can choose to skip this iteration or take other appropriate actions
  })

  }
}

}



Comparison_analysis_of_multiple_datasets <- function(out_file1 = NULL, out_file2 = NULL, 
Genotype = NULL, force = TRUE, pathways.show.list = NULL){
  cellchat.file.name = paste0(out_file1, "cellchat_object.list.rds")
  if (!file.exists (cellchat.file.name) | force) {
    ### Pairwise comparision between genotype
    ### Get celltype list across the sample for liftup cellchat
    level.names <- list()
    for (Genotype.name in unique(srt@meta.data[, Genotype])){
      message(Genotype.name)
      out_file.path <- paste0(out_file1,Genotype.name, '/')
      file.name = readRDS(paste0(out_file.path, Genotype.name, '.cellchat.object.file.rds'))
      group.new = levels(file.name@idents)
      level.names <- c(level.names, group.new)
    }
    level.name.list <- unique(level.names)

    ### Merge all genotype objects into one
    cellchatobj.list <- list()
    for (Genotype.name in unique(srt@meta.data[, Genotype])){
      message(Genotype.name)
      #Genotype.name = 'Nf1'
      out_file.path <- paste0(out_file1,Genotype.name, '/')
      file.name = readRDS(paste0(out_file.path, Genotype.name, '.cellchat.object.file.rds'))
      file.name <- liftCellChat(file.name, level.name.list)
      cellchatobj.list[[Genotype.name]] <- file.name
    }
    print(cellchatobj.list)
    cellchat <- mergeCellChat(cellchatobj.list, add.names = names(cellchatobj.list))


    ## Save cellchat both object list object and merge object
    # Users can now export the merged CellChat object and the list of the two separate objects for later use
    save(cellchatobj.list, file = paste0(out_file1, "cellchat_object.list.rds"))
    save(cellchat, file = paste0(out_file1, "cellchat_merged.rds"))
  } else{
    message('Reading cellchat object')
    load(paste0(out_file1, "cellchat_object.list.rds")) ## load by cellchatobj.list
    load(paste0(out_file1, "cellchat_merged.rds")) ## load by cellchat
  }


  out_file5 <- paste0(Geno_file,'Comparison.analysis.of.multiple.datasets.with.different.cell.type.compositions/')
  dir.create(out_file5)
  # Hierarchy plot
  ## Get the all pathway list in cellchatDB
  # interaction_table <- table(CellChatDB$interaction$pathway_name)
  for (pathways.show in  pathways.show.list){
    tryCatch({
      weight.max <- getMaxWeight(cellchatobj.list, slot.name = c("netP"), attribute = pathways.show)
      plot.list <- list()
      vertex.receiver = seq(1,10)
      ## Visualize the inferred signaling network
      ### Hierarchy plot
      pdf(paste0(out_file5, pathways.show, '_Visualize.the.inferred.signaling.network.at_genotype_level.Hierarchy.plot.pdf'), width = 8, height = 8)
      for (i in 1:length(cellchatobj.list)) {
        p1 <- netVisual_aggregate(cellchatobj.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(cellchatobj.list)[i]))
        grid.text(paste0("Number of interactions - ", pathways.show, " - ",names(cellchatobj.list)[[i]]), x = 0.2, y = 0.88, gp = gpar(fontsize = 7, fontface = "bold"))
        print(p1)
      }
      dev.off()

      # Circle plot
      pdf(paste0(out_file5, pathways.show, '_Visualize.the.inferred.signaling.network.at_genotype_level.Circle.plot.pdf'), width = 8, height = 8)
      for (i in 1:length(cellchatobj.list)) {
        p1 <- netVisual_aggregate(cellchatobj.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(cellchatobj.list)[i]))
        grid.text(paste0("Number of interactions - ", pathways.show, " - ",names(cellchatobj.list)[[i]]), x = 0.2, y = 0.88, gp = gpar(fontsize = 7, fontface = "bold"))
        print(p1)
      }
      dev.off()

      # Chord diagram
      pdf(paste0(out_file5, pathways.show, '_Visualize.the.inferred.signaling.network.at_genotype_level.Chord.diagram.pdf'), width = 8, height = 8)
      for (i in 1:length(cellchatobj.list)) {
        p1 <- netVisual_aggregate(cellchatobj.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(cellchatobj.list)[i]))
        grid.text(paste0("Number of interactions - ", pathways.show, " - ",names(cellchatobj.list)[[i]]), x = 0.2, y = 0.88, gp = gpar(fontsize = 7, fontface = "bold"))
        print(p1)
      }
      dev.off()

    }, error = function(e) {
        # Handle the error gracefully
        cat("Error occurred for source", pathways.show, ":", conditionMessage(e), "\n")
        # You can choose to skip this iteration or take other appropriate actions
    })
  }

}