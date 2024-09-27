
### Color define
color.list <- list()
## Cycling color
non_cycling_color <- "#CCCCCC"  # Light Gray for non-cycling cells
cycling_color <- "#66c2a5"  # Green for cycling cells
cc.color <- c(cycling_color, non_cycling_color)
color.list[['Phase2']] <- cc.color

### Genotype color
egfr.color <- "#0091CA"
nf1.color <- "#D8423D"
pdgfb.color <- "#55AB55"
genotype.color <- c(egfr.color, nf1.color, pdgfb.color)
color.list[['groupID']] <- genotype.color

if (compartment.name == 'High.level'){
    color.list[['celltype']] <- paletteer_d("ggsci::category20_d3")
    color.list[['sub_celltype']] <- paletteer_c("grDevices::rainbow", dim(table(srt@meta.data[,'sub_celltype'])))
    color.list[['high_level_celltype']] <- paletteer_d("RColorBrewer::Set1")
}


if (compartment.name == 'Tumor'){
    color.list[['celltype']] <- paletteer_d("ggsci::default_locuszoom")
    color.list[['sub_celltype']] <- paletteer_c("grDevices::rainbow", 1)
}

if (compartment.name == 'Myeloid'){
    color.list[['celltype']] <- paletteer_d("RColorBrewer::Set3")
}

if (compartment.name == 'Microglia'){
    color.list[['celltype']] <- paletteer_d("ggthemes::Classic_10_Light")
    color.list[['sub_celltype']] <- paletteer_d("ggthemes::Classic_10_Light")
}

if (compartment.name == 'MDM.Monocytes'){
    color.list[['sub_celltype']] <- paletteer_d("tidyquant::tq_dark")
    color.list[['celltype']] <- paletteer_d("tidyquant::tq_dark")[c(length(paletteer_d("tidyquant::tq_dark"))-2, length(paletteer_d("tidyquant::tq_dark"))-1, length(paletteer_d("tidyquant::tq_dark")))]
}

if (compartment.name == 'Neutrophils'){
    color.list[['sub_celltype']] <- paletteer_d("ggsci::default_nejm")
}

if (compartment.name == 'Bplasma'){
    color.list[['celltype']] <- paletteer_d("ggsci::default_jama")
}


if (compartment.name == 'TNK'){
    color.list[['sub_celltype']] <- paletteer_d("ggthemes::stata_s1color")
}


if (compartment.name == 'Stroma'){
    color.list[['sub_celltype']] <- paletteer_d("ggsci::category20b_d3")
}


if (compartment.name == 'Endothelial'){
    color.list[['sub_celltype']] <- paletteer_d("ggthemes::Classic_20")
}

if (compartment.name == 'Astro.Ologi.Neuro'){
    color.list[['celltype']] <- paletteer_d("ggthemes::excel_Marquee")
}

