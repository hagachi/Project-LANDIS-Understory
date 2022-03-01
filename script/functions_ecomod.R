theme_Publication <- function(base_size=20) {
  library(grid)
  library(ggthemes)
  (theme_few(base_size=base_size)
    + theme(plot.title = element_text(
      size = rel(1.2), hjust = 0.5),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = 'black'),
      axis.title = element_text(size = rel(1)),
      axis.title.y = element_text(angle=90,vjust =2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(), 
      axis.line = element_line(colour="black"),
      axis.ticks = element_line(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid = element_blank(),
      plot.margin=unit(c(10,5,5,5),"mm"),
      strip.background=element_rect(colour=NA,fill=NA),
      strip.text = element_text()
    ))
}


ReadYears <- function(x) {
  management = x[1]
  climate = x[2]
  region = x[3]
  nrep = x[4]
  case_fname = file.path(result_dir, paste0('data_000_biom_', management, '_', climate, '_', region, '_', nrep, '.parquet'))
  if (file.exists(case_fname)) {
    case_df <- arrow::read_parquet(case_fname)
  } else {
    case_lst <- list()
    i_biom <- 1
    reg_dir <- file.path(root_dir, paste(wtcase, management, sep = '_'), climate, paste0('reg', region, '_iter', nrep))
    # Read Maps
    # 1. initial total biomass map
    total_biom_init_df <- as.data.frame(raster(file.path(reg_dir, 'OutputMaps/biomass/TotalBiomass-0.tif')))
    biom_init_df <- as.data.frame(stack(file.path(reg_dir, paste0('OutputMaps/biomass/', sppnames, '-0.tif'))))
    ids <- 1:nrow(total_biom_init_df)
    # 2. forest type
    ftype_df <- as.data.frame(raster(file.path(root_dir, paste0('input/mng', region, '_2021-03-09.tif'))))
    # 3. Ecoregions
    eco_df <- as.data.frame(raster(file.path(root_dir, paste0('input/eco', region, '_2021-03-01.tif'))))
    # Identify forest type for each windthrow sites -----
    # 4. Young trees
    young_fnames <- file.path(reg_dir, paste0('OutputMaps/spp-biomass-by-age/', sppnames, '-ageclass1-0.tif'))
    young_df <- as.data.frame(stack(young_fnames))
    for (yr in years) {
      total_biom_df <- as.data.frame(raster(file.path(reg_dir, paste0('OutputMaps/biomass/TotalBiomass-', yr, '.tif'))))
      biom_df <- as.data.frame(stack(file.path(reg_dir, paste0('OutputMaps/biomass/', sppnames, '-', yr, '.tif'))))
      biom_u25_df <- as.data.frame(stack(file.path(reg_dir, paste0('OutputMaps/spp-biomass-by-age/', sppnames, '-ageclass1-', yr, '.tif'))))
      biom_u50_df <- as.data.frame(stack(file.path(reg_dir, paste0('OutputMaps/spp-biomass-by-age/', sppnames, '-ageclass2-', yr, '.tif'))))
      biom_o50_df <- as.data.frame(stack(file.path(reg_dir, paste0('OutputMaps/spp-biomass-by-age/', sppnames, '-ageclass3-', yr, '.tif'))))
      biom_df <- data.frame(id = ids,
                            t1 = total_biom_init_df,
                            t3 = total_biom_df,
                            biom_init_df,
                            biom_df,
                            biom_u25_df, biom_u50_df, biom_o50_df,
                            young_df,
                            ftype = ftype_df, 
                            eco = eco_df)
      colnames(biom_df) <- c("gridID", 't1', 't3', 
                             paste0('byspp_t1_', sppnames), 
                             paste0('byspp_t3_', sppnames), 
                             paste0('u25_', sppnames), paste0('u50_', sppnames), paste0('o50_', sppnames), 
                             paste0('young_', sppnames),
                             'wt', 'eco')
      case_lst[[i_biom]] <- biom_df %>%
        dplyr::filter(wt > 1, t1 > 0) %>%
        dplyr::mutate(wtcase = wtcase,
                      management = management,
                      climate = climate,
                      nrep = nrep,
                      region = region,
                      res_biom = t3/t1,
                      year = yr)
      i_biom <- i_biom + 1
    }
    case_df <- dplyr::bind_rows(case_lst)
    arrow::write_parquet(case_df, case_fname)
  }
  return(case_df)
}


ReadResultsPbapply <- function(managements, climates, regions, nreps, ncl = 1) {
  case_combinations <- expand.grid(managements, climates, regions, c(1:nreps) - 1)
  if (ncl > 1) {
    library(parallel)
    cl <- makeCluster(ncl)
    clusterExport(cl, c('root_dir', 'result_dir', 'wtcase', 'sppnames', 'years',
                        'raster', 'stack', 'as.data.frame', '%>%'))
    cat('\n Started parallel process...\n')
    res <- pbapply(case_combinations, ReadYears, MARGIN = 1, cl = cl) 
    cat('\n Finished parallel process...\n')
    stopCluster(cl)
  } else {
    res <- pbapply(case_combinations, ReadYears, MARGIN = 1) 
  }
  res_df <- bind_rows(res)
  res_df <- res_df %>% 
    mutate(t1_saplings_agb = apply(select(res_df, starts_with('young_'), -'young_sasa_spp'), MARGIN = 1, FUN = sum),
           t1_saplings_agb_cut = cut(t1_saplings_agb, 
                                     breaks = c(seq(0, 10000, by = 1000), 20000), 
                                     labels = c('<1 kg/m2', '<2 kg/m2', '<3 kg/m2', '<4 kg/m2', '<5 kg/m2', 
                                                '<6 kg/m2', '<7 kg/m2', '<8 kg/m2', '<9 kg/m2', '<10 kg/m2', '<20 kg/m2')),
           t1_sasa_cut = cut(byspp_t1_sasa_spp,
                             breaks = seq(0, 5000, by = 1000), 
                             labels = c('<1 kg/m2', '<2 kg/m2', '<3 kg/m2', '<5 kg/m2', '<5 kg/m2')),
                             # labels = c('<2.5 kg/m2', '<5.0 kg/m2')),
           trt = if_else(wt < 100, 'Plantation Forest', 'Natural Forest'),
           year = year + 2014 - 1)
  write_parquet(res_df, resilience_biom_df_name)
  write_feather(res_df, str_replace(resilience_biom_df_name, '.parquet', '.feather'))
  return(res_df)
}


