#' ---
#' title: Analysis for EcolModel
#' date: 2022.01.05
#' author: Haga, Chihiro
#' ---

# 0. load packages ------
# IO
library(feather)
library(arrow)
# Raster data
library(raster)
# Data processing
library(tidyverse)
library(skimr)
library(pbapply)
library(zoo) # rollmean
# Visualization
library(viridis)
library(ggthemes)
library(patchwork)
library(ggalluvial)
library(RColorBrewer)
source('./script/figure_settings.R')


# 0. Define functions ---------
# path
root_dir <- rprojroot::find_root('Alpha_v3.1.0.Rproj')
result_dir <- file.path(root_dir, "result_2022-01-05")
source(file.path(root_dir, 'script', 'functions_ecomod.R'))
if (!dir.exists(result_dir)) dir.create(result_dir)
resilience_biom_df_name <- file.path(result_dir, "data_001_resilience_biom_df.parquet")
# parameters
years <- c(0, 2, seq(8, 38, by = 10), 58, 86) # years <- c(4, 38, 86) # 40 == 2050, 90 == 2099
wtcase <- 'BaU'
wtyear <- 4 # 2017, one year after windthrow
managements <- c('cl', 'slpl1', 'slpl2')
managements_lutable <- data.frame(management = managements,
                                  mng_name = c('CL', 'SL1', 'SL2'))
# Climate names
climates <- c('CSIRO-Mk3-6-0_RCP8.5', 'GFDL-CM3_RCP8.5', 'HadGEM2-ES_RCP8.5', 'MIROC5_RCP2.6', 'noCC')
climate_lookup <- data.frame(climate = climates,
                             climatelabel = c('RCP8.5 CSIRO', 'RCP8.5 GFDL', 'RCP8.5 HadGEM2', 'RCP2.6 MIROC5', 'Current'))
regions <- c(1)
sppnames <- c('betuerma', 'fagucren', 'quercris', 'acermono', 'larikaem', 'abiesach', 'crypjapo', 'sasa_spp')
spplevels <- c('sasa_spp', 'betuerma', 'larikaem', 'crypjapo', 'quercris', 'fagucren', 'acermono', 'abiesach')
eco_lutable <- data.frame(eco = 1:3, eco_name = c('Lowest', 'Middle', 'Highest'))
nreps <- 5
# parallel setting
NCL <- 2^4


# 1. Read data as dataframe ---------
# Current and future AGB for each damaged sites
if (file.exists(resilience_biom_df_name)) {
  resilience_biom_df <- read_parquet(resilience_biom_df_name) %>%
    left_join(eco_lutable, by = 'eco') %>% 
    left_join(managements_lutable, by = 'management')
} else {
  resilience_biom_df <- ReadResultsPbapply(managements, climates, regions, nreps, ncl = NCL) %>%
    left_join(eco_lutable, by = 'eco') %>% 
    left_join(managements_lutable, by = 'management')
}
# Read Sapling biomass for each location
youngs_df_fname <- file.path(result_dir, 'data_002_young_cohorts.parquet')
if (file.exists(youngs_df_fname)) {
  youngs_df <- read_parquet(youngs_df_fname)
} else {
  young_fnames <- paste0(root_dir, '/BaU_cl/noCC/reg1_iter0/OutputMaps/spp-biomass-by-age/', sppnames[1:(length(sppnames) - 1)], '-ageclass1-0.tif')
  youngs <- stack(young_fnames)
  mng <- as.data.frame(raster(paste0(root_dir, '/input/mng1_2021-03-09.tif')))
  colnames(mng) = 'mng'
  youngs_df <- as.data.frame(youngs) %>% 
    bind_cols(mng) %>% 
    pivot_longer(cols = c(-mng)) %>% 
    filter(value > 0, !is.na(mng)) %>% 
    mutate(name = str_remove(name, '.ageclass1.0'))
  write_parquet(youngs_df, youngs_df_fname)
}
head(youngs_df)
(plt001_yng <- youngs_df %>% 
    filter(name != 'sasa_spp') %>% 
    ggplot(aes(x = name, y = value)) +
    geom_boxplot() +
    geom_jitter(color = 'grey', alpha = .3, size = .1) +
    xlab('AGB (g-biomass/m2) at 2016') + ylab('N of grids') +
    # facet_wrap(~name, scales = 'free') +
    theme_Publication(base_size = 20))
ggsave(file.path(result_dir, 'fig000_saplings_AGB.pdf'), plt001_yng, width = 16, height = 9)


# 2. Preprocessing ----------
# Select total biomass in T1 and T3
total_tree_wide <- resilience_biom_df %>%
  left_join(climate_lookup, by = 'climate') %>% 
  dplyr::select(-climate) %>% 
  rename(climate = climatelabel) %>% 
  mutate(t1 = t1 - byspp_t1_sasa_spp,
         t3 = t3 - byspp_t3_sasa_spp,
         res_biom = t3/t1) %>% 
  dplyr::select(-starts_with("byspp_"), -starts_with('young'), -starts_with('u25_'), -starts_with('u50_'), -starts_with('o50_'))
total_tree_long <- total_tree_wide %>% 
  gather(key = time, val = agb, t1, t3)
# Select species biomass in T1 and T3
spp_wide_buf <- resilience_biom_df %>%
  left_join(climate_lookup, by = 'climate') %>% 
  dplyr::select(-climate) %>% 
  rename(climate = climatelabel) %>% 
  mutate(t1 = t1 - byspp_t1_sasa_spp,
         t3 = t3 - byspp_t3_sasa_spp,
         res_biom = t3/t1) %>% 
  dplyr::select(-starts_with('young'), -starts_with('u25_'), -starts_with('u50_'), -starts_with('o50_'))
# Identify dominant species
spp_wide <- spp_wide_buf %>% 
  # t1 is all tree spp, t3 can be a grass species
  mutate(t1_dom = apply(X = spp_wide_buf[paste0('byspp_t1_', sppnames[1:(length(sppnames) - 1)])],
                        MARGIN = 1, which.max),
         t3_dom = apply(X = spp_wide_buf[paste0('byspp_t3_', sppnames[1:(length(sppnames))])],
                        MARGIN = 1, which.max),
         t1_dom = sppnames[t1_dom], t3_dom = sppnames[t3_dom])
# AGB by age classes
ageclass_long <- resilience_biom_df %>% 
  mutate(t1 = t1 - byspp_t1_sasa_spp,
         t3 = t3 - byspp_t3_sasa_spp,
         res_biom = t3/t1) %>% 
  dplyr::select(-starts_with('young'), -starts_with('byspp_')) %>% 
  pivot_longer(names_to = 'ac_spp', values_to = 'agb', c(starts_with('u25_'), starts_with('u50_'), starts_with('o50_'))) %>% 
  separate(col = ac_spp, into = c('ageclass', 'sppname'), sep = 3) %>% 
  mutate(sppname = str_sub(sppname, 2L, 9L))
# Count sites
total_tree_wide %>% 
  filter(year == 2013) %>% 
  group_by(wtcase, management, climate, region, nrep) %>% 
  summarise(n = n())
# Count species
spp_wide %>% 
  filter(year == 2013, management == 'cl', climate == 'Current') %>% 
  group_by(wtcase, management, climate, region, nrep, t1_dom) %>% 
  summarise(n = n())
spp_wide %>% 
  filter(year == 2013, management == 'cl', climate == 'Current') %>% 
  group_by(wtcase, management, climate, region, trt, nrep, t1_dom) %>% 
  summarise(n = n())



# Figure 4. Recovery of AGB at 1753 grid cells -------
## Plantation forest -----
agb_grid_ens_df <- total_tree_wide %>%
  filter(year %in% c(2051)) %>%
  # Grid mean by nreps
  group_by(gridID, mng_name, trt, climate) %>%
  summarise(t1_mean = mean(t1/1000, na.rm = T), t1_se = sd(t1/1000, na.rm = T)/sqrt(nreps),
            t3_mean = mean(t3/1000, na.rm = T), t3_se = sd(t3/1000, na.rm = T)/sqrt(nreps),
            n = n()) %>% ungroup()
agb_trt_ens_df <- total_tree_wide %>%
  filter(year %in% c(2051)) %>%
  # Mean between damaged sites
  group_by(mng_name, trt, climate, year, nrep) %>%
  summarise(t1_spmean = mean(t1/1000, na.rm = T), t1_spsd = sd(t1/1000, na.rm = T),
            t3_spmean = mean(t3/1000, na.rm = T), t3_spsd = sd(t3/1000, na.rm = T),
            n = n()) %>% ungroup() %>% 
  # Mean between 5 replicates
  group_by(mng_name, trt, climate) %>% 
  summarise(t1_ensmean = mean(t1_spmean, na.rm = T), t1_enssd = sd(t1_spmean, na.rm = T),
            t3_ensmean = mean(t3_spmean, na.rm = T), t3_enssd = sd(t3_spmean, na.rm = T),
            n = n()) %>% 
  mutate(t1_ensse = t1_enssd / sqrt(n),
         t3_ensse = t3_enssd / sqrt(n))
write_csv(agb_grid_ens_df, file.path(result_dir, 'csv101_agb_bygrid_overall.csv'))
write_csv(agb_trt_ens_df, file.path(result_dir, 'csv102_agb_bytrt_overall.csv'))


### a. Overall ------
plt101_overall_bytrt <- ggplot() +
  geom_abline(slope = 1, intercept = 0) + 
  # Grid mean and standard-error between ensamble members
  geom_point(data = agb_grid_ens_df,
             aes(x = t1_mean, y = t3_mean, color = climate), size = .001, alpha = .1) +
  # Spatial mean and standard-error between ensamble members
  geom_point(data = agb_trt_ens_df, 
             aes(x = t1_ensmean, y = t3_ensmean, color = climate), size = 2) +
  geom_errorbarh(data = agb_trt_ens_df, 
                 aes(y = t3_ensmean, xmin = t1_ensmean - t1_ensse, xmax = t1_ensmean + t1_ensse, color = climate), height = 0) +
  geom_errorbar(data = agb_trt_ens_df, 
                aes(x = t3_ensmean, ymin = t1_ensmean - t1_ensse, ymax = t1_ensmean + t1_ensse, color = climate), width = 0) +
  scale_color_calc(guide = guide_legend(nrow = 2)) + scale_fill_calc(guide = guide_legend(nrow = 2)) + 
  xlim(0, 35) + ylim(0, 35) + coord_equal(expand = c(0.05, 0.05)) +
  labs(title = 'AGB recovery in 2051', 
       x = expression("AGB (kg-dry weight " ~ m^-2 ~ ") in 2015"),
       y = expression("AGB (kg-dry weight " ~ m^-2 ~ ") in 2051"),
       color = 'Climate scenarios:',
       caption = 'Colored plots and errorbars show mean and standard error.
       The cl and slpl denotes control and salvage logging and plantation cases, respectively.
       The slpl1 plants the same species with damaged trees and slpl2 plants trees which can adapt to climate change.') +
  facet_grid(trt~mng_name) + theme_Publication() +
  theme(legend.position = "bottom")
ggsave(file.path(result_dir, 'Figure4_agb_bytrt_overall.pdf'), plt101_overall_bytrt, width = 16, height = 10)



## Figure 5 Mean AGB recovery at damaged sites ----
# mean AGB of current natural forests and plantation forests
t1agb_natural <- total_tree_wide %>% 
  filter(mng_name == 'CL', trt == 'Natural Forest', climate == 'Current', nrep == 0) %>% 
  pull(t1)/1000
t1agb_plantation <- total_tree_wide %>% 
  filter(mng_name == 'CL', trt == 'Plantation Forest', climate == 'Current', nrep == 0) %>% 
  pull(t1)/1000
# summary plot
data121_summary_plt <- total_tree_wide %>% 
  filter(year %in% c(2051), trt == 'Plantation Forest') %>% 
  mutate(mng_case = case_when(management == 'cl' ~ 'CL', 
                              management == 'slpl1' ~ 'SL1',
                              management == 'slpl2' ~ 'SL2')) %>% 
  # Spatial means by ensamble members
  group_by(t1_saplings_agb_cut, t1_sasa_cut, mng_case, climate, nrep) %>% 
  summarise(t3_mean = mean(t3/1000, na.rm = TRUE), t3_sd = sd(t3/1000, na.rm = TRUE), n = n()) %>% ungroup() %>% 
  # Spatial mean between ensamble members
  group_by(t1_saplings_agb_cut, t1_sasa_cut, mng_case, climate) %>% 
  summarise(min = min(t3_mean), median = median(t3_mean),
            mean = mean(t3_mean), sd = sd(t3_mean), max = max(t3_mean), n = n()) %>% ungroup() %>% 
  mutate(se = sd/sqrt(n))
write_csv(data121_summary, file.path(result_dir, 'csv111_agb_plantation_sapling_agb_2051.csv'))
(plt121_sap_tree_summary_plt_2051 <- data121_summary_plt %>%
    ggplot(aes(x = t1_saplings_agb_cut, y = mean, ymin = mean - se, ymax = mean + se,
               fill = climate, group = interaction(climate, t1_sasa_cut, mng_case))) +
    geom_hline(yintercept = mean(t1agb_plantation), linetype = 'dotted', size = .5) +
    geom_bar(position = position_dodge(width = .8), width = .7, stat = 'identity') +
    geom_errorbar(width = .3, position = position_dodge(width = .8)) +
    scale_color_calc() + scale_fill_calc() +
    labs(title = 'Tree AGB recovery in Plantation Forests',
         subtitle = 'Saplings vs. Tree AGB recovery',
         fill = 'Climate scenarios:',
         x = expression("AGB (kg-dry weight " ~ m^-2 ~ ") of saplings in 2015"),
         y = expression("AGB (kg-dry weight " ~ m^-2 ~ ") of trees in 2051"),
         caption = 'Plots and errorbars show mean and standard error between five replicates.') +
    facet_grid(rows = vars(t1_sasa_cut), cols = vars(mng_case)) +
    theme_Publication(base_size = 20) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = 'bottom'))
data121_summary_nat <- total_tree_wide %>% 
  filter(year %in% c(2051), management != 'slpl2', trt == 'Natural Forest') %>% 
  mutate(mng_case = if_else(management == 'cl', 'CL', 'SL')) %>% 
  # Spatial means by ensamble members
  group_by(t1_saplings_agb_cut, t1_sasa_cut, mng_case, climate, nrep) %>% 
  summarise(t3_mean = mean(t3/1000, na.rm = TRUE), t3_sd = sd(t3/1000, na.rm = TRUE), n = n()) %>% ungroup() %>% 
  # Spatial mean between ensamble members
  group_by(t1_saplings_agb_cut, t1_sasa_cut, mng_case, climate) %>% 
  summarise(min = min(t3_mean), median = median(t3_mean),
            mean = mean(t3_mean), sd = sd(t3_mean), max = max(t3_mean), n = n()) %>% ungroup() %>% 
  mutate(se = sd/sqrt(n))
write_csv(data121_summary, file.path(result_dir, 'csv111_agb_natural_sapling_agb_2051.csv'))
(plt121_sap_tree_summary_nat_2051 <- data121_summary_nat %>%
    ggplot(aes(x = t1_saplings_agb_cut, y = mean, ymin = mean - se, ymax = mean + se,
               fill = climate, group = interaction(climate, t1_sasa_cut, mng_case))) +
    geom_hline(yintercept = mean(t1agb_natural), linetype = 'dotted', size = .5) +
    geom_bar(position = position_dodge(width = .8), width = .7, stat = 'identity') +
    geom_errorbar(width = .3, position = position_dodge(width = .8)) +
    scale_color_calc() + scale_fill_calc() +
    # scale_x_discrete(expand = c(0, 0)) +
    labs(title = 'Tree AGB recovery in CL case',
         subtitle = 'Saplings vs. Tree AGB recovery',
         fill = 'Climate scenarios:',
         x = expression("AGB (kg-dry weight " ~ m^-2 ~ ") of saplings in 2015"),
         y = expression("AGB (kg-dry weight " ~ m^-2 ~ ") of trees in 2051"),
         caption = 'Plots and errorbars show mean and standard error between five replicates.') +
    facet_grid(rows = vars(t1_sasa_cut), cols = vars(mng_case)) +
    theme_Publication(base_size = 20) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          # panel.grid.major = element_line(color = 'grey'),
          legend.position = 'bottom'))
plt121 <- plt121_sap_tree_summary_plt_2051 / plt121_sap_tree_summary_nat_2051
ggsave(file.path(result_dir, 'Figure5_agb_sapling_agb_2051.pdf'), plt121, width = 14, height = 20)


spp_wide %>% 
  filter(year %in% c(2051), management != 'slpl2') %>% 
  mutate(mng_case = if_else(management == 'cl', 'CL', 'SL')) %>% 
  group_by(t1_saplings_agb_cut, t1_sasa_cut, mng_case, climate, t3_dom) %>% 
  summarise(min = min(t3 / 1000), median = median(t3 / 1000),
            mean = mean(t3 / 1000), sd = sd(t3 / 1000), max = max(t3 / 1000),
            n = n()) %>% 
  write_csv(file.path(result_dir, 'csv112_agb_cl_sapling_agb_spp_2051.csv'))

# N of sasa-saplings category
total_tree_wide %>% 
  filter(year %in% c(2051), management != 'slpl2') %>% 
  mutate(mng_case = if_else(management == 'cl', 'CL', 'SL')) %>% 
  # Spatial means by ensamble members
  group_by(t1_saplings_agb_cut, t1_sasa_cut, mng_case, climate, nrep, trt) %>% 
  summarise(t3_mean = mean(t3/1000, na.rm = TRUE), t3_sd = sd(t3/1000, na.rm = TRUE), n = n()) %>% 
  filter(nrep == 0) %>% 
  group_by(trt, t1_sasa_cut, mng_case, climate) %>% 
  summarise(n = sum(n)) %>% 
  filter(climate == 'Current', mng_case == 'CL')


# Figure 6 Dominant Species change ---------
spp_ensemble <- spp_wide %>% 
  filter(year %in% c(2015, 2021, 2031, 2041, 2051, 2071, 2099)) %>%
  mutate(mng_case = case_when(management == 'cl' ~ 'CL',
                              management == 'slpl1' ~ 'SL1',
                              management == 'slpl2' ~ 'SL2')) %>% 
  select(year, gridID, climate, trt, mng_case, t3_dom, nrep) %>% 
  # Check agreement
  pivot_wider(names_from = nrep, values_from = t3_dom) %>% 
  # Identify forest or NON-forest?
  mutate(t3_isgrass = apply(.[, c('0', '1', '2', '3', '4')] == 'sasa_spp', MARGIN = 1, sum),
         t3_dom_ens = apply(.[, c('0', '1', '2', '3', '4')], MARGIN = 1, 
                            FUN = function(x) names(table(x))[which.max(table(x))]),
         t3_dom_agl = apply(.[, c('0', '1', '2', '3', '4')], MARGIN = 1, 
                            FUN = function(x) as.vector(table(x))[which.max(table(x))]/5)) %>% 
  mutate(t3_dom_ens = case_when(t3_dom_agl >= 0.5 ~ t3_dom_ens, 
                                t3_dom_agl < 0.5 & t3_dom_ens == 'sasa_spp' ~ 'uncertain_sasa',
                                t3_dom_agl < 0.5 & t3_dom_ens != 'sasa_spp' ~ 'uncertain_tree')) %>% 
  mutate(t3_dom_ens = factor(t3_dom_ens, levels = c(spplevels, 'uncertain_sasa', 'uncertain_tree'), exclude = F))
count_all <- spp_ensemble %>% 
  to_lodes_form(axes = 1)
fig202_cl_all <- ggplot(count_all,
                        aes(x = stratum, stratum = t3_dom_ens, alluvium = gridID,
                            fill = t3_dom_ens, label = t3_dom_ens)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback") +
  scale_fill_manual(values = c(scales::brewer_pal(palette = 'Paired')(8), 'grey', 'black'), drop = F) +
  geom_stratum() +
  labs(title = 'Dominant Species in CL case',
       x = 'Year', y = 'No. of grids',
       color = 'Species', fill = 'Species') +
  facet_grid(climate~mng_case) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom')
ggsave(file.path(result_dir, 'Figure6_domspp_all.pdf'), fig202_cl_all, width = 14, height = 14)



# Table S11 cross table from 2013 to 2100 -------
for (m in unique(spp_ensemble$mng_case)) {
  for (c in unique(spp_ensemble$climate)) {
    cat('\n', m, c, '\n')
    spp_ensemble_wider <- spp_ensemble %>% 
      filter(mng_case == m, climate == c) %>% 
      select(year, gridID, climate, trt, mng_case, t3_dom_ens) %>% 
      pivot_wider(names_from = year, values_from = t3_dom_ens)
    buf_table <- table(spp_ensemble_wider$`2015`, spp_ensemble_wider$`2099`)
    buf_mat <- matrix(buf_table, nrow(buf_table), ncol(buf_table))
    colnames(buf_mat) <- colnames(buf_table)
    row.names(buf_mat) <- row.names(buf_table)
    write.csv(as.data.frame(buf_mat), file.path(result_dir, paste0('crosstable_', m, '_', c, '.csv')), row.names = T)
    # Check uncertainty
    cat(sum(buf_mat[, c('uncertain_sasa', 'uncertain_tree')]) / sum(buf_mat) * 100,'\n')
  }
}


### Figure S10 --------
agb_spp_plant_df <- spp_wide %>%
  # filter(climate %in% c('GFDL-CM3_RCP8.5', 'MIROC5_RCP2.6', 'noCC')) %>%
  pivot_longer(cols = starts_with('byspp_t3'), names_to = 'sppname', values_to = 'agb') %>%
  mutate(sppname = str_remove(sppname, 'byspp_t3_')) %>%
  mutate(sppname = fct_relevel(sppname, spplevels)) %>%
  # Mean by nreps
  group_by(gridID, year, management, climate, sppname) %>%
  summarise(agb = mean(agb)) %>%
  ungroup() %>% 
  # Mean by cases
  group_by(year, management, climate, sppname) %>%
  summarise(agb_sum = sum(agb),
            agb_mean = mean(agb),
            agb_sd = sd(agb)) %>%
  ungroup()
treeb_prior <- agb_spp_plant_df %>% 
  filter(year == 2013, management == 'cl', climate == 'Current', sppname != 'sasa_spp') %>% 
  pull(agb_sum) %>% sum() / 10^6
(plt202_agb <- ggplot(agb_spp_plant_df,
                      aes(x = factor(year), y = agb_sum/10^6, fill = sppname)) +
    geom_bar(stat = 'identity') +
    geom_hline(yintercept = treeb_prior, linetype = 'dashed') +
    scale_fill_brewer(palette = 'Paired', drop = F) +
    labs(title = 'Mean AGB among damaged sites by species',
         x = 'Year', y = 'Aboveground Biomass (Mg/m2)', fill = 'Species') +
    facet_grid(climate~management) +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
ggsave(file.path(result_dir, 'FigureS10_agb.pdf'), plt202_agb, width = 16, height = 16)




# Figure S9 Read establishment log file --------
ReadEstLog <- function(fname) {
  buf <- str_split(fname, pattern = '/')[[1]]
  foldernames <- buf[(length(buf) - 3):(length(buf) - 1)]
  management <- str_remove(foldernames[1], pattern = 'BaU_')
  climate <- foldernames[2]
  iter <- as.numeric(str_remove(foldernames[3], pattern = 'reg1_iter'))
  # Read file
  est_df <- read_csv(fname, col_types = cols()) %>% 
    pivot_longer(cols = c(-'Time', -'Species', -'ClimateRegion', -'NumberSitesChecked')) %>% 
    mutate(Species = fct_relevel(Species, spplevels),
           climate = climate, management = management, iter = iter)
  return(est_df)
}

est_df_fname <- file.path(result_dir, 'data_101_estlog_all.parquet')
if (file.exists(est_df_fname)) {
  est_df <- read_parquet(est_df_fname)
} else {
  est_fnames <- list.files(root_dir, pattern = 'NECN-prob-establish-log.csv', full.names = TRUE, recursive = TRUE)
  est_df <- map(est_fnames, ReadEstLog) %>% 
    reduce(bind_rows) %>% 
    mutate(management = case_when(management == 'cl' ~ 'CL',
                                  management == 'slpl1' ~ 'SL1',
                                  management == 'slpl2' ~ 'SL2'),
           climate = case_when(climate == "CSIRO-Mk3-6-0_RCP8.5" ~ "RCP8.5 CSIRO",
                               climate == "GFDL-CM3_RCP8.5" ~ 'RCP8.5 GFDL',
                               climate == "HadGEM2-ES_RCP8.5" ~ 'RCP8.5 HadGEM2-ES',
                               climate == "MIROC5_RCP2.6" ~ 'RCP2.6 MIROC5',
                               climate == "noCC" ~ 'Current')) %>% 
    mutate(Species = fct_relevel(Species, spplevels),
           name = fct_relevel(name, 'AvgProbEst', 'AvgTempMult', 'AvgMinJanTempMult', 'AvgSoilMoistureMult'))
  write_parquet(est_df, est_df_fname)
}
# Ensemble mean
est_ens_df <- est_df %>% 
  group_by(climate, management, Time, Species, ClimateRegion, name) %>% 
  summarise(value_mean = mean(value, na.rm = TRUE),
            value_sd = sd(value, na.rm = TRUE),
            value_se = value_sd / sqrt(nreps),
            nsite_mean = mean(NumberSitesChecked, na.rm = TRUE),
            nsite_sd = sd(NumberSitesChecked, na.rm = TRUE),
            nsite_se = nsite_sd / sqrt(nreps)) %>% 
  ungroup()
skimr::skim(est_ens_df)


# climate condition and its prob are same among nreps
est_plt1 <- est_df %>% 
  group_by(Species, ClimateRegion, management, climate, name) %>% 
  mutate(value_rollmean = rollmean(value, 5, na.pad = T)) %>% 
  ungroup() %>% 
  filter(ClimateRegion == 'eco001', management == 'CL') %>% 
  ggplot(aes(x = Time + 2014 - 1, color = Species)) +
  geom_line(aes(y = value), size = .5, alpha = .7) +
  geom_line(aes(y = value_rollmean), size = 1.5) +
  scale_color_brewer(palette = 'Paired', drop = F) +
  facet_grid(name~climate) +
  labs(title = 'Establishment probs in Eco 001',
       x = 'Year', y = 'Probability (-)', color = 'Species') +
  theme_clean(base_size = 20) +
  theme(panel.background = element_rect(color = 'black', fill = NA),
        legend.background = element_rect(color = NA, fill = NA),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13),
        legend.position = 'bottom')
ggsave(file.path(result_dir, 'FigureS9-1_est_eco001_cl.pdf'), est_plt1, width = 16, height = 12)

est_plt2 <- est_df %>% 
  group_by(Species, ClimateRegion, management, climate, name) %>% 
  mutate(value_rollmean = rollmean(value, 5, na.pad = T)) %>% 
  ungroup() %>% 
  filter(ClimateRegion == 'eco002', management == 'CL') %>% 
  ggplot(aes(x = Time + 2014 - 1, color = Species)) +
  geom_line(aes(y = value), size = .5, alpha = .7) +
  geom_line(aes(y = value_rollmean), size = 1.5) +
  scale_color_brewer(palette = 'Paired', drop = F) +
  facet_grid(name~climate) +
  labs(title = 'Establishment probs in Eco 002',
       x = 'Year', y = 'Probability (-)', color = 'Species') +
  theme_clean(base_size = 20) +
  theme(panel.background = element_rect(color = 'black', fill = NA),
        legend.background = element_rect(color = NA, fill = NA),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13),
        legend.position = 'bottom')
ggsave(file.path(result_dir, 'FigureS9-2_est_eco002_cl.pdf'), est_plt2, width = 16, height = 12)

est_plt3 <- est_df %>% 
  group_by(Species, ClimateRegion, management, climate, name) %>% 
  mutate(value_rollmean = rollmean(value, 5, na.pad = T)) %>% 
  ungroup() %>% 
  filter(ClimateRegion == 'eco003', management == 'CL') %>% 
  ggplot(aes(x = Time + 2014 - 1, color = Species)) +
  geom_line(aes(y = value), size = .5, alpha = .7) +
  geom_line(aes(y = value_rollmean), size = 1.5) +
  scale_color_brewer(palette = 'Paired', drop = F) +
  facet_grid(name~climate) +
  labs(title = 'Establishment probs in Eco 003',
       x = 'Year', y = 'Probability (-)', color = 'Species') +
  theme_clean(base_size = 20) +
  theme(panel.background = element_rect(color = 'black', fill = NA),
        legend.background = element_rect(color = NA, fill = NA),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13),
        legend.position = 'bottom')
ggsave(file.path(result_dir, 'FigureS9-3_est_eco003_cl.pdf'), est_plt3, width = 16, height = 12)
