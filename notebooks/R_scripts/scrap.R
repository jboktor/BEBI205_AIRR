


wkdir <- "/central/groups/MazmanianLab/joeB/BEBI205_AIRR"
emb_dir <- glue("{wkdir}/data/input/embeddings")
data_dir <- glue("{wkdir}/data")
source(glue("{wkdir}/notebooks/R_scripts/_misc_functions.R"))

bcr_esm <- data.table::fread(
  glue("{emb_dir}/BCR_heavy_esm2_t30_150M.csv"),
  header=TRUE, stringsAsFactors = FALSE
)

bcr_esm %<>%
  mutate(
    temp = strex::str_after_nth(Label, "-", 2),
    sample_id = strex::str_before_last(temp, "_"),
    case_control_other_latest = strex::str_after_last(temp, "_")
  ) %>% 
  select(-c(temp))

bcr_esm %>% 
  pull(sample_id) %>% 
  unique() %>% length

# bcr_esm_numeric <- bcr_esm %>% 
#   group_by(sample_id) %>% 
#   slice_sample(n=100) %>% 
#   ungroup() %>% 
#   column_to_rownames(var = "Label") %>% 
#   select_if(is.numeric)
# function to slim down data

bcr_esm_numeric <- bcr_esm %>% 
  slice_sample(n=25) %>% 
  column_to_rownames(var = "Label") %>% 
  select_if(is.numeric)
bcr_esm_numeric %>% dim



  # # Initiate future.batchtools backend for parallel processing
  # future::plan(
  #   future.batchtools::batchtools_slurm,
  #   template = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
  #   resources = list(
  #     name = glue("{get_time()}_mRMRe-processing"),
  #     memory = "5G",
  #     ncpus = 1,
  #     walltime = 1200
  #   )
  # )





# Here we are focusing on the BCR heavy chains in our analysis

# here is a breakdown of IG isotypes ---- 

# Other things to note about our dataset is that there is very little redundancy IG sequence detection 
# Plot distribution of sequence counts per sample


# Samples provided longitudinally 





# Exploratory Data Analysis



bcr_finetuned_df <- readRDS(
  glue(
    "{data_dir}/interim/airr/mRMRe-receptor-selection/",
    "2023-06-08_mRMRe-100_finetuned_BCR-heavy_esm2_t30_150M.rds"
  )
)
bcr_base_df <- readRDS(
  glue(
    "{data_dir}/interim/airr/mRMRe-receptor-selection/",
    "2023-06-08_mRMRe-100_BCR_heavy_esm2_t30_150M.rds"
  )
)

# BASE BCR DF

pca_base_all <- plot_pca_ggplot(bcr_base_df)
pca_base_sampleids <- bcr_base_df %>%
  calculate_sample_means() %>%
  plot_pca_ggplot()
pca_base_sampleids_interact <- bcr_base_df %>%
  calculate_sample_means() %>%
  plot_pca()

ggsave(
  glue("{wkdir}/figures/2023-06-08_pca_bcr_base.png"),
  pca_base_all,
  width = 10, height = 8
)
ggsave(
  glue("{wkdir}/figures/2023-06-08_pca_bcr_base_sampleid-mean.png"),
  pca_base_sampleids,
  width = 10, height = 8
)
saveRDS(pca_base_sampleids_interact, glue(
  "{data_dir}/interim/airr/",
  "{Sys.Date()}_PCA-UMAP_BCR-Heavy-base-sampleid-meansummary.rds"
))

readRDS(glue(
  "{data_dir}/interim/airr/",
  "2023-06-08_PCA-UMAP_BCR-Heavy-base-sampleid-meansummary.rds"
))




# dimred <- list()
# dimred[["pca_all_bcr_base_df"]] <- plot_pca(bcr_base_df)
# dimred[["pca_sid_bcr_base_df"]] <- bcr_base_df %>%
#   calculate_sample_means() %>%
#   plot_pca()




# dimred[["umap_all_bcr_base_df"]] <- plot_umap(bcr_base_df)
# dimred[["umap_sid_bcr_base_df"]] <- bcr_base_df %>%
#   calculate_sample_means() %>%
#   plot_umap()

# Finetuned BCR DF
dimred[["pca_all_bcr_ft_df"]] <- plot_pca(bcr_finetuned_df)
dimred[["pca_sid_bcr_ft_df"]] <- bcr_finetuned_df %>%
  calculate_sample_means() %>%
  plot_pca()

# dimred[["umap_all_bcr_ft_df"]] <- plot_umap(bcr_finetuned_df)
# dimred[["umap_sid_bcr_ft_df"]] <- bcr_finetuned_df %>%
#   calculate_sample_means() %>%
#   plot_umap()

saveRDS(dimred, glue(
  "{data_dir}/interim/airr/",
  "{Sys.Date()}_PCA-UMAP_BCR-Heavy-base-and-fine.rds"
))

# ### Finetuned ESM2 Models
# ```{r}
# #| fig-cap: PCA embedding of all BCR-heavy CDR3 chains.
# dimred[["pca_all_bcr_ft_df"]]
# ```

# ```{r}
# #| fig-cap: PCA embedding of all samples colored by cohort, samples are represented by their mean embedding value for BCR-heavy CDR3 chains.
# dimred[["pca_sid_bcr_ft_df"]]
# ```










# long_metdata <-
#   readRDS(
#     glue("{data_dir}/interim/metadata/2023-06-06_longitudinal_metadata.rds")
#   )
# rnaseq_inv <- read.csv(
#   file = glue(
#     "{data_dir}/input/metadata/2021_v2-5release_0510/",
#     "rna_sample_inventory.csv"
#   ),
#   stringsAsFactors = F, header = TRUE
# )
# rnaseq_metadata <- rnaseq_inv %>% left_join(long_metdata)













ttest_results %>% 
  filter(p <= 1e-5) %>% 
  arrange(Difference) %>% 
  pull(ESM2_dim)

bcr_esm_df %>% 
  ggplot(aes(D25, color=case_control_other_latest)) +
  geom_density(aes(D25, after_stat(scaled))) +
  # labs(title = gt) +
  geom_rug(alpha = 0.1) +
  scale_color_d3() +
  theme_bw()


test_col = "mds_updrs_part_i_summary_score"
test_df = bcr_esm_df_case
embedding_ids %>%
  purrr::set_names() %>%
  purrr::map(\(test_df, dim) lmer(dim ~ !!test_col, data = test_df))


lmer_fun = function(response) {
  form = reformulate("group", response = response)
  lm(form, data = dat)
}


as.formula(
  glue(
  "{test_col} ~ {embedding_ids[1]} +",
  " (1|sample_id)"
  )
)


library(easystats)
library(lme4)

do_lmer <- function(test_col,
                    test_df,
                    emdedding_names) {
  testres <-
    emdedding_names %>%
    purrr::set_names() %>%
    purrr::map( ~ lmer(as.formula(
      glue("{.} ~ {test_col} +",
           " (1|sample_id)")
    ), na.action = na.omit, data = test_df) #%>%
      # report %>% as_tibble
    ) #%>%
    # bind_rows(.id = "ESM2_dim") %>%
    # mutate(Group = test_col)
  return(testres)
}

# model <- lmer(
#   D1 ~ mds_updrs_part_i_summary_score + (1 | sample_id), 
#   data = bcr_esm_df_case
# )


clinical_metrics <- bcr_esm_df %>% 
  colnames() %>% 
  keep(grepl("summary_score", .))
bcr_esm_df_case <- bcr_esm_df %>% 
  filter(case_control_other_latest == "Case")

tic()
clincal_ttest_res <-
  clinical_metrics %>%
  purrr::map(
    ~ do_lmer(
      test_col = "mds_updrs_part_i_summary_score",
      test_df = bcr_esm_df_case,
      emdedding_names = embedding_ids
    )) %>% bind_rows()
toc()

saveRDS(
  clincal_ttest_res,
  glue("{data_dir}/interim/projection_pursuit/bcr-heavy_t-tests_clinical-metrics.rds")
)








do_ttest(
  test_col = "mds_updrs_part_i_summary_score",
  test_df = bcr_esm_df_case,
  emdedding_names = paste0("D", 0:639)
)



parameters(model) %>% as_tibble %>% View
model %>% report %>% as_tibble %>% View

contrasts_df <- modelbased::estimate_contrasts(model)
contrasts_df <- modelbased::estimate_slopes(model)

modelbased::estimate_expectation(model)
modelbased::estimate_means(model)

random <- estimate_grouplevel(model)
random %>% 
  arrange(Coefficient) %>% 
  plot()

lm(test_df[[.]] ~ bcr_esm_df_case[["mds_updrs_part_i_summary_score"]], na.action = na.omit)

# for (metric in clinical_metrics) {
#   ttest_results <- do_ttest(
#     test_df = bcr_esm_df_case,
#     test_col = "case_control_other_latest",
#     emdedding_names = embedding_ids
#   )
# }





# ttest_results_mdsupdrs <-
#   paste0("D", 0:639) %>% 
#   purrr::set_names() %>% 
#   purrr::map(
#     ~ t.test(bcr_esm_df[[.]] ~ bcr_esm_df$case_control_other_latest) %>%
#       report %>%
#       as_tibble
#   ) %>% 
#   bind_rows(.id = "ESM2_dim")

bcr_esm_df %>% colnames
bcr_esm_df %>% 
  drop_na(mds_updrs_part_i_summary_score) %>% 
  pull(case_control_other_latest) %>% 
  table()
select(case_control_other_latest, mds_updrs_part_i_summary_score)













#______________________________________________________________________________
# Louvain graph clustering algo
library(igraph)
library(geomnet)
library(visNetwork)

# Data Preparation --------------------------------------------------------

#Load dataset
data(lesmis)

#Edges
edges <- as.data.frame(lesmis[1])
colnames(edges) <- c("from", "to", "weight")

#Create graph for the algorithms
g <- graph_from_data_frame(edges, directed = FALSE)


# Community Detection ----------------------------------------------------
# Louvain
lc <- cluster_louvain(g)
membership(lc)
communities(lc)
plot(lc, g)






bcr_esm <- data.table::fread(
  glue("{wkdir}/data/interim/airr/BCR_heavy_esm2_t30_150M.csv"),
  header=TRUE, stringsAsFactors = FALSE
)
bcr_esm %<>%
  mutate(
    temp = strex::str_after_nth(Label, "-", 2),
    sample_id = strex::str_before_last(temp, "_"),
    case_control_other_latest = strex::str_after_last(temp, "_")
  ) %>% 
  select(-c(temp))

bcr_esm %>% 
  pull(sample_id) %>% 
  unique() %>% length

bcr_esm_numeric <- bcr_esm %>% 
  group_by(sample_id) %>% 
  slice_sample(n=100) %>% 
  ungroup() %>% 
  column_to_rownames(var = "Label") %>% 
  select_if(is.numeric)
# bcr_esm_numeric %>% glimpse

# cosine similarity distance
tic()
cossim <- coop::cosine(as.matrix(t(bcr_esm_numeric)))
toc()

dist_df <- cossim %>% as.data.frame()
graph_df <- dist_df %>% 
  rownames_to_column(var = "from") %>% 
  pivot_longer(!from, names_to = "to", values_to = "weight") %>% 
  filter(weight != 1)

max(graph_df$weight)
min(graph_df$weight)

#Create graph for the algorithms
g <- graph_from_data_frame(graph_df, directed = FALSE)

# Louvain Clustering
tic()
lc <- cluster_louvain(g)
toc()
lc
membership(lc)
communities(lc)



# euclidean distance matrix
distmat <- distances::distances(as.matrix(bcr_esm_numeric))
dist_df <- distance_columns(distmat, 1:nrow(bcr_esm_numeric)) %>% 
  as.data.frame()

colnames(dist_df) <- rownames(bcr_esm_numeric)
rownames(dist_df) <- rownames(bcr_esm_numeric)

graph_df <- dist_df %>% 
  rownames_to_column(var = "from") %>% 
  pivot_longer(!from, names_to = "to", values_to = "weight") %>% 
  filter(weight != 0)

#Create graph for the algorithms
g <- graph_from_data_frame(graph_df, directed = FALSE)

# Louvain Clustering
lc <- cluster_louvain(g)
membership(lc)
communities(lc)
names(communities(lc))





















