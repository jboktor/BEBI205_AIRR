

# options(future.globals.maxSize= 3e+9)


# core_metdata <- 
#   readRDS(glue("{data_dir}/interim/metadata/2023-06-06_core-metadata.rds")) 

long_metdata <- 
  readRDS(glue("{data_dir}/interim/metadata/2023-06-06_longitudinal_metadata.rds")) 
rnaseq_inv <- read.csv(
  file = glue(
    "{data_dir}/input/metadata/2021_v2-5release_0510/",
    "rna_sample_inventory.csv"
  ),
  stringsAsFactors = F, header = TRUE
)
rnaseq_metadata <- rnaseq_inv %>% left_join(long_metdata)
bcr_esm <- data.table::fread(
  glue("{wkdir}/data/interim/airr/BCR_heavy_esm2_t30_150M.csv"),
  header=TRUE, stringsAsFactors = FALSE
)

bcr_esm_df <- bcr_esm %>%
  mutate(
    temp = strex::str_after_nth(Label, "-", 2),
    sample_id = strex::str_before_last(temp, "_"),
    case_control_other_latest = strex::str_after_last(temp, "_")
  ) %>% 
  select(-c(temp)) %>% 
  dplyr::rename_at(vars(`0`:`639`), ~ paste0("D", .)) %>% 
  mutate_if(is.character, factor) %>% 
  left_join(rnaseq_metadata) %>% 
  filter(case_control_other_latest != "Other") %>% 
  as_tibble()


do_ttest <- function(test_col,
                     test_df,
                     emdedding_names) {
  testres <-
    emdedding_names %>%
    purrr::set_names() %>%
    purrr::map( ~ t.test(test_df[[.]] ~ test_df[[test_col]], na.action = na.omit) %>%
                 report %>% as_tibble) %>%
    bind_rows(.id = "ESM2_dim") %>%
    mutate(Group = test_col)
  return(testres)
}

embedding_ids <- paste0("D", 0:639)
tic()
ttest_results <- do_ttest(
  test_df = bcr_esm_df,
  test_col = "case_control_other_latest",
  emdedding_names = embedding_ids
)
toc()
saveRDS(
  ttest_results,
  glue("{data_dir}/interim/projection_pursuit/bcr-heavy_t-tests_case_control.rds")
)

ttest_results %>% 
  arrange(Difference) %>% 
  mutate(ESM2_dim = factor(ESM2_dim, levels = ESM2_dim)) %>% 
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_pointrange(
    aes(
      x = ESM2_dim,
      y = Difference,
      ymin = CI_low,
      ymax = CI_high,
      color = -log10(p)
    ),
    alpha = 0.5
  ) +
  theme_bw() +
  scale_color_viridis_c(option = "F") +
  expand_limits(x = 0, y = 0) +
  labs(x = "Ranked Embedding Dimensions",
       y = "Difference in Mean (PD - Control) \nRepertoire Embedding",
       color = expression(-log[10] ~ "(P-value)")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "#EBEBEB"))



ttest_results %>% 
  filter(p <= 1e-5) %>% 
  arrange(Difference) %>% 
  pull(ESM2_dim)


bcr_esm_df %>% 
  ggplot(aes(D413, color=case_control_other_latest)) +
  geom_density(aes(D413, after_stat(scaled))) +
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



# library(lsa)
library(coop)
library(distances)

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





















