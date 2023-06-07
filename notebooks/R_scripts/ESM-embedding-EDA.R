library(tidyverse)
library(tidymodels)
library(baguette)
library(xgboost)
library(LiblineaR)
library(naivebayes)
library(glmnet)
library(magrittr)
library(janitor)
library(glue)
library(data.table)
library(strex)
library(tidylo)
library(parallel)
library(doParallel)
library(discrim)
library(ggsci)
library(patchwork)
library(plotly)
library(vip)
library(discrim)
library(themis)
library(umap) 

wkdir <- "/Users/josephboktor/Documents/analyses/BEBI205_AIRR"
meta_RNASeq <- readRDS(
  glue("{wkdir}/data/interim/metadata/RNAseq-metadata.rds")
)
meta_RNASeq %>% glimpse

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
  select(-c(temp, Label)) %>% 
  mutate_if(is.character, factor)

bcr_esm %>% glimpse
bcr_esm$case_control_other_latest %>% table()
# on average there are about 200-300 heavy BCR chains measured per donor
bcr_esm$sample_id %>% table() %>%
  as.data.frame() %>% 
  ggplot() +
  geom_histogram(aes(x=Freq), bins = 100)



#__________________________________________________________
# Subsampling df ----


set.seed(123)

rand_sampids <- bcr_esm %>% 
  select(sample_id, case_control_other_latest) %>% 
  distinct() %>% 
  group_by(case_control_other_latest) %>% 
  slice_sample(n=10) %>% 
  pull(sample_id)

bcr_esm_model_df <- bcr_esm %>% 
  filter(sample_id %in% rand_sampids) %>% 
  group_by(sample_id) %>%
  slice_sample(n = 20) %>%
  ungroup() %>% 
  mutate_if(is.character, factor)

#_______________________________________________________________________________
# EDA ----
#_______________________________________________________________________________


rand_subsample_group <-
  function(df,
           seed,
           ngroups = 100,
           nsubgroups = 500) {
    set.seed(seed)
    rand_sampids <- df %>%
      select(sample_id, case_control_other_latest) %>%
      distinct() %>%
      group_by(case_control_other_latest) %>%
      slice_sample(n = ngroups) %>%
      pull(sample_id)
    df_out <- df %>%
      filter(sample_id %in% rand_sampids) %>%
      group_by(sample_id) %>%
      slice_sample(n = nsubgroups) %>%
      ungroup()
    return(df_out)
  }

get_summary_stats <- function(df) {
  df %>% 
    group_by(case_control_other_latest) %>%
    summarise_if(is.numeric, median) %>%
    column_to_rownames(var ="case_control_other_latest") %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "embedding_dim")
}


individual_avgs_df <-
  1:50 %>%  # 50 df subsamples
  purrr::map(~ rand_subsample_group(bcr_esm, seed = .)) %>% 
  purrr::map(., get_summary_stats) %>% 
  bind_rows(.id = "subsample_n")

subsample_avgs_df <-
  individual_avgs_df %>%
  mutate(case_control = Case - Control,
         case_other = Case - Other) %>%
  group_by(embedding_dim) %>%
  summarize_if(is.numeric,
               c(
                 "mean" = mean,
                 "median" = median,
                 "sd" = sd
               ))

average_embeddings_long <- subsample_avgs_df %>% 
  pivot_longer(!embedding_dim)

plot_df <- average_embeddings_long %>% 
  filter(grepl(c("case_control"), name) | 
           grepl(c("case_other"), name))
  # mutate(embedding_dim = factor(embedding_dim, levels = rank_case_control)) %>% 


rank_case_control <- subsample_avgs_df %>% 
  arrange(case_control_mean) %>% 
  pull(embedding_dim)

df_ranked_case_control <- subsample_avgs_df %>% 
  arrange(case_control_mean) %>% 
  mutate(embedding_dim = factor(embedding_dim, levels = embedding_dim))
df_ranked_case_other <- subsample_avgs_df %>% 
  arrange(case_other_mean) %>% 
  mutate(embedding_dim = factor(embedding_dim, levels = embedding_dim))


p_subsamp1 <- df_ranked_case_control %>% 
  ggplot() +
  geom_pointrange(
    aes(
      x = embedding_dim,
      y = case_control_mean,
      ymin = case_control_mean-case_control_sd,
      ymax = case_control_mean+case_control_sd
    ),
    shape = 21,
    color = "blue",
    alpha = 0.3
  ) +
  geom_point(
    aes(
      x = embedding_dim,
      y = case_other_mean
    ),
    
    color = "orange",
    alpha =  0.6
  ) +
  theme_bw() +
  labs(x ="Embedding Dimension (Ranked by median Case - median Control)", 
       y = expression(Delta ~ " Repertoire Average by Group")) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "#EBEBEB"))
  

p_subsamp2 <- df_ranked_case_other %>% 
  ggplot() +
  geom_pointrange(
    aes(
      x = embedding_dim,
      y = case_other_mean,
      ymin = case_other_mean-case_other_sd,
      ymax = case_other_mean+case_other_sd
    ),
    color = "orange",
    alpha = 0.3
  ) +
  geom_point(
    aes(
      x = embedding_dim,
      y = case_control_mean,
    ),
    color = "blue",
    alpha = 0.6
  ) +
  theme_bw() +
  labs(x ="Embedding Dimension (Ranked by median Case - median Control)", 
       y = expression(Delta ~ " Repertoire Average by Group")) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "#EBEBEB"))


patch_subsamp <- p_subsamp1 / p_subsamp2 + plot_layout(guides = 'collect')
ggsave(glue("{wkdir}/figures/heavy-BCR-embeddings_subsampled-averages-median-embedding-comparisons.png"),
       patch_subsamp,
       width = 7, height = 7)


#______________________________________________________________________________
# direct analysis

average_embeddings <- bcr_esm %>%
  group_by(case_control_other_latest) %>%
  summarise_if(is.numeric, median) %>%
  column_to_rownames(var ="case_control_other_latest") %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(case_control = Case - Control, case_other = Case - Other) %>% 
  rownames_to_column(var = "embedding_dim")

average_embeddings_long <- average_embeddings %>% 
  pivot_longer(!embedding_dim)

average_embeddings %>% 
  ggplot(aes(case_control, case_other)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() + 
  theme()

rank_case_control <- average_embeddings_long %>% 
  filter(name %in% c("case_control")) %>% 
  arrange(value) %>% 
  pull(embedding_dim)
rank_case_other <- average_embeddings_long %>% 
  filter(name %in% c("case_other")) %>% 
  arrange(value) %>% 
  pull(embedding_dim)

# bootstrap resample this data over different sequences for each 
p1 <- average_embeddings_long %>% 
  filter(name %in% c("case_control", "case_other")) %>% 
  mutate(embedding_dim = factor(embedding_dim, levels = rank_case_control)) %>% 
  ggplot(aes(
    x=embedding_dim, 
    y=value, 
    color = name)) +
  geom_point() +
  theme_bw() +
  expand_limits(x = 0, y = 0) +
  labs(x ="Embedding Dimension (Ranked by median Case - median Control)", 
       y = expression(Delta ~ " Repertoire Average by Group"), 
       color = NULL) +
  scale_color_d3() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "#EBEBEB"))

p2 <- average_embeddings_long %>% 
  filter(name %in% c("case_control", "case_other")) %>% 
  mutate(embedding_dim = factor(embedding_dim, levels = rank_case_other)) %>% 
  ggplot(aes(
    x=embedding_dim, 
    y=value, 
    color = name)) +
  geom_point() +
  theme_bw() +
  expand_limits(x = 0, y = 0) +
  labs(x = "Embedding Dimension (Ranked by median Case - median Other)",  
       y = expression(Delta ~ " Repertoire Average by Group"), 
       color = NULL) +
  scale_color_d3() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "#EBEBEB"))

patch1 <- p1 / p2 + plot_layout(guides = 'collect')
ggsave(glue("{wkdir}/figures/heavy-BCR-embeddings_summary-stat-comparison.png"),
       width = 7, height = 7)


#_______________________________________________________________________________
# ML MODELING SETUP ----
#_______________________________________________________________________________

# group_split <- group_initial_split(
#   bcr_esm_model_df, 
#   group = sample_id,
#   strata = case_control_other_latest
#   )
# vb_train <- training(group_split)
# vb_test <- testing(group_split)

embedding_avg_tbl <- bcr_esm %>%
  group_by(sample_id, case_control_other_latest) %>%
  summarise_if(is.numeric, median) %>%
  as_tibble() %>% 
  dplyr::rename_at(vars(`0`:`639`), ~ paste0("D", .) )

bcr_split <- initial_split(
  embedding_avg_tbl, 
  strata = case_control_other_latest
)
bcr_train <- training(bcr_split)
bcr_test <- testing(bcr_split)
bcr_folds <- vfold_cv(bcr_train, strata = case_control_other_latest)


set.seed(234)
comp_folds <- bootstraps(bcr_test, strata = case_control_other_latest)

#_______________________________________________________________________________
#  MODELS ----
#_______________________________________________________________________________

# building a data recipe
rec_all <- 
  bcr_train %>% 
  recipe(case_control_other_latest ~ .) %>%
  step_rm(sample_id) %>% 
  step_zv(all_numeric(), -all_outcomes())
# rec_all_nrom <- rec_all %>% 
#   step_normalize(all_numeric(), -all_outcomes())

# rec_all %>% prep() %>% bake(new_data = NULL) %>% colnames
# rec_all_nrom %>% prep() %>% bake(new_data = NULL) %>% colnames

svm_spec <-
  svm_linear() %>%
  set_mode("classification") %>%
  set_engine("LiblineaR")
nb_spec <-
  naive_Bayes() %>%
  set_mode("classification") %>%
  set_engine("naivebayes")
mlr_spec <- multinom_reg(penalty = tune(), mixture = tune()) %>%  # lasso
  set_mode("classification") %>%
  set_engine("glmnet")
# trees ------------------
bag_spec <- bag_tree(
  # tree_depth = tune(),
  min_n = tune()
) %>%
  set_mode("classification") %>%
  set_engine("C5.0")
bart_spec <-
  parsnip::bart() %>%
  set_mode("classification") %>%
  set_engine("dbarts")
xgb_spec <- boost_tree(
  trees = 1000,
  sample_size = 0.8,
  min_n = 2,
  tree_depth = tune(),
  loss_reduction = tune(),
  # mtry = tune(),
  learn_rate = tune()
) %>%
  set_mode("classification") %>%
  set_engine("xgboost")

# mlr_spec <- multinom_reg(penalty = 0.1, mixture = 1) %>%  # lasso
#   set_mode("classification") %>%
#   set_engine("glmnet")
# # trees ------------------
# bag_spec <- bag_tree(
#   tree_depth = 5,
#   min_n = 2
# ) %>%
#   set_mode("classification") %>%
#   set_engine("C5.0")
# bart_spec <-
#   parsnip::bart() %>%
#   set_mode("classification") %>%
#   set_engine("dbarts")
# xgb_spec <- boost_tree(
#   trees = 1000,
#   tree_depth = 5,
#   min_n = 2,
#   loss_reduction = 0.2,
#   sample_size = 0.8,
#   mtry = 2,
#   learn_rate = 0.2
# ) %>%
#   set_mode("classification") %>%
#   set_engine("xgboost")


comp_models <-
  workflow_set(
    preproc = list(unnormalized = rec_all,
                   normlized = rec_all_nrom),
    models = list(
      nb = nb_spec,
      svm = svm_spec,
      mlr = mlr_spec,
      bart = bart_spec,
      xgb = xgb_spec,
      bag = bag_spec
    ),
    cross = TRUE
  )
comp_models

set.seed(123)
cl <- parallel::makePSOCKcluster(14)
doParallel::registerDoParallel(cl)
computer_rs <-
  comp_models %>%
  workflow_map(
    "tune_bayes",
    resamples = comp_folds,
    metrics = metric_set(accuracy, sensitivity, specificity)
  )
stopCluster(cl)


rank_results(computer_rs) %>%
  filter(.metric == "accuracy") %>% View

p_panelres <- rank_results(computer_rs) %>%
  separate("wflow_id", into = c("norm", "model")) %>% 
  ggplot() +
  geom_pointrange(aes(x=rank, y=mean, ymax=mean+std_err, ymin=mean-std_err, color = model, shape = norm)) +
  facet_wrap(~.metric)
p_panelres  

ggplotly(p_panelres)

# rank_results(computer_rs) %>%
#   separate("wflow_id", into = c("norm", "model")) %>% 
#   saveRDS(glue("{wkdir}/data/2023-06-01_heavyBCR-tuned-model-panel.rds"))
# 
# computer_rs %>% 
#   saveRDS(glue("{wkdir}/data/2023-06-01_heavyBCR-tuned-model-panel_FULL.rds"))


# Train and evaluate final SVM model

comp_wf <- workflow(rec_all, svm_spec)

comp_fitted <- last_fit(
  comp_wf, 
  bcr_split,
  metrics = metric_set(accuracy, sensitivity, specificity)
)

collect_metrics(comp_fitted)

collect_predictions(comp_fitted) %>% 
  conf_mat(case_control_other_latest, .pred_class) %>% 
  autoplot()

#_______________________________________________________________________________
# LASSO MODEL ----
#_______________________________________________________________________________

embedding_avg_tbl <- bcr_esm %>%
  group_by(sample_id, case_control_other_latest) %>%
  summarise_if(is.numeric, median) %>%
  as_tibble() %>% 
  dplyr::rename_at(vars(`0`:`639`), ~ paste0("D", .) )

bcr_esm$sample_id %>% unique() %>% length

# bcr_split <- initial_split(
#   embedding_avg_tbl, 
#   strata = case_control_other_latest
# )
# bcr_train <- training(bcr_split)
# bcr_test <- testing(bcr_split)
# bcr_folds <- vfold_cv(bcr_train, strata = case_control_other_latest)

bcr_split <- rsample::group_initial_split(
  bcr_esm, 
  group = sample_id,
  strata = case_control_other_latest
)

bcr_train <- training(bcr_split)
bcr_test <- testing(bcr_split)
# bcr_folds <- vfold_cv(bcr_train, strata = case_control_other_latest)
bcr_folds <- rsample::group_vfold_cv(bcr_train, 
                                     strata = case_control_other_latest, 
                                     group = sample_id,
                                     v = 10)

# upsample control and other case samples to match case #'s
bcr_rec <- 
  bcr_train %>% 
  recipe(case_control_other_latest ~ .) %>%
  update_role(sample_id, new_role = "id") %>%
  step_downsample(case_control_other_latest)

# View group numbers
prep(bcr_rec) %>% bake(new_data = NULL) %>% pull(case_control_other_latest) %>% table

multi_spec <- multinom_reg(penalty = tune(), mixture = 1) %>% 
  set_mode("classification") %>% 
  set_engine("glmnet")

# tune hyperparameters using cfv
bcr_wf <- workflow(bcr_rec, multi_spec)
bcr_tune_grid <- grid_regular(penalty(range = c(-5, 0)), levels = 5)
bcr_tune_grid

set.seed(42)
cl <- parallel::makePSOCKcluster(14)
doParallel::registerDoParallel()
bcr_rs <-
  tune_grid(
    bcr_wf,
    bcr_folds,
    grid = bcr_tune_grid
  )
stopCluster(cl)

# evaluate tune params
autoplot(bcr_rs)
show_best(bcr_rs)
show_best(bcr_rs, metric = "accuracy")

# choose and evaluate final model
# select_best(bcr_rs)
final_params <- bcr_rs %>% 
  select_by_one_std_err(metric = "roc_auc", desc(penalty))

final_res <- bcr_wf %>% 
  finalize_workflow(final_params) %>% 
  last_fit(bcr_split)
final_res

collect_metrics(final_res)
collect_predictions(final_res)

collect_predictions(final_res) %>% 
  conf_mat(case_control_other_latest, .pred_class) %>% 
  autoplot()

bcr_test$case_control_other_latest %>% table


#_______________________________________________________________________________
# XGBOOST MODEL ----
#_______________________________________________________________________________

# xgb_spec <- boost_tree(
#   trees = 1000,
#   tree_depth = 5,
#   min_n = 2,
#   loss_reduction = 0.2,
#   sample_size = 0.8,
#   mtry = 2,
#   learn_rate = 0.2
# ) %>%
#   set_mode("classification") %>%
#   set_engine("xgboost")

xgb_spec <- boost_tree(
  trees = 100,
  tree_depth = 5,
  min_n = 3,
  loss_reduction = tune(),
  sample_size = tune(),
  mtry = tune(),
  learn_rate = tune()                          ## step size
) %>%
  set_engine("xgboost") %>%
  set_mode("classification")
xgb_spec

xgb_grid <- grid_latin_hypercube(
  # trees = 100,
  # tree_depth = 5,
  # min_n = 3,
  loss_reduction(),
  sample_size = sample_prop(),
  finalize(mtry(), bcr_train),
  learn_rate(),
  size = 5
)
xgb_grid

xgb_wf <- workflow() %>%
  add_formula(case_control_other_latest ~ `0`:`639`) %>%
  add_model(xgb_spec)
xgb_wf

set.seed(123)
vb_folds <- group_vfold_cv(
  bcr_train,
  group = sample_id,
  strata = case_control_other_latest,
  v = 4)
vb_folds


cl <- parallel::makePSOCKcluster(10)
doParallel::registerDoParallel(cl)
set.seed(234)
xgb_res <- tune_grid(
  xgb_wf,
  resamples = vb_folds,
  grid = xgb_grid,
  control = control_grid(save_pred = TRUE)
)
stopCluster(cl)
xgb_res


collect_metrics(xgb_res)

xgb_res %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  select(mean, mtry:sample_size) %>%
  pivot_longer(mtry:sample_size,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "AUC")

final_xgb_params <- select_best(xgb_res, "roc_auc")
final_xgb_params

final_res <- xgb_wf %>% 
  finalize_workflow(final_xgb_params) %>% 
  last_fit(bcr_split)
final_res

collect_metrics(final_res)
collect_predictions(final_res)

collect_predictions(final_res) %>% 
  conf_mat(case_control_other_latest, .pred_class) %>% 
  autoplot()
# # ROC curve
# collect_predictions(final_res) %>% 
#   roc_curve(case_control_other_latest, .pred_class) %>% 
#   autoplot()

bcr_test$case_control_other_latest %>% table











#_______________________________________________________________________________

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
  select(-c(temp, Label)) %>% 
  mutate_if(is.character, factor)

bcr_esm_df <- bcr_esm %>% 
  group_by(sample_id) %>% 
  summarise_if(is.numeric, c("sum" = sum, "mean" = mean)) %>% 
  ungroup()

bcr_data <- bcr_esm_df %>% select(contains("_mean"))

bcr_labels <- bcr_esm_df %>% 
  select(sample_id) %>% 
  left_join(meta_RNASeq)

pca_results <- prcomp(bcr_data, scale = TRUE)
pca_results$rotation <- -1*pca_results$rotation
pca_results$x <- -1*pca_results$x
percent_var <- pca_results$sdev^2*100 / sum(pca_results$sdev^2)

plot_df <- pca_results$x %>% 
  as_tibble() %>% 
  bind_cols(bcr_labels)
fig_pca <- plot_ly(
  plot_df,
  x = ~ PC1,
  y = ~ PC2,
  z = ~ PC3,
  color = ~ plot_df$case_control_other_latest
) %>%
  add_markers() %>%
  layout(scene = list(
    xaxis = list(title = '0'),
    yaxis = list(title = '1'),
    zaxis = list(title = '2')
  ))

#_______________________________________________________________________________
# UMAP on avgs

umap_results = umap(bcr_data, n_components = 3, random_state = 15) 
results_df <- 
  data.frame(umap_results[["layout"]]) %>% 
  bind_cols(bcr_labels)
fig_umap <- plot_ly(
  results_df,
  x = ~ X1,
  y = ~ X2,
  z = ~ X3,
  color = ~ results_df$case_control_other_latest
) %>%
  add_markers() %>%
  layout(scene = list(
    xaxis = list(title = '0'),
    yaxis = list(title = '1'),
    zaxis = list(title = '2')
  ))
fig_umap


#_______________________________________________________________________________
# Cosine Dissimilarity on avgs

library(lsa)
cossim <- cosine(as.matrix(t(bcr_data)))
pca_mds <- cmdscale(1-as.matrix(cossim), k=3)
pca_mds %>% glimpse

cosine_df <- 
  data.frame(pca_mds) %>% 
  bind_cols(bcr_labels)
fig_cos <- plot_ly(
  cosine_df,
  x = ~ X1,
  y = ~ X2,
  z = ~ X3,
  color = ~ cosine_df$case_control_other_latest
) %>%
  add_markers() %>%
  layout(scene = list(
    xaxis = list(title = '0'),
    yaxis = list(title = '1'),
    zaxis = list(title = '2')
  ))
fig_cos


#_______________________________________________________________________________
#  Calculate Cosine Similarity on entire dataset

# Pass 1-cosine simililarity through a PERMANOVA looping through meta data of interest 

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
  select(-c(Label)) %>% 
  mutate_if(is.character, factor) %>% 
  left_join(meta_RNASeq)

meta_RNASeq %>% glimpse
# meta_RNASeq %>% colnames
# meta_RNASeq$race %>% table

set.seed(123)
rand_sampids <- bcr_esm %>% 
  select(sample_id, case_control_other_latest) %>% 
  distinct() %>% 
  group_by(case_control_other_latest) %>% 
  slice_sample(n=50) %>% 
  pull(sample_id)

bcr_esm_model_df <- bcr_esm %>% 
  filter(sample_id %in% rand_sampids) %>% 
  group_by(sample_id) %>%
  slice_sample(n = 20) %>%
  ungroup() %>% 
  # mutate_if(is.character, factor) %>% 
  dplyr::rename_at(vars(`0`:`639`), ~ paste0("D", .) )

library(vegan)
library(BiodiversityR) 
embedding_mat <- bcr_esm_model_df[,1:640]

dbrda_res <- dbrda(
  formula = embedding_mat ~  visit_month + age_at_baseline + Condition(sample_id),
  data = bcr_esm_model_df,
  distance = "euclidean",
  dfun = vegdist,
  metaMDSdist = F,
  na.action = na.exclude
)
dbrda_res

var_dbrda_res <- anova(dbrda_res, by="terms", permu=200) # test for sign. environ. variables
var_dbrda_res

# scores(dbrda_res)
anova(dbrda_res) # overall test of the significant of the analysis
anova(dbrda_res, by="axis", perm.max=500) # test axes for significance

dbrda_res_summary <- summary(dbrda_res)
dbrda_res_summary %>% View

dbbiplot <- ordiplot(dbrda_res, choices=c(1,2))

 # as.data.frame(dbbiplot$default)

dbbiplot.long <- sites.long(dbbiplot, env.data=bcr_esm_model_df)
head(sites.long1)

dbrda_res$Ybar


data(varespec)
data(varechem)
## Basic Analysis
vare.cap <- capscale(varespec ~ N + P + K + Condition(Al), varechem,
                     dist="bray")
vare.cap
plot(vare.cap)
anova(vare.cap)
## Avoid negative eigenvalues with additive constant
capscale(varespec ~ N + P + K + Condition(Al), varechem,
         dist="bray", add =TRUE)
## Avoid negative eigenvalues by taking square roots of dissimilarities
capscale(varespec ~ N + P + K + Condition(Al), varechem,
         dist = "bray", sqrt.dist= TRUE)
## Principal coordinates analysis with extended dissimilarities
capscale(varespec ~ 1, dist="bray", metaMDS = TRUE)
## dbrda
dbrda(varespec ~ N + P + K + Condition(Al), varechem,
      dist="bray")
## avoid negative eigenvalues also with Jaccard distances
dbrda(varespec ~ N + P + K + Condition(Al), varechem,
      dist="jaccard")





# library(mRMRe)
# dd <- mRMR.data(data = bcr_data)
# ensemble <- mRMR.ensemble(
#   data = dd,
#   target_indices = 7,
#   solution_count = 1,
#   feature_count = 100
# )
# 
# bcr_data_select <- bcr_data[, solutions(ensemble)$`7`[,1]]
# pca_results <- prcomp(bcr_data_select, scale = TRUE)
# pca_results$rotation <- -1*pca_results$rotation
# pca_results$rotation
# 
# pca_results$x <- -1*pca_results$x
# head(pca_results$x)
# #calculate total variance explained by each principal component
# pca_results$sdev^2*100 / sum(pca_results$sdev^2)
# pca_results$x %>% glimpse
# 
# # pca_results$x %>% 
# #   as_tibble() %>% 
# #   bind_cols(bcr_labels) %>% 
# #   ggplot(aes(x=PC1, y=PC2)) +
# #   geom_point(aes(color=case_control_other_latest), alpha = 0.6)
# plot_df <- pca_results$x %>% 
#   as_tibble() %>% 
#   bind_cols(bcr_labels)
# 
# fig_pca <- plot_ly(
#   plot_df,
#   x = ~ PC1,
#   y = ~ PC2,
#   z = ~ PC3,
#   color = ~ plot_df$case_control_other_latest
# ) %>%
#   add_markers() %>%
#   layout(scene = list(
#     xaxis = list(title = '0'),
#     yaxis = list(title = '1'),
#     zaxis = list(title = '2')
#   ))
# 
# fig_pca






