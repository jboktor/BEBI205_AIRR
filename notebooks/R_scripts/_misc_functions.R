# Miscellanous functions

`%nin%` <- Negate(`%in%`)

get_time <- function(){
  print(format(Sys.time(), "%Y-%m-%d_%H:%M:%S"))
}
chunk_func <- function(x, n) {
  split(x, cut(seq_along(x), n, labels = FALSE))
}
mrmre_features <- function(df) {
  dd <- df %>%
    dplyr::select_if(is.numeric) %>%
    t() %>%
    as.data.frame() %>%
    mRMRe::mRMR.data()
  ensemble <- mRMRe::mRMR.ensemble(
    data = dd,
    target_indices = 1,
    solution_count = 1,
    feature_count = 100
  )
  selected_feats <- unique(mRMRe::solutions(ensemble)$`1`[, 1])
  return(df[selected_feats, ])
}


#' batchtools job submission function
#' walltime is in seconds
#' memory is in MB
submit_batchtools_slurm <- function(stdout_dir, batch_df, run_func,
                                    walltime = 3600, memory = 1000, ncpus = 1) {
  message("\n\nCreating Registry:  ", stdout_dir, "\n")
  breg <- batchtools::makeRegistry(
    file.dir = glue::glue("{wkdir}/.cluster_runs/", stdout_dir),
    seed = 42
  )
  breg$cluster.functions <- batchtools::makeClusterFunctionsSlurm(
    template = glue::glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
    scheduler.latency = 0.05,
    fs.latency = 65
  )
  jobs <- batchMap(
    fun = run_func,
    args = batch_df,
    reg = breg
  )
  batchtools::submitJobs(jobs,
    resources = list(
      walltime = walltime,
      memory = memory,
      ncpus = ncpus,
      max.concurrent.jobs = 9999
    )
  )
}

plot_pca <- function(df) {
  data_df <- df %>% select_if(is.numeric)
  labels_df <- df %>%
    select_if(is.character) %>%
    left_join(rnaseq_metadata)
  pca_results <- prcomp(data_df, scale = TRUE)
  pca_results$rotation <- -1 * pca_results$rotation
  pca_results$x <- -1 * pca_results$x
  percent_var <- pca_results$sdev^2 * 100 / sum(pca_results$sdev^2)
  plot_df <- pca_results$x %>%
    as_tibble() %>%
    bind_cols(labels_df)

  fig_pca <- plot_ly(
    plot_df,
    x = ~PC1,
    y = ~PC2,
    z = ~PC3,
    color = ~ plot_df$case_control_other_latest
  ) %>%
    add_markers() %>%
    layout(
      scene = list(
        xaxis = list(title = glue("{round(percent_var[1], 2)}%")),
        yaxis = list(title = glue("{round(percent_var[2], 2)}%")),
        zaxis = list(title = glue("{round(percent_var[3], 2)}%"))
      )
    )
  return(fig_pca)
}

plot_pca_ggplot <- function(df) {
  data_df <- df %>% select_if(is.numeric)
  labels_df <- df %>%
    select_if(is.character) %>%
    left_join(rnaseq_metadata)
  pca_results <- prcomp(data_df, scale = TRUE)
  pca_results$rotation <- -1 * pca_results$rotation
  pca_results$x <- -1 * pca_results$x
  percent_var <- pca_results$sdev^2 * 100 / sum(pca_results$sdev^2)
  plot_df <- pca_results$x %>%
    as_tibble() %>%
    bind_cols(labels_df)

  fig_pca <- plot_df %>%
    ggplot(aes(x = PC1, y = PC2, color = case_control_other_latest)) +
    geom_point(size = 2) +
    labs(
      x = glue("{round(percent_var[1], 2)}%"),
      y = glue("{round(percent_var[2], 2)}%")
    ) +
    theme_bw() +
    scale_color_d3()

  return(fig_pca)
}

plot_umap <- function(df) {
  data_df <- df %>% select_if(is.numeric)
  labels_df <- df %>%
    select_if(is.character) %>%
    left_join(rnaseq_metadata)
  umap_results = umap(data_df, n_components = 3, random_state = 15)
  results_df <-
    data.frame(umap_results[["layout"]]) %>%
    bind_cols(labels_df)
  fig_umap <- plot_ly(
    results_df,
    x = ~X1,
    y = ~X2,
    z = ~X3,
    color = ~ results_df$case_control_other_latest
    ) %>%
    add_markers() %>%
      layout(scene = list(
        xaxis = list(title = "0"),
        yaxis = list(title = "1"),
        zaxis = list(title = "2")
      ))
  fig_umap
}

calculate_sample_means <- function(df){
  df %>%
  group_by(sample_id) %>%
  summarise_if(is.numeric, c("mean" = mean)) %>%
  ungroup()
}