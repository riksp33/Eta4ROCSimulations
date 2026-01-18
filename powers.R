# Setup needed to install the Eta4ROC package from github
library(remotes)
remotes::install_github("riksp33/Eta4ROC")
library(Eta4ROC)

# Setup needed to generate the files
library(jsonlite)
library(here)
library(parallel)


get_power = function(file_name,
                     use_box_cox_in_parametric = FALSE,
                     use_box_cox_in_kernel = FALSE,
                     param_adjuster_function = function(x, sd, auc, tol = 1e-3){x},
                     case = c("gaussian", "lognormal", "gamma"),
                     controls_params = list(param1 = 1, param2 = 1), # params depend on case
                     cases_params = list(param1 = 1.1, param2 = 1.1),
                     MC = 1000,
                     BootstrapSize = 500,
                     alpha = 0.05
){
  set.seed(1)

  AUCs = c(0.6, 0.75, 0.9)
  ns = c(20, 50, 100)
  t0s = c(0.2, 0.4, 0.8, 1)

  # Sample generator
  if(case == "gaussian") {
    sample_distribution = function(n, mu, sigma) rnorm(n, mu, sigma)
  } else if(case == "lognormal") {
    sample_distribution = function(n, meanlog, sdlog) rlnorm(n, meanlog, sdlog)
  } else if(case == "gamma") {
    sample_distribution = function(n, shape, rate) rgamma(n, shape = shape, rate = rate)
  } else {
    stop("Parameter 'case' must be 'gaussian', 'lognormal' or 'gamma'")
  }

  results_summary = list(header = paste("Simulation results for", file_name))

  for(auc in AUCs){
    auc_list = list()

    # obtain missing param (e.g. mu for cases) so that theoretical AUC matches
    test_sample = sample_distribution(1e5, controls_params$param1, controls_params$param2)
    missing_param = param_adjuster_function(test_sample, cases_params$param2, auc, tol = 0.001)
    auc_list$missing_param = missing_param

    for(t0 in t0s){
      t0_list = list()

      true_eta = analytical_eta(controls_params$param1,
                                controls_params$param2,
                                missing_param,
                                cases_params$param2,
                                case = case,
                                t0 = t0)

      t0_list$true_eta = true_eta

      for(n in ns){
        # containers
        parametric_base = numeric(MC)
        kernel_hscv_base = numeric(MC)
        kernel_opt_base = numeric(MC)
        kernel_iqr_base = numeric(MC)
        auc_base = numeric(MC)

        parametric_boot = matrix(NA_real_, nrow = BootstrapSize, ncol = MC)
        kernel_hscv_boot = matrix(NA_real_, nrow = BootstrapSize, ncol = MC)
        kernel_opt_boot = matrix(NA_real_, nrow = BootstrapSize, ncol = MC)
        kernel_iqr_boot = matrix(NA_real_, nrow = BootstrapSize, ncol = MC)
        auc_boot = matrix(NA_real_, nrow = BootstrapSize, ncol = MC)

        parametric_p = numeric(MC)
        kernel_hscv_p = numeric(MC)
        kernel_opt_p = numeric(MC)
        kernel_iqr_p = numeric(MC)
        auc_p = numeric(MC)

        for(mc_it in seq_len(MC)){
          controls = sample_distribution(n, controls_params$param1, controls_params$param2)
          cases = sample_distribution(n, missing_param, cases_params$param2)

          # base estimates on original sample
          parametric_base[mc_it] = parametric_eta(controls, cases, t0, box_cox = use_box_cox_in_parametric)
          kernel_hscv_base[mc_it] = kernel_eta(controls, cases, "hscv", t0, box_cox = use_box_cox_in_kernel)
          kernel_opt_base[mc_it] = kernel_eta(controls, cases, "optimal", t0, box_cox = use_box_cox_in_kernel)
          kernel_iqr_base[mc_it] = kernel_eta(controls, cases, "iqr", t0, box_cox = use_box_cox_in_kernel)
          auc_base[mc_it] = max(calculate_auc_normal(cases, controls), calculate_auc_normal(controls, cases))

          # bootstrap
          for(bc_it in seq_len(BootstrapSize)){
            combined_boot = sample(c(controls, cases), replace = TRUE)
            controls_b = combined_boot[1:n]
            cases_b = combined_boot[(n + 1):(2 * n)]

            parametric_boot[bc_it, mc_it] = parametric_eta(controls_b, cases_b, t0, box_cox = use_box_cox_in_parametric)
            kernel_hscv_boot[bc_it, mc_it] = kernel_eta(controls_b, cases_b, "hscv", t0, box_cox = use_box_cox_in_kernel)
            kernel_opt_boot[bc_it, mc_it] = kernel_eta(controls_b, cases_b, "optimal", t0, box_cox = use_box_cox_in_kernel)
            kernel_iqr_boot[bc_it, mc_it] = kernel_eta(controls_b, cases_b, "iqr", t0, box_cox = use_box_cox_in_kernel)
            auc_boot[bc_it, mc_it] = max(calculate_auc_normal(cases_b, controls_b), calculate_auc_normal(controls_b, cases_b))
          }

          # p-values from bootstrap
          parametric_p[mc_it] = mean(parametric_boot[, mc_it] >= parametric_base[mc_it])
          kernel_hscv_p[mc_it] = mean(kernel_hscv_boot[, mc_it] >= kernel_hscv_base[mc_it])
          kernel_opt_p[mc_it] = mean(kernel_opt_boot[, mc_it] >= kernel_opt_base[mc_it])
          kernel_iqr_p[mc_it] = mean(kernel_iqr_boot[, mc_it] >= kernel_iqr_base[mc_it])
          auc_p[mc_it] = mean(auc_boot[, mc_it] >= auc_base[mc_it])
        }

        # summary metrics (power at alpha). bias/rmse removed (computed elsewhere).
        summary_result = list(
          parametric = list(
            power = mean(parametric_p < alpha)
          ),
          kernel_hscv = list(
            power = mean(kernel_hscv_p < alpha)
          ),
          kernel_opt = list(
            power = mean(kernel_opt_p < alpha)
          ),
          kernel_iqr = list(
            power = mean(kernel_iqr_p < alpha)
          ),
          auc = list(
            power = mean(auc_p < alpha)
          )
        )

        # prepare heavy detailed data and save as compressed RDS in results_powers_sim
        dir_path = here("results_powers_sim")
        if(!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
        rds_path = file.path(dir_path, paste0(file_name, "_n", n, "_t0", gsub("\\.", "_", as.character(t0)), "_auc", gsub("\\.", "_", as.character(auc)), ".rds"))

        detailed = list(
          meta = list(file = file_name, n = n, t0 = t0, auc = auc, case = case, true_eta = true_eta),
          base = list(parametric = parametric_base, kernel_hscv = kernel_hscv_base, kernel_opt = kernel_opt_base, kernel_iqr = kernel_iqr_base, auc = auc_base),
          bootstrap = list(parametric = parametric_boot, kernel_hscv = kernel_hscv_boot, kernel_opt = kernel_opt_boot, kernel_iqr = kernel_iqr_boot, auc = auc_boot),
          p_values = list(parametric = parametric_p, kernel_hscv = kernel_hscv_p, kernel_opt = kernel_opt_p, kernel_iqr = kernel_iqr_p, auc = auc_p)
        )

        saveRDS(detailed, file = rds_path, compress = "xz")

        # attach summary and pointer to the detailed file
        summary_result$rds = basename(rds_path)
        t0_list[[paste0("size:", n)]] = summary_result
      }
      auc_list[[paste0("t0:", t0)]] = t0_list
    }
    results_summary[[paste0("AUC:", auc)]] = auc_list
  }

  # write summary JSON
  json = toJSON(results_summary, pretty = TRUE, digits = NA)
  meta_path = file.path(here("results_powers_sim"), paste0(file_name, "_summary.json"))
  if(!dir.exists(here("results_powers_sim"))) dir.create(here("results_powers_sim"), recursive = TRUE)
  write(json, file = meta_path)
  return(list(summary = meta_path, detailed_dir = here("results_powers_sim")))
}



run_parallel_simulations = function(category_configs, max_cores = NULL) {

  available_cores = detectCores()
  
  if (is.null(max_cores)) {
    max_cores = max(1, floor(available_cores * 0.75))
  } else {
    max_cores = min(max_cores, available_cores)
  }
  
  cat(paste("Using", max_cores, "cores out of", available_cores, "available"))
  cat("\n\n\n")
  
  for (category_name in names(category_configs)) {
    cat(paste("\nRunning scenarios for:", category_name, "\n"))
    
    configs = category_configs[[category_name]]
  
    cl = makeCluster(max_cores)
    
    clusterEvalQ(cl, {
      library(remotes)
      if (!requireNamespace("Eta4ROC", quietly = TRUE)) {
        remotes::install_github("riksp33/Eta4ROC")
      }
      library(Eta4ROC)
      library(jsonlite)
      library(here)
    })
    clusterExport(cl, c("get_power"))
    results = parLapply(cl, configs, function(config) {
      do.call(get_power, config)
    })
    
    stopCluster(cl)
    
    cat(paste("Simulations finished for:", category_name, "\n"))
    cat(paste("Output files :", paste(basename(unlist(results)), collapse=", "), "\n"))
  }
  
  cat("\nAll simulations executed successfully.\n")
}

gaussian_configs = list(
  normal_1 = list(
    file_name = "normal_1_box_cox_parametric_and_kernel",
    case = "gaussian",
    use_box_cox_in_parametric = TRUE,
    use_box_cox_in_kernel = TRUE,
    param_adjuster_function = get_mux_bisection,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 1)
  ),
  normal_2 = list(
    file_name = "normal_2_box_cox_parametric_and_kernel",
    case = "gaussian",
    use_box_cox_in_parametric = TRUE,
    use_box_cox_in_kernel = TRUE,
    param_adjuster_function = get_mux_bisection,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 1.4)
  ),
  normal_3 = list(
    file_name = "normal_3_box_cox_parametric_and_kernel",
    case = "gaussian",
    use_box_cox_in_parametric = TRUE,
    use_box_cox_in_kernel = TRUE,
    param_adjuster_function = get_mux_bisection,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 3)
  )
)

lognormal_configs = list(
  lognormal_1 = list(
    file_name = "lognormal_1_box_cox_parametric",
    case = "lognormal",
    param_adjuster_function = get_mux_bisection,
    use_box_cox_in_parametric = TRUE,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 0.5)
  ),
  lognormal_2 = list(
    file_name = "lognormal_2_box_cox_parametric",
    case = "lognormal",
    param_adjuster_function = get_mux_bisection,
    use_box_cox_in_parametric = TRUE,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 3/2)
  ),
  lognormal_3 = list(
    file_name = "lognormal_3_box_cox_parametric",
    case = "lognormal",
    param_adjuster_function = get_mux_bisection,
    use_box_cox_in_parametric = TRUE,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 0.2)
  ),
  lognormal_4 = list(
    file_name = "lognormal_4_box_cox_parametric",
    case = "lognormal",
    param_adjuster_function = get_mux_bisection,
    use_box_cox_in_parametric = TRUE,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 2)
  ),
  lognormal_1_bc = list(
    file_name = "lognormal_1_box_cox_parametric_and_kernel",
    case = "lognormal",
    param_adjuster_function = get_mux_bisection,
    use_box_cox_in_parametric = TRUE,
    use_box_cox_in_kernel = TRUE,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 0.5)
  ),
  lognormal_2_bc = list(
    file_name = "lognormal_2_box_cox_parametric_and_kernel",
    case = "lognormal",
    param_adjuster_function = get_mux_bisection,
    use_box_cox_in_parametric = TRUE,
    use_box_cox_in_kernel = TRUE,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 3/2)
  ),
  lognormal_3_bc = list(
    file_name = "lognormal_3_box_cox_parametric_and_kernel",
    case = "lognormal",
    param_adjuster_function = get_mux_bisection,
    use_box_cox_in_parametric = TRUE,
    use_box_cox_in_kernel = TRUE,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 0.2)
  ),
  lognormal_4_bc = list(
    file_name = "lognormal_4_box_cox_parametric_and_kernel",
    case = "lognormal",
    param_adjuster_function = get_mux_bisection,
    use_box_cox_in_parametric = TRUE,
    use_box_cox_in_kernel = TRUE,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 2)
  )
)

gamma_configs = list(
  gamma_1 = list(
    file_name = "gamma_1_box_cox_parametric",
    case = "gamma",
    param_adjuster_function = get_gamma_rate,
    use_box_cox_in_parametric = TRUE,
    controls_params = list(param1 = 0.5, param2 = 0.5),
    cases_params = list(param1 = 0, param2 = 1)
  ),
  gamma_2 = list(
    file_name = "gamma_2_box_cox_parametric",
    case = "gamma",
    param_adjuster_function = get_gamma_rate,
    use_box_cox_in_parametric = TRUE,
    controls_params = list(param1 = 0.5, param2 = 0.5),
    cases_params = list(param1 = 0, param2 = 4)
  ),
  gamma_3 = list(
    file_name = "gamma_3_box_cox_parametric",
    case = "gamma",
    param_adjuster_function = get_gamma_rate,
    use_box_cox_in_parametric = TRUE,
    controls_params = list(param1 = 0.5, param2 = 0.5),
    cases_params = list(param1 = 0, param2 = 1/8)
  ),
  gamma_1_bc = list(
    file_name = "gamma_1_box_cox_parametric_and_kernel",
    case = "gamma",
    param_adjuster_function = get_gamma_rate,
    use_box_cox_in_parametric = TRUE,
    use_box_cox_in_kernel = TRUE,
    controls_params = list(param1 = 0.5, param2 = 0.5),
    cases_params = list(param1 = 0, param2 = 1)
  ),
  gamma_2_bc = list(
    file_name = "gamma_2_box_cox_parametric_and_kernel",
    case = "gamma",
    param_adjuster_function = get_gamma_rate,
    use_box_cox_in_parametric = TRUE,
    use_box_cox_in_kernel = TRUE,
    controls_params = list(param1 = 0.5, param2 = 0.5),
    cases_params = list(param1 = 0, param2 = 4)
  ),
  gamma_3_bc = list(
    file_name = "gamma_3_box_cox_parametric_and_kernel",
    case = "gamma",
    param_adjuster_function = get_gamma_rate,
    use_box_cox_in_parametric = TRUE,
    use_box_cox_in_kernel = TRUE,
    controls_params = list(param1 = 0.5, param2 = 0.5),
    cases_params = list(param1 = 0, param2 = 1/8)
  )
)


# Scenarios to simulate
all_configs = list(
  gaussian = gaussian_configs
)



run_parallel_simulations(all_configs, max_cores = 3)