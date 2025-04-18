# Setup needed to install the Eta4ROC package from github
library(remotes)
remotes::install_github("riksp33/Eta4ROC")
library(Eta4ROC)

# Setup needed to generate the files
library(jsonlite)
library(here)
library(parallel)


simulate_bias_rmse = function(file_name,
                              use_box_cox_in_eta = FALSE,
                              use_box_cox_before_eta = FALSE,
                              param_adjuster_function = function(mu){mu},
                              case = c("gaussian", "lognormal", "gamma"),
                              controls_params = list(param1 = 1, param2 = 1), #Params for cases (mean, std) for gaussian & lognormal, (shape, rate) for gamma 
                              cases_params = list(param1 = 1.1, param2 = 1.1) #Params for controls (mean, std) for gaussian & lognormal, (shape, rate) for gamma 
                              ){
  # Set the seed
  set.seed(1)
  
  # Output JSON
  ret_json = list()
  ret_json[["header"]] = paste("Simulation results for ", file_name)

  # Simulation params
  MC = 1
  AUCs = c(0.6, 0.75, 0.9)
  ns = c(20, 50, 100)
  t0s = c(0.2, 0.4, 0.8, 1)
  
  # Sample generator function
  if(case == "gaussian") {
    sample_distribution = function(n, mu, sigma) {rnorm(n, mu, sigma)}
  } else if(case == "lognormal") {
    sample_distribution = function(n, mu, sigma) {rlnorm(n, mu, sigma)}
  } else if(case == "gamma") {
    sample_distribution = function(n, shape, rate) {rgamma(n, shape, rate = rate)}
  } else {
    stop("Non valid input, the param case must be gaussian, lognormal or gamma")
  }
  
  for(auc in AUCs){
    ret_aucs = list()    
    
    test = sample_distribution(1e5, controls_params$param1, controls_params$param2)
    missing_param = param_adjuster_function(test,
                                            cases_params$param2,
                                            auc,
                                            tol = 0.001)
    
    ret_aucs[["missing_param"]] = missing_param
    
    for(t0 in t0s){
      ret_t0 = list()

      # True eta
      true_eta = analytical_eta(controls_params$param1,
                                controls_params$param2,
                                missing_param,
                                cases_params$param2,
                                case = case,
                                t0 = t0)
      
      ret_t0[["true_eta"]] = true_eta
      
      for(n in ns){
        
        parametric = numeric(MC)
        kernel_hscv = numeric(MC)
        kernel_opt = numeric(MC)
        kernel_iqr = numeric(MC)
        
        for(i in 1:MC){
          
          # Generate samples
          controls = sample_distribution(n, controls_params$param1, controls_params$param2)
          cases = sample_distribution(n, missing_param, cases_params$param2)
          
          if(use_box_cox_before_eta){
            transformed = apply_box_cox(controls, cases)
            controls = transformed$transformed_x
            cases = transformed$transformed_y
          }
          
          parametric[i] = parametric_eta(controls, cases, t0, box_cox = use_box_cox_in_eta)
          kernel_hscv[i] = kernel_eta(controls, cases, "hscv", t0, box_cox = use_box_cox_in_eta)
          kernel_opt[i] = kernel_eta(controls, cases, "optimal", t0, box_cox = use_box_cox_in_eta)
          kernel_iqr[i] = kernel_eta(controls, cases, "iqr", t0, box_cox = use_box_cox_in_eta)
        }
        
        result = list(
          parametric = list(
            bias = get_bias(parametric, true_eta),
            rmse = get_rmse(parametric, true_eta)
          ),
          kernel_hscv = list(
            bias = get_bias(kernel_hscv, true_eta),
            rmse = get_rmse(kernel_hscv, true_eta)
          ),
          kernel_opt = list(
            bias = get_bias(kernel_opt, true_eta),
            rmse = get_rmse(kernel_opt, true_eta)
          ),
          kernel_iqr = list(
            bias = get_bias(kernel_iqr, true_eta),
            rmse = get_rmse(kernel_iqr, true_eta)
          )
        )
        ret_t0[[paste("size:", n)]] = result
      }
      ret_aucs[[paste("t0:", t0)]] = ret_t0
    }
    ret_json[[paste("AUC:", auc)]] = ret_aucs
  }
  json = toJSON(ret_json, pretty = TRUE, digits = NA)
  dir_path = here("results")
  full_path = file.path(dir_path, paste0(file_name, ".json"))
  write(json, file = full_path)
  return(full_path)
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
    clusterExport(cl, c("simulate_bias_rmse"))
    results = parLapply(cl, configs, function(config) {
      do.call(simulate_bias_rmse, config)
    })
    
    stopCluster(cl)
    
    cat(paste("Simulations finished for:", category_name, "\n"))
    cat(paste("Output files :", paste(basename(unlist(results)), collapse=", "), "\n"))
  }
  
  cat("\nAll simulations executed successfully.\n")
}

gaussian_configs = list(
  normal_1 = list(
    file_name = "normal_1",
    case = "gaussian",
    param_adjuster_function = get_mux_bisection,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 1)
  ),
  normal_2 = list(
    file_name = "normal_2",
    case = "gaussian",
    param_adjuster_function = get_mux_bisection,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 1.4)
  ),
  normal_3 = list(
    file_name = "normal_3",
    case = "gaussian",
    param_adjuster_function = get_mux_bisection,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 3)
  )
)

lognormal_configs = list(
  lognormal_1 = list(
    file_name = "lognormal_1_box_cox_inside_eta",
    case = "lognormal",
    param_adjuster_function = get_mux_bisection,
    use_box_cox_in_eta = TRUE,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 0.5)
  ),
  lognormal_2 = list(
    file_name = "lognormal_2_box_cox_inside_eta",
    case = "lognormal",
    param_adjuster_function = get_mux_bisection,
    use_box_cox_in_eta = TRUE,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 3/2)
  ),
  lognormal_3 = list(
    file_name = "lognormal_3_box_cox_inside_eta",
    case = "lognormal",
    param_adjuster_function = get_mux_bisection,
    use_box_cox_in_eta = TRUE,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 0.2)
  ),
  lognormal_4 = list(
    file_name = "lognormal_4_box_cox_inside_eta",
    case = "lognormal",
    param_adjuster_function = get_mux_bisection,
    use_box_cox_in_eta = TRUE,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 2)
  ),
  lognormal_1_bc = list(
    file_name = "lognormal_1_box_cox_inside_eta_box_cox_before_eta",
    case = "lognormal",
    param_adjuster_function = get_mux_bisection,
    use_box_cox_before_eta = TRUE,
    use_box_cox_in_eta = TRUE,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 0.5)
  ),
  lognormal_2_bc = list(
    file_name = "lognormal_2_box_cox_inside_eta_box_cox_before_eta",
    case = "lognormal",
    param_adjuster_function = get_mux_bisection,
    use_box_cox_before_eta = TRUE,
    use_box_cox_in_eta = TRUE,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 3/2)
  ),
  lognormal_3_bc = list(
    file_name = "lognormal_3_box_cox_inside_eta_box_cox_before_eta",
    case = "lognormal",
    param_adjuster_function = get_mux_bisection,
    use_box_cox_before_eta = TRUE,
    use_box_cox_in_eta = TRUE,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 0.2)
  ),
  lognormal_4_bc = list(
    file_name = "lognormal_4_box_cox_inside_eta_box_cox_before_eta",
    case = "lognormal",
    param_adjuster_function = get_mux_bisection,
    use_box_cox_before_eta = TRUE,
    use_box_cox_in_eta = TRUE,
    controls_params = list(param1 = 0, param2 = 1),
    cases_params = list(param1 = 0, param2 = 2)
  )
)

gamma_configs = list(
  gamma_1 = list(
    file_name = "gamma_1_box_cox_inside_eta",
    case = "gamma",
    param_adjuster_function = get_gamma_rate,
    use_box_cox_in_eta = TRUE,
    controls_params = list(param1 = 0.5, param2 = 0.5),
    cases_params = list(param1 = 0, param2 = 1)
  ),
  gamma_2 = list(
    file_name = "gamma_2_box_cox_inside_eta",
    case = "gamma",
    param_adjuster_function = get_gamma_rate,
    use_box_cox_in_eta = TRUE,
    controls_params = list(param1 = 0.5, param2 = 0.5),
    cases_params = list(param1 = 0, param2 = 4)
  ),
  gamma_3 = list(
    file_name = "gamma_3_box_cox_inside_eta",
    case = "gamma",
    param_adjuster_function = get_gamma_rate,
    use_box_cox_in_eta = TRUE,
    controls_params = list(param1 = 0.5, param2 = 0.5),
    cases_params = list(param1 = 0, param2 = 1/8)
  ),
  gamma_1_bc = list(
    file_name = "gamma_1_box_cox_inside_eta_box_cox_before_eta",
    case = "gamma",
    param_adjuster_function = get_gamma_rate,
    use_box_cox_before_eta = TRUE,
    use_box_cox_in_eta = TRUE,
    controls_params = list(param1 = 0.5, param2 = 0.5),
    cases_params = list(param1 = 0, param2 = 1)
  ),
  gamma_2_bc = list(
    file_name = "gamma_2_box_cox_inside_eta_box_cox_before_eta",
    case = "gamma",
    param_adjuster_function = get_gamma_rate,
    use_box_cox_before_eta = TRUE,
    use_box_cox_in_eta = TRUE,
    controls_params = list(param1 = 0.5, param2 = 0.5),
    cases_params = list(param1 = 0, param2 = 4)
  ),
  gamma_3_bc = list(
    file_name = "gamma_3_box_cox_inside_eta_box_cox_before_eta",
    case = "gamma",
    param_adjuster_function = get_gamma_rate,
    use_box_cox_before_eta = TRUE,
    use_box_cox_in_eta = TRUE,
    controls_params = list(param1 = 0.5, param2 = 0.5),
    cases_params = list(param1 = 0, param2 = 1/8)
  )
)


# Scenarios to simulate
all_configs = list(
  gaussian = gaussian_configs,
  lognormal = lognormal_configs,
  gamma = gamma_configs
)

run_parallel_simulations(all_configs, max_cores = NULL)  # NULL 


# 
# # Gaussian samples
# 
# simulate_bias_rmse("normal_1",
#                    case = "gaussian",
#                    param_adjuster_function = get_mux_bisection,
#                    controls_params = list(param1 = 0, param2 = 1),
#                    cases_params = list(param1 = 0, param2 = 1))
# 
# simulate_bias_rmse("normal_2",
#                    case = "gaussian",
#                    param_adjuster_function = get_mux_bisection,
#                    controls_params = list(param1 = 0, param2 = 1),
#                    cases_params = list(param1 = 0, param2 = 1.4))
# 
# 
# simulate_bias_rmse("normal_3",
#                    case = "gaussian",
#                    param_adjuster_function = get_mux_bisection,
#                    controls_params = list(param1 = 0, param2 = 1),
#                    cases_params = list(param1 = 0, param2 = 3))
# 
# # Lognormal samples
# 
# simulate_bias_rmse("lognormal_1_box_cox_inside_eta",
#                    case = "lognormal",
#                    param_adjuster_function = get_mux_bisection,
#                    use_box_cox_in_eta = TRUE,
#                    controls_params = list(param1 = 0, param2 = 1),
#                    cases_params = list(param1 = 0, param2 = 0.5))
# 
# simulate_bias_rmse("lognormal_2_box_cox_inside_eta",
#                    case = "lognormal",
#                    param_adjuster_function = get_mux_bisection,
#                    use_box_cox_in_eta = TRUE,
#                    controls_params = list(param1 = 0, param2 = 1),
#                    cases_params = list(param1 = 0, param2 = 3/2))
# 
# simulate_bias_rmse("lognormal_3_box_cox_inside_eta",
#                    case = "lognormal",
#                    param_adjuster_function = get_mux_bisection,
#                    use_box_cox_in_eta = TRUE,
#                    controls_params = list(param1 = 0, param2 = 1),
#                    cases_params = list(param1 = 0, param2 = 0.2))
# 
# simulate_bias_rmse("lognormal_4_box_cox_inside_eta",
#                    case = "lognormal",
#                    param_adjuster_function = get_mux_bisection,
#                    use_box_cox_in_eta = TRUE,
#                    controls_params = list(param1 = 0, param2 = 1),
#                    cases_params = list(param1 = 0, param2 = 2))
# 
# 
# 
# simulate_bias_rmse("lognormal_1_box_cox_inside_eta_box_cox_before_eta",
#                    case = "lognormal",
#                    param_adjuster_function = get_mux_bisection,
#                    use_box_cox_before_eta = TRUE,
#                    use_box_cox_in_eta = TRUE,
#                    controls_params = list(param1 = 0, param2 = 1),
#                    cases_params = list(param1 = 0, param2 = 0.5))
# 
# simulate_bias_rmse("lognormal_2_box_cox_inside_eta_box_cox_before_eta",
#                    case = "lognormal",
#                    param_adjuster_function = get_mux_bisection,
#                    use_box_cox_before_eta = TRUE,
#                    use_box_cox_in_eta = TRUE,
#                    controls_params = list(param1 = 0, param2 = 1),
#                    cases_params = list(param1 = 0, param2 = 3/2))
# 
# simulate_bias_rmse("lognormal_3_box_cox_inside_eta_box_cox_before_eta",
#                    case = "lognormal",
#                    param_adjuster_function = get_mux_bisection,
#                    use_box_cox_before_eta = TRUE,
#                    use_box_cox_in_eta = TRUE,
#                    controls_params = list(param1 = 0, param2 = 1),
#                    cases_params = list(param1 = 0, param2 = 0.2))
# 
# simulate_bias_rmse("lognormal_4_box_cox_inside_eta_box_cox_before_eta",
#                    case = "lognormal",
#                    param_adjuster_function = get_mux_bisection,
#                    use_box_cox_before_eta = TRUE,
#                    use_box_cox_in_eta = TRUE,
#                    controls_params = list(param1 = 0, param2 = 1),
#                    cases_params = list(param1 = 0, param2 = 2))
# 
# # Gamma samples
# 
# simulate_bias_rmse("gamma_1_box_cox_inside_eta",
#                    case = "gamma",
#                    param_adjuster_function = get_gamma_rate,
#                    use_box_cox_in_eta = TRUE,
#                    controls_params = list(param1 = 0.5, param2 = 0.5),
#                    cases_params = list(param1 = 0, param2 = 1))
# 
# simulate_bias_rmse("gamma_2_box_cox_inside_eta",
#                    case = "gamma",
#                    param_adjuster_function = get_gamma_rate,
#                    use_box_cox_in_eta = TRUE,
#                    controls_params = list(param1 = 0.5, param2 = 0.5),
#                    cases_params = list(param1 = 0, param2 = 4))
# 
# simulate_bias_rmse("gamma_3_box_cox_inside_eta",
#                    case = "gamma",
#                    param_adjuster_function = get_gamma_rate,
#                    use_box_cox_in_eta = TRUE,
#                    controls_params = list(param1 = 0.5, param2 = 0.5),
#                    cases_params = list(param1 = 0, param2 = 1/8))
# 
# simulate_bias_rmse("gamma_1_box_cox_inside_eta_box_cox_before_eta",
#                    case = "gamma",
#                    param_adjuster_function = get_gamma_rate,
#                    use_box_cox_before_eta = TRUE,
#                    use_box_cox_in_eta = TRUE,
#                    controls_params = list(param1 = 0.5, param2 = 0.5),
#                    cases_params = list(param1 = 0, param2 = 1))
# 
# simulate_bias_rmse("gamma_2_box_cox_inside_eta_box_cox_before_eta",
#                    case = "gamma",
#                    param_adjuster_function = get_gamma_rate,
#                    use_box_cox_before_eta = TRUE,
#                    use_box_cox_in_eta = TRUE,
#                    controls_params = list(param1 = 0.5, param2 = 0.5),
#                    cases_params = list(param1 = 0, param2 = 4))
# 
# simulate_bias_rmse("gamma_3_box_cox_inside_eta_box_cox_before_eta",
#                    case = "gamma",
#                    param_adjuster_function = get_gamma_rate,
#                    use_box_cox_before_eta = TRUE,
#                    use_box_cox_in_eta = TRUE,
#                    controls_params = list(param1 = 0.5, param2 = 0.5),
#                    cases_params = list(param1 = 0, param2 = 1/8))
# 
# 
# 
# 
# 
# 
# 


