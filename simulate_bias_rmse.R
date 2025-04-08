# Setup needed to install the Eta4ROC package from github
library(remotes)
remotes::install_github("riksp33/Eta4ROC")
library(Eta4ROC)

# Setup needed to generate the files
library(jsonlite)
library(here)


simulate_bias_rmse = function(file_name,
                              use_box_cox_in_eta = FALSE,
                              use_box_cox_before_eta = FALSE,
                              param_adjuster_function = function(mu){mu},
                              controls_params = list(param1 = 1, param2 = 1), #Params for cases (mean, std) for gaussian & lognormal, (shape, rate) for gamma 
                              cases_params = list(param1 = 1.1, param2 = 1,1) #Params for controls (mean, std) for gaussian & lognormal, (shape, rate) for gamma 
                              ){

  # Set the seed
  set.seed(1)
  
  # Output JSON
  ret_json = list()
  ret_json[["header"]] = paste("Simulation results for ", file_name)

  # Simulation params
  MC = 1000
  AUCs = c(0.6, 0.75, 0.9)
  ns = c(20, 50, 100)
  t0s = c(0.2, 0.4, 0.8, 1)
  
  # Sample generator function
  switch(case,
          normal = {
            sample_distribution = function(n, mu , sigma) {rnorm(n, mu, sigma)}
          },
          lognormal = {
            sample_distribution = function(n, mu, sigma) {rlnorm(n, mu, sigma)}
          },
          gamma ={
            sample_distribution = function(n, shape, rate) {rgamma(n, shape, rate = rate)}
          })

  for(auc in AUCs){
    ret_aucs = list()    

    missing_param = param_adjuster_function(sample_distribution(1e5, controls_params$param1, controls_params$param2),
                                            cases_params$param2,
                                            auc)
    
    for(t0 in t0s){
      ret_t0 = list()

      # True eta
      true_eta = analytical_eta(controls_params&param1,
                                controls_params$param2,
                                cases_params$param1,
                                cases_params$param2,
                                case = case,
                                t0)
      
      for(n in ns){
        
        parametric = numeric(MC)
        kernel_hscv = numeric(MC)
        kernel_opt = numeric(MC)
        kernel_iqr = numeric(MC)
        
        for(i in 1:MC){
          
          # Generate samples
          controls = sample_distribution(n, controls_params$param1, controls_params$param2)
          cases = sample_distribution(n, cases_params$param1, cases_params$param2)
          
          if(use_box_cox_before_eta){
            transformed = apply_box_cox(controls, cases)
            controls = transformed$transformed_x
            cases = transfromed$transformed_y
          }
          
          # NOTA: Cambiar documentacion para hablar sobre la opcion iqr
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
      ret_aucs[[paste("t0:", t0)]]
    }
    ret_json[[paste("AUC:", auc)]]
  }
  json = toJSON(ret_json, pretty = TRUE, digits = NA)
  dir_path = here("results")
  full_path = file.path(dir_path, file_name)
  write(json, file = full_path)
}


simulate_bias_rmse("pruebas",
                  param_adjuster_function = get_mux_bisection)