library(remotes)
remotes::install_github("riksp33/Eta4ROC")
library(Eta4ROC)
library(jsonlite)

print_test = function(name, observed, pvalue, null_desc) {
  cat(sprintf(
    "%s\n  Observed statistic: %.6f\n  Null hypothesis: %s\n  Permutation p-value: %.6f\n\n",
    name, observed, null_desc, pvalue
  ))
}

print_ci = function(name, observed, ci, conf_level) {
  cat(sprintf(
    "%s\n  Observed statistic: %.6f\n  %.1f%% percentile CI: [%.6f, %.6f]\n\n",
    name, observed, 100 * conf_level, ci[1], ci[2]
  ))
}

# Accumulator for all simulation results
results_all <- list()

# Helper to add results for a file path and file name
add_results <- function(file_path, p_res, ci_res) {
  fname <- basename(file_path)
  p_list <- list()
  ci_list <- list()

  for (nm in names(p_res$observed)) {
    p_list[[nm]] <- p_res$pvalues[[nm]]
    ci_vals <- ci_res$ci[[nm]]
    ci_list[[nm]] <- list(
      value = ci_res$observed[[nm]],
      lower = ci_vals[1],
      upper = ci_vals[2]
    )
  }

  if (is.null(results_all[[file_path]])) {
    results_all[[file_path]] <<- list()
  }

  results_all[[file_path]][[fname]] <<- list(
    p_values = p_list,
    confidence_intervals = ci_list
  )
}


############################################
## Permutation tests
############################################
get_p_values = function(controls, cases, n_perm = 2) {

    set.seed(1)

    n_controls = length(controls)
    n_cases = length(cases)
    combined = c(controls, cases)

    observed = list(
    auc = max(calculate_auc_normal(cases, controls), calculate_auc_normal(controls, cases)),
    youden = max(calculate_youden_normal(cases, controls), calculate_youden_normal(controls, cases)),
    eta_par_no_bc = parametric_eta(controls, cases, 1, box_cox = FALSE),
    eta_par_yes_bc = parametric_eta(controls, cases, 1, box_cox = TRUE),
    eta_hscv_no_bc = kernel_eta(controls, cases, "hscv", 1, box_cox = FALSE),
    eta_hscv_yes_bc = kernel_eta(controls, cases, "hscv", 1, box_cox = TRUE),
    eta_opt_no_bc = kernel_eta(controls, cases, "optimal", 1, box_cox = FALSE),
    eta_opt_yes_bc = kernel_eta(controls, cases, "optimal", 1, box_cox = TRUE),
    eta_iqr_no_bc = kernel_eta(controls, cases, "iqr", 1, box_cox = FALSE),
    eta_iqr_yes_bc = kernel_eta(controls, cases, "iqr", 1, box_cox = TRUE)
    )

    # Generates the empty dataframe of numerics(n_perm)
    perm = lapply(observed, function(x) numeric(n_perm))

    for (b in seq_len(n_perm)) {

    perm_sample = sample(combined, replace = FALSE)
    controls_p = perm_sample[1:n_controls]
    cases_p = perm_sample[(n_controls + 1):(n_controls + n_cases)]

    perm$auc[b] = max(calculate_auc_normal(cases_p, controls_p), calculate_auc_normal(controls_p, cases_p))
    perm$youden[b] = max(calculate_youden_normal(cases_p, controls_p), calculate_youden_normal(controls_p, cases_p))
    perm$eta_par_no_bc[b] = parametric_eta(controls_p, cases_p, 1, box_cox = FALSE)
    perm$eta_par_yes_bc[b] = parametric_eta(controls_p, cases_p, 1, box_cox = TRUE)
    perm$eta_hscv_no_bc[b] = kernel_eta(controls_p, cases_p, "hscv", 1, box_cox = FALSE)
    perm$eta_hscv_yes_bc[b] = kernel_eta(controls_p, cases_p, "hscv", 1, box_cox = TRUE)
    perm$eta_opt_no_bc[b] = kernel_eta(controls_p, cases_p, "optimal", 1, box_cox = FALSE)
    perm$eta_opt_yes_bc[b] = kernel_eta(controls_p, cases_p, "optimal", 1, box_cox = TRUE)
    perm$eta_iqr_no_bc[b] = kernel_eta(controls_p, cases_p, "iqr", 1, box_cox = FALSE)
    perm$eta_iqr_yes_bc[b] = kernel_eta(controls_p, cases_p, "iqr", 1, box_cox = TRUE)
    }

    pvalues = mapply(
    function(obs, perm_dist) {
        sum(perm_dist >= obs) / n_perm
    },
    observed,
    perm,
    SIMPLIFY = FALSE
    )

    cat("============================================\n")
    cat(" PERMUTATION TESTS (NO DISCRIMINATION)\n")
    cat("============================================\n\n")

    print_test("AUC", observed$auc, pvalues$auc, "AUC = 0.5")
    print_test("Youden", observed$youden, pvalues$youden, "Youden = 0")
    print_test("Parametric ETA (no BC)", observed$eta_par_no_bc, pvalues$eta_par_no_bc, "ETA = 0")
    print_test("Parametric ETA (yes BC)", observed$eta_par_yes_bc, pvalues$eta_par_yes_bc, "ETA = 0")
    print_test("Kernel ETA hscv (no BC)", observed$eta_hscv_no_bc, pvalues$eta_hscv_no_bc, "ETA = 0")
    print_test("Kernel ETA hscv (yes BC)", observed$eta_hscv_yes_bc, pvalues$eta_hscv_yes_bc, "ETA = 0")
    print_test("Kernel ETA optimal (no BC)", observed$eta_opt_no_bc, pvalues$eta_opt_no_bc, "ETA = 0")
    print_test("Kernel ETA optimal (yes BC)", observed$eta_opt_yes_bc, pvalues$eta_opt_yes_bc, "ETA = 0")
    print_test("Kernel ETA iqr (no BC)", observed$eta_iqr_no_bc, pvalues$eta_iqr_no_bc, "ETA = 0")
    print_test("Kernel ETA iqr (yes BC)", observed$eta_iqr_yes_bc, pvalues$eta_iqr_yes_bc, "ETA = 0")

    invisible(list(observed = observed, pvalues = pvalues))
}


############################################
## Bootstrap confidence intervals
############################################
get_confidence_interval = function(controls, cases, conf_level = 0.95, n_boot = 2) {

  set.seed(1)

  alpha = 1 - conf_level

  observed = list(
    auc = max(calculate_auc_normal(cases, controls), calculate_auc_normal(controls, cases)),
    youden = max(calculate_youden_normal(cases, controls), calculate_youden_normal(controls, cases)),
    eta_par_no_bc = parametric_eta(controls, cases, 1, box_cox = FALSE),
    eta_par_yes_bc = parametric_eta(controls, cases, 1, box_cox = TRUE),
    eta_hscv_no_bc = kernel_eta(controls, cases, "hscv", 1, box_cox = FALSE),
    eta_hscv_yes_bc = kernel_eta(controls, cases, "hscv", 1, box_cox = TRUE),
    eta_opt_no_bc = kernel_eta(controls, cases, "optimal", 1, box_cox = FALSE),
    eta_opt_yes_bc = kernel_eta(controls, cases, "optimal", 1, box_cox = TRUE),
    eta_iqr_no_bc = kernel_eta(controls, cases, "iqr", 1, box_cox = FALSE),
    eta_iqr_yes_bc = kernel_eta(controls, cases, "iqr", 1, box_cox = TRUE)
  )

  boot = lapply(observed, function(x) numeric(n_boot))

  for (b in seq_len(n_boot)) {

    controls_b = sample(controls, replace = TRUE)
    cases_b = sample(cases, replace = TRUE)

    boot$auc[b] = max(calculate_auc_normal(cases_b, controls_b), calculate_auc_normal(controls_b, cases_b))
    boot$youden[b] = max(calculate_youden_normal(cases_b, controls_b), calculate_youden_normal(controls_b, cases_b))
    boot$eta_par_no_bc[b] = parametric_eta(controls_b, cases_b, 1, box_cox = FALSE)
    boot$eta_par_yes_bc[b] = parametric_eta(controls_b, cases_b, 1, box_cox = TRUE)
    boot$eta_hscv_no_bc[b] = kernel_eta(controls_b, cases_b, "hscv", 1, box_cox = FALSE)
    boot$eta_hscv_yes_bc[b] = kernel_eta(controls_b, cases_b, "hscv", 1, box_cox = TRUE)
    boot$eta_opt_no_bc[b] = kernel_eta(controls_b, cases_b, "optimal", 1, box_cox = FALSE)
    boot$eta_opt_yes_bc[b] = kernel_eta(controls_b, cases_b, "optimal", 1, box_cox = TRUE)
    boot$eta_iqr_no_bc[b] = kernel_eta(controls_b, cases_b, "iqr", 1, box_cox = FALSE)
    boot$eta_iqr_yes_bc[b] = kernel_eta(controls_b, cases_b, "iqr", 1, box_cox = TRUE)
  }

  ci = lapply(
    boot,
    function(x) quantile(x, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  )

  cat("============================================\n")
  cat(" BOOTSTRAP CONFIDENCE INTERVALS\n")
  cat("============================================\n\n")

  for (name in names(observed)) {
    print_ci(name, observed[[name]], ci[[name]], conf_level)
  }

  invisible(list(observed = observed, ci = ci))
}


############################################
## DATASET 207039N
############################################
data_207039N = read.table(
  "/Users/Riki/Desktop/ucm/articulo/Simulciones/Ejemplo Bantis0 compartido/53,53/data207039.txt",
  header = FALSE
)

controls_207039N = data_207039N$V1[data_207039N$V2 == 0]
cases_207039N = data_207039N$V1[data_207039N$V2 == 1]

cat("\n========== DATASET 207039N ==========\n\n")
cat("/Users/Riki/Desktop/ucm/articulo/Simulciones/Ejemplo Bantis0 compartido/53,53/data207039.txt \n \n")
file_path <- "/Users/Riki/Desktop/ucm/articulo/Simulciones/Ejemplo Bantis0 compartido/53,53/data207039.txt"
p_res <- get_p_values(controls_207039N, cases_207039N)
ci_res <- get_confidence_interval(controls_207039N, cases_207039N)
add_results(file_path, p_res, ci_res)


############################################
## DATASET 207039N
############################################
data_207039N = read.table(
  "/Users/Riki/Desktop/ucm/articulo/Simulciones/Ejemplo Bantis0 compartido/data207039.txt",
  header = FALSE
)

controls_207039N = data_207039N$V1[data_207039N$V2 == 0]
cases_207039N = data_207039N$V1[data_207039N$V2 == 1]

cat("\n========== DATASET 207039N ==========\n\n")
cat("/Users/Riki/Desktop/ucm/articulo/Simulciones/Ejemplo Bantis0 compartido/data207039.txt \n \n")
file_path <- "/Users/Riki/Desktop/ucm/articulo/Simulciones/Ejemplo Bantis0 compartido/data207039.txt"
p_res <- get_p_values(controls_207039N, cases_207039N)
ci_res <- get_confidence_interval(controls_207039N, cases_207039N)
add_results(file_path, p_res, ci_res)


############################################
## DATASET 209644N
############################################
data_209644N = read.table(
  "/Users/Riki/Desktop/ucm/articulo/Simulciones/Ejemplo Bantis0 compartido/53,53/data209644.txt",
  header = FALSE
)

controls_209644N = data_209644N$V1[data_209644N$V2 == 0]
cases_209644N = data_209644N$V1[data_209644N$V2 == 1]

cat("\n========== DATASET 209644N ==========\n\n")
file_path <- "/Users/Riki/Desktop/ucm/articulo/Simulciones/Ejemplo Bantis0 compartido/53,53/data209644.txt"
p_res <- get_p_values(controls_209644N, cases_209644N)
ci_res <- get_confidence_interval(controls_209644N, cases_209644N)
add_results(file_path, p_res, ci_res)

############################################
## DATASET 209644N
############################################
data_209644N = read.table(
  "/Users/Riki/Desktop/ucm/articulo/Simulciones/Ejemplo Bantis0 compartido/data209644.txt",
  header = FALSE
)

controls_209644N = data_209644N$V1[data_209644N$V2 == 0]
cases_209644N = data_209644N$V1[data_209644N$V2 == 1]

cat("\n========== DATASET 209644N ==========\n\n")
cat("/Users/Riki/Desktop/ucm/articulo/Simulciones/Ejemplo Bantis0 compartido/data209644.txt \n \n")
file_path <- "/Users/Riki/Desktop/ucm/articulo/Simulciones/Ejemplo Bantis0 compartido/data209644.txt"
p_res <- get_p_values(controls_209644N, cases_209644N)
ci_res <- get_confidence_interval(controls_209644N, cases_209644N)
add_results(file_path, p_res, ci_res)

# Write accumulated results to JSON file
write_json(results_all, "results/simulations_all.json", pretty = TRUE, auto_unbox = TRUE)
