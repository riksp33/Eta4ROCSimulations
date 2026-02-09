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
get_p_values = function(controls, cases, n_perm = 500) {

    set.seed(1)

    n_controls = length(controls)
    n_cases = length(cases)
    combined = c(controls, cases)

    transformed_original = apply_box_cox(controls, cases)
    transformed_controls = transformed_original$transformed_x
    transformed_cases = transformed_original$transformed_y
    observed = list(
    auc_par_no_bc = max(calculate_auc_normal(cases, controls), calculate_auc_normal(controls, cases)),
    auc_par_yes_bc = max(calculate_auc_normal(transformed_cases, transformed_controls), calculate_auc_normal(transformed_controls, transformed_cases)),
    auc_hscv_no_bc = max(calculate_auc_kernel(cases, controls, "hscv", box_cox = FALSE), calculate_auc_kernel(controls, cases, "hscv", box_cox = FALSE)),
    auc_hscv_yes_bc = max(calculate_auc_kernel(cases, controls, "hscv", box_cox = TRUE), calculate_auc_kernel(controls, cases, "hscv", box_cox = TRUE)),
    auc_opt_no_bc = max(calculate_auc_kernel(cases, controls, "optimal", box_cox = FALSE), calculate_auc_kernel(controls, cases, "optimal", box_cox = FALSE)),
    auc_opt_yes_bc = max(calculate_auc_kernel(cases, controls, "optimal", box_cox = TRUE), calculate_auc_kernel(controls, cases, "optimal", box_cox = TRUE)),
    auc_iqr_no_bc = max(calculate_auc_kernel(cases, controls, "iqr", box_cox = FALSE), calculate_auc_kernel(controls, cases, "iqr", box_cox = FALSE)),
    auc_iqr_yes_bc = max(calculate_auc_kernel(cases, controls, "iqr", box_cox = TRUE), calculate_auc_kernel(controls, cases, "iqr", box_cox = TRUE)),
    youden_par_no_bc = max(calculate_youden_normal(cases, controls), calculate_youden_normal(controls, cases)),
    youden_par_yes_bc = max(calculate_youden_normal(transformed_cases, transformed_controls), calculate_youden_normal(transformed_controls, transformed_cases)),
    youden_hscv_no_bc = max(calculate_youden_kernel(cases, controls, "hscv", box_cox = FALSE), calculate_youden_kernel(controls, cases, "hscv", box_cox = FALSE)),
    youden_hscv_yes_bc = max(calculate_youden_kernel(cases, controls, "hscv", box_cox = TRUE), calculate_youden_kernel(controls, cases, "hscv", box_cox = TRUE)),
    youden_opt_no_bc = max(calculate_youden_kernel(cases, controls, "optimal", box_cox = FALSE), calculate_youden_kernel(controls, cases, "optimal", box_cox = FALSE)),
    youden_opt_yes_bc = max(calculate_youden_kernel(cases, controls, "optimal", box_cox = TRUE), calculate_youden_kernel(controls, cases, "optimal", box_cox = TRUE)),
    youden_iqr_no_bc = max(calculate_youden_kernel(cases, controls, "iqr", box_cox = FALSE), calculate_youden_kernel(controls, cases, "iqr", box_cox = FALSE)),
    youden_iqr_yes_bc = max(calculate_youden_kernel(cases, controls, "iqr", box_cox = TRUE), calculate_youden_kernel(controls, cases, "iqr", box_cox = TRUE)),
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

      perm_sample = sample(combined, replace = TRUE)
      controls_p = perm_sample[1:n_controls]
      cases_p = perm_sample[(n_controls + 1):(n_controls + n_cases)]

      transformed = apply_box_cox(controls_p, cases_p)
      transformed_controls_p = transformed$transformed_x
      transformed_cases_p = transformed$transformed_y



      perm$auc_par_yes_bc[b] = max(calculate_auc_normal(transformed_cases_p, transformed_controls_p), calculate_auc_normal(transformed_controls_p, transformed_cases_p))
      perm$auc_par_no_bc[b] = max(calculate_auc_normal(controls_p, cases_p), calculate_auc_normal(controls_p, cases_p))
      perm$auc_hscv_no_bc[b] = max(calculate_auc_kernel(cases_p, controls_p, "hscv", box_cox = FALSE), calculate_auc_kernel(controls_p, cases_p, "hscv", box_cox = FALSE))
      perm$auc_hscv_yes_bc[b] = max(calculate_auc_kernel(cases_p, controls_p, "hscv", box_cox = TRUE), calculate_auc_kernel(controls_p, cases_p, "hscv", box_cox = TRUE))
      perm$auc_opt_no_bc[b] = max(calculate_auc_kernel(cases_p, controls_p, "optimal", box_cox = FALSE), calculate_auc_kernel(controls_p, cases_p, "optimal", box_cox = FALSE))
      perm$auc_opt_yes_bc[b] = max(calculate_auc_kernel(cases_p, controls_p, "optimal", box_cox = TRUE), calculate_auc_kernel(controls_p, cases_p, "optimal", box_cox = TRUE))
      perm$auc_iqr_no_bc[b] = max(calculate_auc_kernel(cases_p, controls_p, "iqr", box_cox = FALSE), calculate_auc_kernel(controls_p, cases_p, "iqr", box_cox = FALSE))
      perm$auc_iqr_yes_bc[b] = max(calculate_auc_kernel(cases_p, controls_p, "iqr", box_cox = TRUE), calculate_auc_kernel(controls_p, cases_p, "iqr", box_cox = TRUE))
      perm$youden_par_yes_bc[b] = max(calculate_youden_normal(transformed_cases_p, transformed_controls_p), calculate_youden_normal(transformed_controls_p, transformed_cases_p))
      perm$youden_par_no_bc[b] = max(calculate_youden_normal(controls_p, cases_p), calculate_youden_normal(controls_p, cases_p))
      perm$youden_hscv_no_bc[b] = max(calculate_youden_kernel(cases_p, controls_p, "hscv", box_cox = FALSE), calculate_youden_kernel(controls_p, cases_p, "hscv", box_cox = FALSE))
      perm$youden_hscv_yes_bc[b] = max(calculate_youden_kernel(cases_p, controls_p, "hscv", box_cox = TRUE), calculate_youden_kernel(controls_p, cases_p, "hscv", box_cox = TRUE))
      perm$youden_opt_no_bc[b] = max(calculate_youden_kernel(cases_p, controls_p, "optimal", box_cox = FALSE), calculate_youden_kernel(controls_p, cases_p, "optimal", box_cox = FALSE))
      perm$youden_opt_yes_bc[b] = max(calculate_youden_kernel(cases_p, controls_p, "optimal", box_cox = TRUE), calculate_youden_kernel(controls_p, cases_p, "optimal", box_cox = TRUE))
      perm$youden_iqr_no_bc[b] = max(calculate_youden_kernel(cases_p, controls_p, "iqr", box_cox = FALSE), calculate_youden_kernel(controls_p, cases_p, "iqr", box_cox = FALSE))
      perm$youden_iqr_yes_bc[b] = max(calculate_youden_kernel(cases_p, controls_p, "iqr", box_cox = TRUE), calculate_youden_kernel(controls_p, cases_p, "iqr", box_cox = TRUE))
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

    print_test("AUC", observed$auc_gaussian, pvalues$auc_gaussian, "AUC = 0.5")
    print_test("AUC (hscv, no BC)", observed$auc_hscv_no_bc, pvalues$auc_hscv_no_bc, "AUC = 0.5")
    print_test("AUC (hscv, yes BC)", observed$auc_hscv_yes_bc, pvalues$auc_hscv_yes_bc, "AUC = 0.5")
    print_test("AUC (optimal, no BC)", observed$auc_opt_no_bc, pvalues$auc_opt_no_bc, "AUC = 0.5")
    print_test("AUC (optimal, yes BC)", observed$auc_opt_yes_bc, pvalues$auc_opt_yes_bc, "AUC = 0.5")
    print_test("AUC (iqr, no BC)", observed$auc_iqr_no_bc, pvalues$auc_iqr_no_bc, "AUC = 0.5")
    print_test("AUC (iqr, yes BC)", observed$auc_iqr_yes_bc, pvalues$auc_iqr_yes_bc, "AUC = 0.5")
    print_test("Youden", observed$youden_gaussian, pvalues$youden_gaussian, "Youden = 0")
    print_test("Youden (hscv, no BC)", observed$youden_hscv_no_bc, pvalues$youden_hscv_no_bc, "Youden = 0")
    print_test("Youden (hscv, yes BC)", observed$youden_hscv_yes_bc, pvalues$youden_hscv_yes_bc, "Youden = 0")
    print_test("Youden (optimal, no BC)", observed$youden_opt_no_bc, pvalues$youden_opt_no_bc, "Youden = 0")
    print_test("Youden (optimal, yes BC)", observed$youden_opt_yes_bc, pvalues$youden_opt_yes_bc, "Youden = 0")
    print_test("Youden (iqr, no BC)", observed$youden_iqr_no_bc, pvalues$youden_iqr_no_bc, "Youden = 0")
    print_test("Youden (iqr, yes BC)", observed$youden_iqr_yes_bc, pvalues$youden_iqr_yes_bc, "Youden = 0")
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
get_confidence_interval = function(controls, cases, conf_level = 0.95, n_boot = 500) {

  set.seed(1)

  alpha = 1 - conf_level

  transformed_original = apply_box_cox(controls, cases)
  transformed_controls = transformed_original$transformed_x
  transformed_cases = transformed_original$transformed_y

  observed = list(
    auc_par_no_bc = max(calculate_auc_normal(cases, controls), calculate_auc_normal(controls, cases)),
    auc_par_yes_bc = max(calculate_auc_normal(transformed_cases, transformed_controls), calculate_auc_normal(transformed_controls, transformed_cases)),
    auc_hscv_no_bc = max(calculate_auc_kernel(cases, controls, "hscv", box_cox = FALSE), calculate_auc_kernel(controls, cases, "hscv", box_cox = FALSE)),
    auc_hscv_yes_bc = max(calculate_auc_kernel(cases, controls, "hscv", box_cox = TRUE), calculate_auc_kernel(controls, cases, "hscv", box_cox = TRUE)),
    auc_opt_no_bc = max(calculate_auc_kernel(cases, controls, "optimal", box_cox = FALSE), calculate_auc_kernel(controls, cases, "optimal", box_cox = FALSE)),
    auc_opt_yes_bc = max(calculate_auc_kernel(cases, controls, "optimal", box_cox = TRUE), calculate_auc_kernel(controls, cases, "optimal", box_cox = TRUE)),
    auc_iqr_no_bc = max(calculate_auc_kernel(cases, controls, "iqr", box_cox = FALSE), calculate_auc_kernel(controls, cases, "iqr", box_cox = FALSE)),
    auc_iqr_yes_bc = max(calculate_auc_kernel(cases, controls, "iqr", box_cox = TRUE), calculate_auc_kernel(controls, cases, "iqr", box_cox = TRUE)),
    youden_par_no_bc = max(calculate_youden_normal(cases, controls), calculate_youden_normal(controls, cases)),
    youden_par_yes_bc = max(calculate_youden_normal(transformed_cases, transformed_controls), calculate_youden_normal(transformed_controls, transformed_cases)),
    youden_hscv_no_bc = max(calculate_youden_kernel(cases, controls, "hscv", box_cox = FALSE), calculate_youden_kernel(controls, cases, "hscv", box_cox = FALSE)),
    youden_hscv_yes_bc = max(calculate_youden_kernel(cases, controls, "hscv", box_cox = TRUE), calculate_youden_kernel(controls, cases, "hscv", box_cox = TRUE)),
    youden_opt_no_bc = max(calculate_youden_kernel(cases, controls, "optimal", box_cox = FALSE), calculate_youden_kernel(controls, cases, "optimal", box_cox = FALSE)),
    youden_opt_yes_bc = max(calculate_youden_kernel(cases, controls, "optimal", box_cox = TRUE), calculate_youden_kernel(controls, cases, "optimal", box_cox = TRUE)),
    youden_iqr_no_bc = max(calculate_youden_kernel(cases, controls, "iqr", box_cox = FALSE), calculate_youden_kernel(controls, cases, "iqr", box_cox = FALSE)),
    youden_iqr_yes_bc = max(calculate_youden_kernel(cases, controls, "iqr", box_cox = TRUE), calculate_youden_kernel(controls, cases, "iqr", box_cox = TRUE)),
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

    transformed_cases_b = apply_box_cox(cases_b, controls_b)$transformed_y
    transformed_controls_b = apply_box_cox(cases_b, controls_b)$transformed_x


    boot$auc_par_no_bc[b] = max(calculate_auc_normal(cases_b, controls_b), calculate_auc_normal(controls_b, cases_b))
    boot$auc_par_yes_bc[b] = max(calculate_auc_normal(transformed_cases_b, transformed_controls_b), calculate_auc_normal(transformed_controls_b, transformed_cases_b))
    boot$auc_hscv_no_bc[b] = max(calculate_auc_kernel(cases_b, controls_b, "hscv", box_cox = FALSE), calculate_auc_kernel(controls_b, cases_b, "hscv", box_cox = FALSE))
    boot$auc_hscv_yes_bc[b] = max(calculate_auc_kernel(cases_b, controls_b, "hscv", box_cox = TRUE), calculate_auc_kernel(controls_b, cases_b, "hscv", box_cox = TRUE))
    boot$auc_opt_no_bc[b] = max(calculate_auc_kernel(cases_b, controls_b, "optimal", box_cox = FALSE), calculate_auc_kernel(controls_b, cases_b, "optimal", box_cox = FALSE))
    boot$auc_opt_yes_bc[b] = max(calculate_auc_kernel(cases_b, controls_b, "optimal", box_cox = TRUE), calculate_auc_kernel(controls_b, cases_b, "optimal", box_cox = TRUE))
    boot$auc_iqr_no_bc[b] = max(calculate_auc_kernel(cases_b, controls_b, "iqr", box_cox = FALSE), calculate_auc_kernel(controls_b, cases_b, "iqr", box_cox = FALSE))
    boot$auc_iqr_yes_bc[b] = max(calculate_auc_kernel(cases_b, controls_b, "iqr", box_cox = TRUE), calculate_auc_kernel(controls_b, cases_b, "iqr", box_cox = TRUE))
    boot$youden_par_no_bc[b] = max(calculate_youden_normal(cases_b, controls_b), calculate_youden_normal(controls_b, cases_b))
    boot$youden_par_yes_bc[b] = max(calculate_youden_normal(transformed_cases_b, transformed_controls_b), calculate_youden_normal(transformed_controls_b, transformed_cases_b))
    boot$youden_hscv_no_bc[b] = max(calculate_youden_kernel(cases_b, controls_b, "hscv", box_cox = FALSE), calculate_youden_kernel(controls_b, cases_b, "hscv", box_cox = FALSE))
    boot$youden_hscv_yes_bc[b] = max(calculate_youden_kernel(cases_b, controls_b, "hscv", box_cox = TRUE), calculate_youden_kernel(controls_b, cases_b, "hscv", box_cox = TRUE))
    boot$youden_opt_no_bc[b] = max(calculate_youden_kernel(cases_b, controls_b, "optimal", box_cox = FALSE), calculate_youden_kernel(controls_b, cases_b, "optimal", box_cox = FALSE))
    boot$youden_opt_yes_bc[b] = max(calculate_youden_kernel(cases_b, controls_b, "optimal", box_cox = TRUE), calculate_youden_kernel(controls_b, cases_b, "optimal", box_cox = TRUE))
    boot$youden_iqr_no_bc[b] = max(calculate_youden_kernel(cases_b, controls_b, "iqr", box_cox = FALSE), calculate_youden_kernel(controls_b, cases_b, "iqr", box_cox = FALSE))
    boot$youden_iqr_yes_bc[b] = max(calculate_youden_kernel(cases_b, controls_b, "iqr", box_cox = TRUE), calculate_youden_kernel(controls_b, cases_b, "iqr", box_cox = TRUE))
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
write_json(results_all, "results/simulations_all_kernel_for_all_v2.json", pretty = TRUE, auto_unbox = TRUE)
