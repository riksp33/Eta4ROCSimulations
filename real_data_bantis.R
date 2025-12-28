# Setup needed to install the Eta4ROC package from github
library(remotes)
remotes::install_github("riksp33/Eta4ROC")
library(Eta4ROC)


get_confidence_interval = function(controls, cases, conf_level = 0.95, n_boot = 500) {

    combined_sample = c(controls, cases)

    auc = numeric(n_boot)
    youden = numeric(n_boot)
    parametric_eta_no_bc = numeric(n_boot)
    parametric_eta_yes_bc = numeric(n_boot)
    kernel_eta_hscv_no_bc = numeric(n_boot)
    kernel_eta_hscv_yes_bc = numeric(n_boot)
    kernel_eta_optim_no_bc = numeric(n_boot)
    kernel_eta_optim_yes_bc = numeric(n_boot)
    kernel_eta_iqr_no_bc = numeric(n_boot)
    kernel_eta_iqr_yes_bc = numeric(n_boot)

    for (b_it in 1:n_boot){
        controls_size = length(controls)
        cases_size = length(cases)
        bootstrap_sample = sample(combined_sample, replace = FALSE)
        bootstrap_controls = bootstrap_sample[1:controls_size]
        bootstrap_cases = bootstrap_sample[controls_size: length(combined_sample)]

        auc[b_it] = max(
            calculate_auc_normal(bootstrap_cases, bootstrap_controls),
            calculate_auc_normal(bootstrap_controls, bootstrap_cases)
        )
        youden[b_it] = max(
            calculate_youden_normal(bootstrap_cases, bootstrap_controls),
            calculate_youden_normal(bootstrap_controls, bootstrap_cases)
        )

        parametric_eta_no_bc[b_it] = parametric_eta(bootstrap_controls, bootstrap_cases, 1, box_cox = FALSE)
        parametric_eta_yes_bc[b_it] = parametric_eta(bootstrap_controls, bootstrap_cases, 1, box_cox = TRUE)
        kernel_eta_hscv_no_bc[b_it] = kernel_eta(bootstrap_controls, bootstrap_cases, "hscv", 1, box_cox = FALSE)
        kernel_eta_hscv_yes_bc[b_it] = kernel_eta(bootstrap_controls, bootstrap_cases, "hscv", 1, box_cox = TRUE)
        kernel_eta_optim_no_bc[b_it] = kernel_eta(bootstrap_controls, bootstrap_cases, "optimal", 1, box_cox = FALSE)
        kernel_eta_optim_yes_bc[b_it] = kernel_eta(bootstrap_controls, bootstrap_cases, "optimal", 1, box_cox = TRUE)
        kernel_eta_iqr_no_bc[b_it] = kernel_eta(bootstrap_controls, bootstrap_cases, "iqr", 1, box_cox = FALSE)
        kernel_eta_iqr_yes_bc[b_it] = kernel_eta(bootstrap_controls, bootstrap_cases, "iqr", 1, box_cox = TRUE)
    }

    original_auc = max(
        calculate_auc_normal(cases, controls),
        calculate_auc_normal(controls, cases)
    )
    original_youden = max(
        calculate_youden_normal(cases, controls),
        calculate_youden_normal(controls, cases)
    )
    original_parametric_eta_no_bc = parametric_eta(controls, cases, 1, box_cox = FALSE)
    original_parametric_eta_yes_bc = parametric_eta(controls, cases, 1, box_cox = TRUE)
    original_kernel_eta_hscv_no_bc = kernel_eta(controls, cases, "hscv", 1, box_cox = FALSE)
    original_kernel_eta_hscv_yes_bc = kernel_eta(controls, cases, "hscv", 1, box_cox = TRUE)
    original_kernel_eta_optim_no_bc = kernel_eta(controls, cases, "optimal", 1, box_cox = FALSE)
    original_kernel_eta_optim_yes_bc = kernel_eta(controls, cases, "optimal", 1, box_cox = TRUE)
    original_kernel_eta_iqr_no_bc = kernel_eta(controls, cases, "iqr", 1, box_cox = FALSE)
    original_kernel_eta_iqr_yes_bc = kernel_eta(controls, cases, "iqr", 1, box_cox = TRUE)


    pvalue_auc = sum(auc > original_auc) / n_boot
    pvalue_youden = sum(youden > original_youden) / n_boot
    pvalue_parametric_eta_no_bc = sum(parametric_eta_no_bc > original_parametric_eta_no_bc) / n_boot
    pvalue_parametric_eta_yes_bc = sum(parametric_eta_yes_bc > original_parametric_eta_yes_bc) / n_boot
    pvalue_kernel_eta_hscv_no_bc = sum(kernel_eta_hscv_no_bc > original_kernel_eta_hscv_no_bc) / n_boot
    pvalue_kernel_eta_hscv_yes_bc = sum(kernel_eta_hscv_yes_bc > original_kernel_eta_hscv_yes_bc) / n_boot
    pvalue_kernel_eta_optim_no_bc = sum(kernel_eta_optim_no_bc > original_kernel_eta_optim_no_bc) / n_boot
    pvalue_kernel_eta_optim_yes_bc = sum(kernel_eta_optim_yes_bc > original_kernel_eta_optim_yes_bc) / n_boot
    pvalue_kernel_eta_iqr_no_bc = sum(kernel_eta_iqr_no_bc > original_kernel_eta_iqr_no_bc) / n_boot
    pvalue_kernel_eta_iqr_yes_bc = sum(kernel_eta_iqr_yes_bc > original_kernel_eta_iqr_yes_bc) / n_boot

    sprintf("pvalue auc: %f", pvalue_auc)
    sprintf("pvalue youden: %f", pvalue_youden)
    sprintf("pvalue parametric eta no bc: %f", pvalue_parametric_eta_no_bc)
    sprintf("pvalue parametric eta yes bc: %f", pvalue_parametric_eta_yes_bc)
    sprintf("pvalue kernel eta hscv no bc: %f", pvalue_kernel_eta_hscv_no_bc)
    sprintf("pvalue kernel eta hscv yes bc: %f", pvalue_kernel_eta_hscv_yes_bc)
    sprintf("pvalue kernel eta optim no bc: %f", pvalue_kernel_eta_optim_no_bc)
    sprintf("pvalue kernel eta optim yes bc: %f", pvalue_kernel_eta_optim_yes_bc)
    sprintf("pvalue kernel eta iqr no bc: %f", pvalue_kernel_eta_iqr_no_bc)
    sprintf("pvalue kernel eta iqr yes bc: %f", pvalue_kernel_eta_iqr_yes_bc)

    lower_ci_auc = quantile(auc, probs = (1 - conf_level) / 2)
    upper_ci_auc = quantile(auc, probs = 1 - (1 - conf_level) / 2)
    sprintf("AUC %f CI: [%f, %f ]", conf_level, lower_ci_auc, upper_ci_auc) 
    lower_ci_youden = quantile(youden, probs = (1 - conf_level) / 2)
    upper_ci_youden = quantile(youden, probs = 1 - (1 - conf_level) / 2)
    sprintf("Youden %f CI: [%f, %f ]", conf_level, lower_ci_youden, upper_ci_youden)
    lower_ci_parametric_eta_no_bc = quantile(parametric_eta_no_bc, probs = (1 - conf_level) / 2)
    upper_ci_parametric_eta_no_bc = quantile(parametric_eta_no_bc, probs = 1 - (1 - conf_level) / 2)
    sprintf("Parametric ETA no BC %f CI: [%f, %f ]", conf_level, lower_ci_parametric_eta_no_bc, upper_ci_parametric_eta_no_bc)
    lower_ci_parametric_eta_yes_bc = quantile(parametric_eta_yes_bc, probs = (1 - conf_level) / 2)
    upper_ci_parametric_eta_yes_bc = quantile(parametric_eta_yes_bc, probs = 1 - (1 - conf_level) / 2)
    sprintf("Parametric ETA yes BC %f CI: [%f, %f ]", conf_level, lower_ci_parametric_eta_yes_bc, upper_ci_parametric_eta_yes_bc)
    lower_ci_kernel_eta_hscv_no_bc = quantile(kernel_eta_hscv_no_bc, probs = (1 - conf_level) / 2)
    upper_ci_kernel_eta_hscv_no_bc = quantile(kernel_eta_hscv_no_bc, probs = 1 - (1 - conf_level) / 2)
    sprintf("Kernel ETA hscv no BC %f CI: [%f, %f ]", conf_level, lower_ci_kernel_eta_hscv_no_bc, upper_ci_kernel_eta_hscv_no_bc)
    lower_ci_kernel_eta_hscv_yes_bc = quantile(kernel_eta_hscv_yes_bc, probs = (1 - conf_level) / 2)
    upper_ci_kernel_eta_hscv_yes_bc = quantile(kernel_eta_hscv_yes_bc, probs = 1 - (1 - conf_level) / 2)
    sprintf("Kernel ETA hscv yes BC %f CI: [%f, %f ]", conf_level, lower_ci_kernel_eta_hscv_yes_bc, upper_ci_kernel_eta_hscv_yes_bc)
    lower_ci_kernel_eta_optim_no_bc = quantile(kernel_eta_optim_no_bc, probs = (1 - conf_level) / 2)
    upper_ci_kernel_eta_optim_no_bc = quantile(kernel_eta_optim_no_bc, probs = 1 - (1 - conf_level) / 2)
    sprintf("Kernel ETA optim no BC %f CI: [%f, %f ]", conf_level, lower_ci_kernel_eta_optim_no_bc, upper_ci_kernel_eta_optim_no_bc)
    lower_ci_kernel_eta_optim_yes_bc = quantile(kernel_eta_optim_yes_bc, probs = (1 - conf_level) / 2)
    upper_ci_kernel_eta_optim_yes_bc = quantile(kernel_eta_optim_yes_bc, probs = 1 - (1 - conf_level) / 2)
    sprintf("Kernel ETA optim yes BC %f CI: [%f, %f ]", conf_level, lower_ci_kernel_eta_optim_yes_bc, upper_ci_kernel_eta_optim_yes_bc)
    lower_ci_kernel_eta_iqr_no_bc = quantile(kernel_eta_iqr_no_bc, probs = (1 - conf_level) / 2)
    upper_ci_kernel_eta_iqr_no_bc = quantile(kernel_eta_iqr_no_bc, probs = 1 - (1 - conf_level) / 2)
    sprintf("Kernel ETA iqr no BC %f CI: [%f, %f ]", conf_level, lower_ci_kernel_eta_iqr_no_bc, upper_ci_kernel_eta_iqr_no_bc)
    lower_ci_kernel_eta_iqr_yes_bc = quantile(kernel_eta_iqr_yes_bc, probs = (1 - conf_level) / 2)
    upper_ci_kernel_eta_iqr_yes_bc = quantile(kernel_eta_iqr_yes_bc, probs = 1 - (1 - conf_level) / 2)
    sprintf("Kernel ETA iqr yes BC %f CI: [%f, %f ]", conf_level, lower_ci_kernel_eta_iqr_yes_bc, upper_ci_kernel_eta_iqr_yes_bc)


    return(list(
        pvalue_auc = pvalue_auc,
        pvalue_youden = pvalue_youden,
        pvalue_parametric_eta_no_bc = pvalue_parametric_eta_no_bc,
        pvalue_parametric_eta_yes_bc = pvalue_parametric_eta_yes_bc,
        pvalue_kernel_eta_hscv_no_bc = pvalue_kernel_eta_hscv_no_bc,
        pvalue_kernel_eta_hscv_yes_bc = pvalue_kernel_eta_hscv_yes_bc,
        pvalue_kernel_eta_optim_no_bc = pvalue_kernel_eta_optim_no_bc,
        pvalue_kernel_eta_optim_yes_bc = pvalue_kernel_eta_optim_yes_bc,
        pvalue_kernel_eta_iqr_no_bc = pvalue_kernel_eta_iqr_no_bc,
        pvalue_kernel_eta_iqr_yes_bc = pvalue_kernel_eta_iqr_yes_bc,
        ci_auc = c(lower_ci_auc, upper_ci_auc),
        ci_youden = c(lower_ci_youden, upper_ci_youden),
        ci_parametric_eta_no_bc = c(lower_ci_parametric_eta_no_bc, upper_ci_parametric_eta_no_bc),
        ci_parametric_eta_yes_bc = c(lower_ci_parametric_eta_yes_bc, upper_ci_parametric_eta_yes_bc),
        ci_kernel_eta_hscv_no_bc = c(lower_ci_kernel_eta_hscv_no_bc, upper_ci_kernel_eta_hscv_no_bc),
        ci_kernel_eta_hscv_yes_bc = c(lower_ci_kernel_eta_hscv_yes_bc, upper_ci_kernel_eta_hscv_yes_bc),
        ci_kernel_eta_optim_no_bc = c(lower_ci_kernel_eta_optim_no_bc, upper_ci_kernel_eta_optim_no_bc),
        ci_kernel_eta_optim_yes_bc = c(lower_ci_kernel_eta_optim_yes_bc, upper_ci_kernel_eta_optim_yes_bc),
        ci_kernel_eta_iqr_no_bc = c(lower_ci_kernel_eta_iqr_no_bc, upper_ci_kernel_eta_iqr_no_bc),
        ci_kernel_eta_iqr_yes_bc = c(lower_ci_kernel_eta_iqr_yes_bc, upper_ci_kernel_eta_iqr_yes_bc)
    ))
}

format_value <- function(x) {
    if (is.null(x)) return("NULL")
    if (is.atomic(x)) return(paste(round(as.numeric(x), 4), collapse = ", "))
    paste(capture.output(print(x)), collapse = " | ")
}


## DATA 207039N
###################################
data_207039N = read.table("/Users/Riki/Desktop/ucm/articulo/Simulciones/Ejemplo Bantis0 compartido/53,53/data207039.txt", header = FALSE, sep = "", dec = ".")
controls_207039N = data_207039N$V1[data_207039N$V2==0]
cases_207039N = data_207039N$V1[data_207039N$V2==1]

auc_207039N_cases_controls = calculate_auc_normal(cases_207039N, controls_207039N)
auc_207039N_controls_cases = calculate_auc_normal(controls_207039N, cases_207039N)
auc_207039N = max(auc_207039N_cases_controls, auc_207039N_controls_cases)

youden_207039N_cases_controls = calculate_youden_normal(cases_207039N, controls_207039N)
youden_207039N_controls_cases = calculate_youden_normal(controls_207039N, cases_207039N)
youden_207039N = max(youden_207039N_cases_controls, youden_207039N_controls_cases)

parametric_eta_no_bc_207039N = parametric_eta(controls_207039N, cases_207039N, 1, box_cox = FALSE)
parametric_eta_yes_bc_207039N = parametric_eta(controls_207039N, cases_207039N, 1, box_cox = TRUE)
kernel_eta_hscv_no_bc_207039N= kernel_eta(controls_207039N, cases_207039N, "hscv", 1, box_cox = FALSE)
kernel_eta_hscv_yes_bc_207039N= kernel_eta(controls_207039N, cases_207039N, "hscv", 1, box_cox = TRUE)
kernel_eta_optim_no_bc_207039N = kernel_eta(controls_207039N, cases_207039N, "optimal", 1, box_cox = FALSE)
kernel_eta_optim_yes_bc_207039N = kernel_eta(controls_207039N, cases_207039N, "optimal", 1, box_cox = TRUE)
kernal_eta_iqr_no_bc_207039N = kernel_eta(controls_207039N, cases_207039N, "iqr", 1, box_cox = FALSE)
kernal_eta_iqr_yes_bc_207039N = kernel_eta(controls_207039N, cases_207039N, "iqr", 1, box_cox = TRUE)


lines <- c(
    "===== RESULTS for 207039N =====",
    sprintf("AUC (cases vs controls): %s", format_value(auc_207039N_cases_controls)),
    sprintf("AUC (controls vs cases): %s", format_value(auc_207039N_controls_cases)),
    sprintf("AUC chosen: %s", format_value(auc_207039N)),
    "",
    sprintf("Youden (cases vs controls): %s", format_value(youden_207039N_cases_controls)),
    sprintf("Youden (controls vs cases): %s", format_value(youden_207039N_controls_cases)),
    sprintf("Youden chosen: %s", format_value(youden_207039N)),
    "",
    "---- ETA estimates ----",
    sprintf("Parametric ETA (no BC): %s", format_value(parametric_eta_no_bc_207039N)),
    sprintf("Parametric ETA (yes BC): %s", format_value(parametric_eta_yes_bc_207039N)),
    sprintf("Kernel ETA (hscv, no BC): %s", format_value(kernel_eta_hscv_no_bc_207039N)),
    sprintf("Kernel ETA (hscv, yes BC): %s", format_value(kernel_eta_hscv_yes_bc_207039N)),
    sprintf("Kernel ETA (optimal, no BC): %s", format_value(kernel_eta_optim_no_bc_207039N)),
    sprintf("Kernel ETA (optimal, yes BC): %s", format_value(kernel_eta_optim_yes_bc_207039N)),
    sprintf("Kernel ETA (iqr, no BC): %s", format_value(kernal_eta_iqr_no_bc_207039N)),
    sprintf("Kernel ETA (iqr, yes BC): %s", format_value(kernal_eta_iqr_yes_bc_207039N)),
    "================================"
)

cat(paste(lines, collapse = "\n"), "\n")
invisible(lines)


p_values_207039N = get_confidence_interval(controls_207039N, cases_207039N, conf_level = 0.95, n_boot = 500)


## DATA 209644N
#########################################
data_209644N = read.table("/Users/Riki/Desktop/ucm/articulo/Simulciones/Ejemplo Bantis0 compartido/53,53/data209644.txt", header = FALSE, sep = "", dec = ".")
controls_209644N = data_209644N$V1[data_209644N$V2==0]
cases_209644N = data_209644N$V1[data_209644N$V2==1]
shapiro.test(controls_209644N)
shapiro.test(cases_209644N)


auc_209644N_cases_controls = calculate_auc_normal(cases_209644N, controls_209644N)
auc_209644N_controls_cases = calculate_auc_normal(controls_209644N, cases_209644N)
auc_209644N = max(auc_209644N_cases_controls, auc_209644N_controls_cases)
youden_209644N_cases_controls = calculate_youden_normal(cases_209644N, controls_209644N)
youden_209644N_controls_cases = calculate_youden_normal(controls_209644N, cases_209644N)
youden_209644N = max(youden_209644N_cases_controls, youden_209644N_controls_cases)

parametric_eta_no_bc_209644N = parametric_eta(controls_209644N, cases_209644N, 1, box_cox = FALSE)
parametric_eta_yes_bc_209644N = parametric_eta(controls_209644N, cases_209644N, 1, box_cox = TRUE)
kernel_eta_hscv_no_bc_209644N= kernel_eta(controls_209644N, cases_209644N, "hscv", 1, box_cox = FALSE)
kernel_eta_hscv_yes_bc_209644N= kernel_eta(controls_209644N, cases_209644N, "hscv", 1, box_cox = TRUE)
kernel_eta_optim_no_bc_209644N = kernel_eta(controls_209644N, cases_209644N, "optimal", 1, box_cox = FALSE)
kernel_eta_optim_yes_bc_209644N = kernel_eta(controls_209644N, cases_209644N, "optimal", 1, box_cox = TRUE)
kernal_eta_iqr_no_bc_209644N = kernel_eta(controls_209644N, cases_209644N, "iqr", 1, box_cox = FALSE)
kernal_eta_iqr_yes_bc_209644N = kernel_eta(controls_209644N, cases_209644N, "iqr", 1, box_cox = TRUE)

lines <- c(
    "===== RESULTS for 209644N =====",
    sprintf("AUC (cases vs controls): %s", format_value(auc_209644N_cases_controls)),
    sprintf("AUC (controls vs cases): %s", format_value(auc_209644N_controls_cases)),
    sprintf("AUC chosen: %s", format_value(auc_209644N)),
    "",
    sprintf("Youden (cases vs controls): %s", format_value(youden_209644N_cases_controls)),
    sprintf("Youden (controls vs cases): %s", format_value(youden_209644N_controls_cases)),
    sprintf("Youden chosen: %s", format_value(youden_209644N)),
    "",
    "---- ETA estimates ----",
    sprintf("Parametric ETA (no BC): %s", format_value(parametric_eta_no_bc_209644N)),
    sprintf("Parametric ETA (yes BC): %s", format_value(parametric_eta_yes_bc_209644N)),
    sprintf("Kernel ETA (hscv, no BC): %s", format_value(kernel_eta_hscv_no_bc_209644N)),
    sprintf("Kernel ETA (hscv, yes BC): %s", format_value(kernel_eta_hscv_yes_bc_209644N)),
    sprintf("Kernel ETA (optimal, no BC): %s", format_value(kernel_eta_optim_no_bc_209644N)),
    sprintf("Kernel ETA (optimal, yes BC): %s", format_value(kernel_eta_optim_yes_bc_209644N)),
    sprintf("Kernel ETA (iqr, no BC): %s", format_value(kernal_eta_iqr_no_bc_209644N)),
    sprintf("Kernel ETA (iqr, yes BC): %s", format_value(kernal_eta_iqr_yes_bc_209644N)),
    "================================"
)

cat(paste(lines, collapse = "\n"), "\n")
invisible(lines)

p_values_209644N = get_confidence_interval(controls_209644N, cases_209644N, conf_level = 0.95, n_boot = 500)