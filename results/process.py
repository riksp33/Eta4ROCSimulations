#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_tables.py

Genera tablas LaTeX de BIAS y RMSE para simulaciones ROC,
formateadas según el estilo especificado.
"""

import os
import json
from textwrap import indent

# # === CONFIGURACIÓN GLOBAL ===
os.chdir("/Users/Riki/Desktop/ucm/articulo/Simulciones/results/")
SCENARIO_TYPE = 'normal'  # 'normal', 'lognormal' o 'gamma'
NO_BC_FILES = [
    "normal_1.json",
    "normal_2.json",
    "normal_3.json"
    # Agrega más archivos aquí...
]
BC_FILES = [
    "normal_1_box_cox_parametric_and_kernel.json",
    "normal_2_box_cox_parametric_and_kernel.json",
    "normal_3_box_cox_parametric_and_kernel.json",
    # Agrega más archivos aquí...
]
OUTPUT_FILE = "nuevas_tablas_normal.txt"

# # === CONFIGURACIÓN GLOBAL ===
# os.chdir("/Users/Riki/Desktop/ucm/articulo/Simulciones/results/")
# SCENARIO_TYPE = 'lognormal'  # 'normal', 'lognormal' o 'gamma'
# NO_BC_FILES = [
#     "lognormal_1_box_cox_parametric.json",
#     "lognormal_2_box_cox_parametric.json",
#     "lognormal_3_box_cox_parametric.json",
#     "lognormal_4_box_cox_parametric.json",

#     # Agrega más archivos aquí...
# ]
# BC_FILES = [
#     "lognormal_1_box_cox_parametric_and_kernel.json",
#     "lognormal_2_box_cox_parametric_and_kernel.json",
#     "lognormal_3_box_cox_parametric_and_kernel.json",
#     "lognormal_4_box_cox_parametric_and_kernel.json",

# ]
# OUTPUT_FILE = "nuevas_tablas_lognormal.txt"

# === CONFIGURACIÓN GLOBAL ===
# os.chdir("/Users/Riki/Desktop/ucm/articulo/Simulciones/results/")
# SCENARIO_TYPE = 'gamma'  # 'normal', 'lognormal' o 'gamma'
# NO_BC_FILES = [
#     "gamma_1_box_cox_parametric.json",
#     "gamma_2_box_cox_parametric.json",
#     "gamma_3_box_cox_parametric.json",

#     # Agrega más archivos aquí...
# ]
# BC_FILES = [
#     "gamma_1_box_cox_parametric_and_kernel.json",
#     "gamma_2_box_cox_parametric_and_kernel.json",
#     "gamma_3_box_cox_parametric_and_kernel.json",

# ]
# OUTPUT_FILE = "nuevas_tablas_gamma.txt"

def process_files(no_bc_files, bc_files, scenario_type, output_file):
    """
    Procesa archivos JSON para crear tablas LaTeX de BIAS y RMSE
    """
    if len(no_bc_files) != len(bc_files) and len(bc_files) > 0:
        raise ValueError("Las listas de archivos deben tener la misma longitud")
    
    all_tables = []
    
    for i in range(len(no_bc_files)):
        no_bc_file = no_bc_files[i]
        bc_file = bc_files[i] if bc_files else None
        
        # Cargar datos
        with open(no_bc_file, 'r') as f:
            no_bc_data = json.load(f)
            
        bc_data = None
        if bc_file:
            with open(bc_file, 'r') as f:
                bc_data = json.load(f)
        
        # Extraer información de archivo
        file_id = os.path.splitext(os.path.basename(no_bc_file))[0]
        
        # Limpiar file_id para extraer solo el identificador base
        if scenario_type == 'normal':
            base_id = file_id
        else:
            parts = file_id.split('_')
            base_id = f"{parts[0]}_{parts[1]}"
        
        # Obtener valores de AUC y t0
        auc_values = ["0.6", "0.75", "0.9"]
        t0_values = sorted({
            t0 for auc in no_bc_data if auc.startswith("AUC") for t0 in no_bc_data[auc] if t0.startswith("t0")
        }, key=lambda x: float(x.split(": ")[1]))
        
        # Crear tablas
        latex_bias = create_latex_table(no_bc_data, bc_data, auc_values, t0_values, "BIAS", scenario_type, base_id)
        latex_rmse = create_latex_table(no_bc_data, bc_data, auc_values, t0_values, "RMSE", scenario_type, base_id)
        
        all_tables.extend([latex_bias, "", latex_rmse, ""])
    
    with open(output_file, 'w') as f:
        f.write("\n".join(all_tables))
    
    print(f"Todas las tablas fueron escritas en {output_file}")

def create_latex_table(no_bc_data, bc_data, auc_values, t0_values, metric_type, scenario_type, file_id):
    """
    Crea una tabla LaTeX para un tipo de métrica específica
    """
    sizes = ["50", "100"]
    metric_type_lower = metric_type.lower()
    
    if scenario_type == 'normal':
        estimadores = [
            "parametric", "parametric_bc", "kernel_hscv", "kernel_hscv_bc",
            "kernel_opt", "kernel_opt_bc", "kernel_iqr", "kernel_iqr_bc"
        ]
        nombres_est = {
            "parametric": "$\\widehat{\\eta}_{p}^{N}$",
            "parametric_bc": "$\\widehat{\\eta}_{p}^{N_{T}}$",
            "kernel_hscv": "$\\widehat{\\eta}_{p}^{K_{HSCV}}$",
            "kernel_hscv_bc": "$\\widehat{\\eta}_{p}^{K_{HSCV_{T}}}$",
            "kernel_opt": "$\\widehat{\\eta}_{p}^{K_{OPT}}$",
            "kernel_opt_bc": "$\\widehat{\\eta}_{p}^{K_{OPT_{T}}}$",
            "kernel_iqr": "$\\widehat{\\eta}_{p}^{K_{IQR}}$",
            "kernel_iqr_bc": "$\\widehat{\\eta}_{p}^{K_{IQR_{T}}}$"
        }
    else:
        estimadores = [
            "parametric_bc", "kernel_hscv", "kernel_hscv_bc",
            "kernel_opt", "kernel_opt_bc", "kernel_iqr", "kernel_iqr_bc"
        ]
        nombres_est = {
            "parametric_bc": "$\\widehat{\\eta}_{p}^{N_{T}}$",
            "kernel_hscv": "$\\widehat{\\eta}_{p}^{K_{HSCV}}$",
            "kernel_hscv_bc": "$\\widehat{\\eta}_{p}^{K_{HSCV_{T}}}$",
            "kernel_opt": "$\\widehat{\\eta}_{p}^{K_{OPT}}$",
            "kernel_opt_bc": "$\\widehat{\\eta}_{p}^{K_{OPT_{T}}}$",
            "kernel_iqr": "$\\widehat{\\eta}_{p}^{K_{IQR}}$",
            "kernel_iqr_bc": "$\\widehat{\\eta}_{p}^{K_{IQR_{T}}}$"

        }
    
    dist_names = {
        'normal': '\\text{N}',
        'lognormal': '\\text{LogN}',
        'gamma': '\\text{Gamma}'
    }
    
    if scenario_type == 'normal':
        if file_id == 'normal_1':
            dist_title = "\\( X\\sim\\text{N}(0, 1), \\, Y\\sim\\text{N}(\\mu, 1) \\)"
        elif file_id == 'normal_2':
            dist_title = "\\( X\\sim\\text{N}(0, 1), \\, Y\\sim\\text{N}(\\mu, \\sigma^2) \\)"
        elif file_id == 'normal_3':
            dist_title = "\\( X\\sim\\text{N}(\\mu_X, \\sigma_X^2), \\, Y\\sim\\text{N}(\\mu_Y, \\sigma_Y^2) \\)"
        else:
            dist_title = f"\\( {file_id} \\)"
    elif scenario_type == 'lognormal':
        if file_id == 'lognormal_1':
            dist_title = "\\( X\\sim\\text{LogN}(0, 0.25), \\, Y\\sim\\text{LogN}(\\mu, 0.25) \\)"
        elif file_id == 'lognormal_2':
            dist_title = "\\( X\\sim\\text{LogN}(0, 0.25), \\, Y\\sim\\text{LogN}(\\mu, \\sigma^2) \\)"
        elif file_id == 'lognormal_3':
            dist_title = "\\( X\\sim\\text{LogN}(0, 1), \\, Y\\sim\\text{LogN}(\\mu, 1) \\)"
        elif file_id == 'lognormal_4':
            dist_title = "\\( X\\sim\\text{LogN}(0, 1), \\, Y\\sim\\text{LogN}(\\mu, \\sigma^2) \\)"
        else:
            dist_title = f"\\( {file_id} \\)"
    elif scenario_type == 'gamma':
        if file_id == 'gamma_1':
            dist_title = "\\( X\\sim\\text{Gamma}(2, 1), \\, Y\\sim\\text{Gamma}(\\alpha, 1) \\)"
        elif file_id == 'gamma_2':
            dist_title = "\\( X\\sim\\text{Gamma}(2, 1), \\, Y\\sim\\text{Gamma}(\\alpha, \\beta) \\)"
        elif file_id == 'gamma_3':
            dist_title = "\\( X\\sim\\text{Gamma}(1, 2), \\, Y\\sim\\text{Gamma}(\\alpha, \\beta) \\)"
        else:
            dist_title = f"\\( {file_id} \\)"
    
    tabla = []
    tabla.append("\\begin{table}[H]")
    tabla.append("\\centering")
    tabla.append(f"\\caption{{{metric_type}: {dist_title}}}")
    tabla.append("\\vspace{1em}")
    tabla.append(f"\\label{{tabla_{metric_type.lower()}_{file_id}}}")
    tabla.append("\\renewcommand{\\arraystretch}{1.2}")
    tabla.append("\\resizebox{\\textwidth}{!}{%")
    tabla.append("\\begin{tabular}{c c *{3}{c} c *{3}{c}}")
    tabla.append("\\toprule")
    tabla.append("\\multicolumn{2}{c}{} & \\multicolumn{3}{c}{\\textbf{n = m = 50}} & & \\multicolumn{3}{c}{\\textbf{n = m = 100}} \\\\")
    tabla.append("\\cmidrule(lr){3-5} \\cmidrule(lr){7-9}")
    tabla.append("\\multicolumn{2}{c}{} & " + " & ".join([f"\\textbf{{AUC = {a}}}" for a in auc_values]) +
                " & & " + " & ".join([f"\\textbf{{AUC = {a}}}" for a in auc_values]) + " \\\\")
    tabla.append("\\midrule")
    
    for t0 in t0_values:
        t0_val = t0.split(": ")[1]
        tabla.append("")
        tabla.append(f"\\textbf{{$t_0 = {t0_val}$}} ")
        
        for est in estimadores:
            if not est.endswith("_bc"):
                if est == estimadores[0] or est.endswith("_bc"):
                    row_prefix = "&"
                else:
                    row_prefix = " &"

                tabla.append(f"{row_prefix} {nombres_est[est]} & "  +
                           " & ".join([
                               format_value(no_bc_data, f"AUC: {auc}", t0, "50", est, metric_type_lower)
                               for auc in auc_values
                           ]) +
                           " &  & " +
                           " & ".join([
                               format_value(no_bc_data, f"AUC: {auc}", t0, "100", est, metric_type_lower)
                               for auc in auc_values
                           ]) + " \\\\")
            else:
                if bc_data:
                    base_est = est[:-3]
                    tabla.append(f" & {nombres_est[est]} & " +
                               " & ".join([
                                   format_value(bc_data, f"AUC: {auc}", t0, "50", base_est, metric_type_lower)
                                   for auc in auc_values
                               ]) +
                               " &  & " +
                               " & ".join([
                                   format_value(bc_data, f"AUC: {auc}", t0, "100", base_est, metric_type_lower)
                                   for auc in auc_values
                               ]) + " \\\\")
                else:
                    tabla.append(f" & {nombres_est[est]} & --- & --- & --- &  & --- & --- & --- \\\\")
        
        tabla.append("\\midrule")
    
    tabla.append("")
    tabla.append("\\bottomrule")
    tabla.append("\\end{tabular}")
    tabla.append("}")
    tabla.append("\\end{table}")
    
    return "\n".join(tabla)

def format_value(data, auc_key, t0, size, est, metric_type):
    """
    Formatea un valor para la tabla
    """
    try:
        val = data[auc_key][t0][f"size: {size}"][est][metric_type][0]
        if isinstance(val, str):
            return "-NaN-"
        return f"{val:.4f}"
    except KeyError:
        return "---"

if __name__ == "__main__":
    process_files(NO_BC_FILES, BC_FILES, SCENARIO_TYPE, OUTPUT_FILE)
