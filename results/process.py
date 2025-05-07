import json
import os

def process_multiple_json_files(json_files, output_file):
    all_tables = []

    for json_file in json_files:
        print(json_file)
        with open(json_file, 'r') as f:
            data = json.load(f)

        file_id = os.path.splitext(os.path.basename(json_file))[0]
        auc_values = ["0.6", "0.75", "0.9"]
        t0_values = sorted({
            t0 for auc in data if auc.startswith("AUC") for t0 in data[auc] if t0.startswith("t0")
        }, key=lambda x: float(x.split(": ")[1]))

        latex_bias = create_latex_table(data, auc_values, t0_values, "bias", file_id)
        latex_rmse = create_latex_table(data, auc_values, t0_values, "rmse", file_id)

        all_tables.extend([latex_bias, "", latex_rmse, ""])

    with open(output_file, 'w') as f:
        f.write("\n".join(all_tables))

    print(f"Todas las tablas fueron escritas en {output_file}")

def create_latex_table(data, auc_values, t0_values, metric_type, file_id):
    sizes = ["50", "100"]
    estimadores = ["parametric", "kernel_hscv", "kernel_opt", "kernel_iqr"]
    nombres_est = {
        "parametric": "$\\widehat{\\eta}_{p}^{N}$",
        "kernel_hscv": "$\\widehat{\\eta}_{p}^{K_{hscv}}$",
        "kernel_opt": "$\\widehat{\\eta}_{p}^{K_{h^*}}$",
        "kernel_iqr": "$\\widehat{\\eta}_{p}^{K_{IQR}}$"
    }

    tabla = []
    tabla.append("\\begin{table}[H]")
    tabla.append("\\centering")
    tabla.append(f"\\caption{{Resultados de \\texttt{{{file_id}}} ({metric_type.upper()})}}")
    tabla.append("\\vspace{1em}")
    tabla.append(f"\\label{{tabla_{metric_type}_{file_id}}}")
    tabla.append("\\renewcommand{\\arraystretch}{1.2}")
    tabla.append("\\resizebox{\\textwidth}{!}{%")
    tabla.append("\\begin{tabular}{c c *{3}{c} c *{3}{c}}")
    tabla.append("\\toprule")
    tabla.append("\\multicolumn{2}{c}{} & \\multicolumn{3}{c}{\\textbf{n = m = 50}} & & \\multicolumn{3}{c}{\\textbf{n = m = 100}} \\\\")
    tabla.append("\\cmidrule(lr){3-5} \\cmidrule(lr){7-9}")
    tabla.append("\\multicolumn{2}{c}{} & " + " & ".join([f"\\textbf{{AUC = {a}}}" for a in auc_values[:3]]) +
                 " & & " + " & ".join([f"\\textbf{{AUC = {a}}}" for a in auc_values[:3]]) + " \\\\")
    tabla.append("\\midrule")

    for t0 in t0_values:
        t0_val = t0.split(": ")[1]
        tabla.append(f"\\textbf{{$t_0 = {t0_val}$}} & {nombres_est[estimadores[0]]}" + " & " +
                     " & ".join([
                         format_value(data, f"AUC: {auc}", t0, size, estimadores[0], metric_type)
                         for size in sizes for auc in auc_values
                     ][:3] + [""] + [
                         format_value(data, f"AUC: {auc}", t0, size, estimadores[0], metric_type)
                         for size in sizes for auc in auc_values
                     ][3:]) + " \\\\")
        for est in estimadores[1:]:
            tabla.append(f" & {nombres_est[est]}" + " & " +
                         " & ".join([
                             format_value(data, f"AUC: {auc}", t0, size, est, metric_type)
                             for size in sizes for auc in auc_values
                         ][:3] + [""] + [
                             format_value(data, f"AUC: {auc}", t0, size, est, metric_type)
                             for size in sizes for auc in auc_values
                         ][3:]) + " \\\\")
        tabla.append("\\midrule")

    tabla.append("\\bottomrule")
    tabla.append("\\end{tabular}")
    tabla.append("}")
    tabla.append("\\end{table}")

    return "\n".join(tabla)

def format_value(data, auc_key, t0, size, est, metric_type):
    try:
        val = data[auc_key][t0][f"size: {size}"][est][metric_type][0]
        if type(val) == str:
            print( auc_key, t0, size, est, metric_type)
            return "-NaN-"

        return f"{val:.4f}"
    except KeyError:
        return "-"

def main():
    os.chdir("/Users/Riki/Desktop/ucm/articulo/Simulciones/results/")

    files = {
        "normales" : ["normal_1.json", "normal_2.json", "normal_3.json"],
        "lognormales_no_bc" : ["lognormal_1_box_cox_parametric.json",
                               "lognormal_2_box_cox_parametric.json",
                               "lognormal_3_box_cox_parametric.json",
                               "lognormal_4_box_cox_parametric.json"],

        "lognormales_si_bc" : [ "lognormal_1_box_cox_parametric_and_kernel.json",
                                "lognormal_2_box_cox_parametric_and_kernel.json",
                                "lognormal_3_box_cox_parametric_and_kernel.json",
                                "lognormal_4_box_cox_parametric_and_kernel.json"],

        "gamma_no_bc" :       ["gamma_1_box_cox_parametric.json",
                               "gamma_2_box_cox_parametric.json",
                               "gamma_3_box_cox_parametric.json"],

        "gamma_si_bc" :      [  "gamma_1_box_cox_parametric_and_kernel.json",
                                "gamma_2_box_cox_parametric_and_kernel.json",
                                "gamma_3_box_cox_parametric_and_kernel.json"]
        
    }

    for key, value in files.items():
        print(key)
        process_multiple_json_files(value, key+".txt")




    # json_files = ["normal_1.json", "normal_2.json", "normal_3.json"]  # Aquí defines tus archivos
    # output_file = "metrics_tables.txt"
    # process_multiple_json_files(json_files, output_file)

if __name__ == "__main__":
    main()
