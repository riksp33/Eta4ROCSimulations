import os
import json
from pathlib import Path

def process_json_files(file_paths):
    all_data = []
    
    for file_path in file_paths:
        with open(file_path, 'r') as f:
            data = json.load(f)
            # Extract filename without extension to use as scenario name
            scenario_name = Path(file_path).stem
            data['scenario_name'] = scenario_name
            all_data.append(data)
    
    return all_data

def create_latex_tables(all_data):
    # Extract AUC values and t0 values from the first file to set up the structure
    auc_values = []
    t0_values = []
    
    for auc_key in all_data[0].keys():
        if auc_key.startswith("AUC:"):
            auc_values.append(auc_key)
            # Extract t0 values from the first AUC
            if not t0_values and auc_key != "header" and auc_key != "scenario_name":
                for t0_key in all_data[0][auc_key].keys():
                    if t0_key.startswith("t0:"):
                        t0_values.append(t0_key)
    
    # Sort the values
    auc_values.sort(key=lambda x: float(x.split(":")[1].strip()))
    t0_values.sort(key=lambda x: float(x.split(":")[1].strip()))
    
    # Create bias table
    bias_table = generate_table(all_data, auc_values, t0_values, "bias", "Bias")
    
    # Create RMSE table
    rmse_table = generate_table(all_data, auc_values, t0_values, "rmse", "RMSE")
    
    # Combine tables
    latex_output = bias_table + "\n\n" + rmse_table
    
    return latex_output

def generate_table(all_data, auc_values, t0_values, metric_type, metric_name):
    # Define the methods to include
    methods = ["parametric", "kernel_hscv", "kernel_opt"]
    method_names = {
        "parametric": r"$\hat{\eta}_{log}^{N}$",
        "kernel_hscv": r"$\hat{\eta}_{log}^{K = g, h=csv}$",
        "kernel_opt": r"$\hat{\eta}_{log}^{K = g, h = h*}$"
    }
    
    # Define sample sizes to include
    sample_sizes = ["20", "50", "100"]
    
    # Start building the LaTeX table
    table = r"\begin{table}[H]" + "\n"
    table += r"\centering" + "\n"
    table += r"\caption{Estimación de " + metric_name + r" para diferentes tamaños muestrales}" + "\n"
    table += r"\setlength{\tabcolsep}{5pt} % Ajusta el espacio entre las columnas" + "\n"
    
    # Calculate number of columns for all t0 values and methods
    num_methods = len(methods)
    num_t0_values = len(t0_values)
    total_data_columns = num_methods * num_t0_values
    
    # Define the table format
    table += r"\begin{tabularx}{0.95\textwidth}{c c *{" + str(total_data_columns) + "}{>{\centering\\arraybackslash}X}}" + "\n"
    table += r"\toprule" + "\n"
    
    # First row of headers - t0 values
    table += r"$\eta_{log}$ & n"
    for t0 in t0_values:
        t0_value = t0.split(":")[1].strip()
        for _ in range(num_methods):
            table += f" & $t_0 = {t0_value}$"
    table += r" \\" + "\n"
    
    # Second row of headers - methods
    table += r" & "
    for _ in range(num_t0_values):
        for method in methods:
            table += f" & {method_names[method]}"
    table += r" \\" + "\n"
    
    table += r"\midrule" + "\n"
    
    # Data rows grouped by AUC
    for auc in auc_values:
        auc_value = auc.split(":")[1].strip()
        
        # Get true_eta value for this AUC (using the first t0, assuming it's similar across t0s)
        first_t0 = t0_values[0]
        true_eta = "N/A"
        for data in all_data:
            if auc in data and first_t0 in data[auc]:
                if "true_eta" in data[auc][first_t0]:
                    true_eta = data[auc][first_t0]["true_eta"][0]
                    break
        
        # First row - sample size 20
        table += f" & {sample_sizes[0]}"
        
        # Add data for each t0 and method
        for t0 in t0_values:
            for method in methods:
                # Average over all files for this configuration
                metric_values = []
                for data in all_data:
                    if auc in data and t0 in data[auc]:
                        size_key = f"size: {sample_sizes[0]}"
                        if size_key in data[auc][t0]:
                            if method in data[auc][t0][size_key]:
                                if metric_type in data[auc][t0][size_key][method]:
                                    metric_values.append(data[auc][t0][size_key][method][metric_type][0])
                
                if metric_values:
                    avg_metric = sum(metric_values) / len(metric_values)
                    table += f" & {avg_metric:.4f}"
                else:
                    table += " & -"
        
        table += r" \\" + "\n"
        
        # Second row - sample size 50
        table += f"{true_eta} & {sample_sizes[1]}"
        
        # Add data for each t0 and method
        for t0 in t0_values:
            for method in methods:
                # Average over all files for this configuration
                metric_values = []
                for data in all_data:
                    if auc in data and t0 in data[auc]:
                        size_key = f"size: {sample_sizes[1]}"
                        if size_key in data[auc][t0]:
                            if method in data[auc][t0][size_key]:
                                if metric_type in data[auc][t0][size_key][method]:
                                    metric_values.append(data[auc][t0][size_key][method][metric_type][0])
                
                if metric_values:
                    avg_metric = sum(metric_values) / len(metric_values)
                    table += f" & {avg_metric:.4f}"
                else:
                    table += " & -"
        
        table += r" \\" + "\n"
        
        # Third row - sample size 100
        mu_d_placeholder = f"($\\mu_D = {auc_value}$)"
        table += f"{mu_d_placeholder} & {sample_sizes[2]}"
        
        # Add data for each t0 and method
        for t0 in t0_values:
            for method in methods:
                # Average over all files for this configuration
                metric_values = []
                for data in all_data:
                    if auc in data and t0 in data[auc]:
                        size_key = f"size: {sample_sizes[2]}"
                        if size_key in data[auc][t0]:
                            if method in data[auc][t0][size_key]:
                                if metric_type in data[auc][t0][size_key][method]:
                                    metric_values.append(data[auc][t0][size_key][method][metric_type][0])
                
                if metric_values:
                    avg_metric = sum(metric_values) / len(metric_values)
                    table += f" & {avg_metric:.4f}"
                else:
                    table += " & -"
        
        table += r" \\" + "\n"
        
        # Add midrule between AUC sections
        table += r"\midrule" + "\n"
    
    # Close the table
    table += r"\end{tabularx}" + "\n"
    table += r"\label{table:" + metric_type + "}" + "\n"
    table += r"\end{table}"
    
    return table

def main():
    # Directory where the script is located
    script_dir = Path(__file__).parent.absolute()
    
    # Buscar automáticamente todos los archivos JSON en el directorio
    # json_files = list(script_dir.glob("*.json"))
    
    # Alternativamente, puedes especificar manualmente los archivos
    json_files = [
        script_dir / "normal_1.json",
        script_dir / "normal_2.json",
        script_dir / "normal_3.json",
        # Agrega más archivos según sea necesario
    ]
    
    # Make sure all files exist
    valid_files = []
    for file_path in json_files:
        if file_path.exists():
            valid_files.append(file_path)
        else:
            print(f"Warning: File {file_path} not found")
    
    # Process the files
    if valid_files:
        all_data = process_json_files(valid_files)
        latex_tables = create_latex_tables(all_data)
        
        # Write output to file
        output_file = script_dir / "latex_tables.txt"
        with open(output_file, 'w') as f:
            f.write(latex_tables)
        
        print(f"LaTeX tables have been written to {output_file}")
        print(f"Processed {len(valid_files)} JSON files: {[f.name for f in valid_files]}")
    else:
        print("No valid files to process")

if __name__ == "__main__":
    main()