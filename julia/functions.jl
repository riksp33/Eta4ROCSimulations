
using Optim
using SpecialFunctions: erfc, erfcinv
using Distributions
using Statistics
using StatsBase
using JSON
using Random
using LinearAlgebra

# Helpers for JSON diagnostics
function stringify_keys(x)
    if isa(x, Dict)
        return Dict(string(k) => stringify_keys(v) for (k,v) in x)
    elseif isa(x, AbstractVector)
        return [stringify_keys(v) for v in x]
    else
        return x
    end
end

function json_type_map(x)
    if isa(x, Dict)
        return Dict(string(k) => json_type_map(v) for (k,v) in x)
    elseif isa(x, AbstractVector)
        return [json_type_map(v) for v in x]
    else
        return string(typeof(x))
    end
end

# Callable struct for likelihood function. This avoids allocations when calling optimizer
struct LikBox!
    log_x::Vector{Float64}
    log_y::Vector{Float64}
    xh::Vector{Float64}
    yh::Vector{Float64}
    n::Int
    m::Int
    inv_n::Float64
    inv_m::Float64
    n_half::Float64
    m_half::Float64
    sum_log_xy::Float64
end

function (L::LikBox!)(h::Float64)
    log_x = L.log_x
    log_y = L.log_y
    xh    = L.xh
    yh    = L.yh
    n     = L.n
    m     = L.m

    # Case where h = 0. The transformation is log. 
    if abs(h) < 1e-5
        @inbounds  for i in 1:n
            xh[i] = log_x[i]
        end
        @inbounds  for i in 1:m
            yh[i] = log_y[i]
        end
    # General case, transformation is (exp(h*log(x)) -1)/h. 
    else
        h_inv = 1.0 / h
        @inbounds  for i in 1:n
            xh[i] = (exp(h * log_x[i]) - 1.0) * h_inv
        end
        @inbounds  for i in 1:m
            yh[i] = (exp(h * log_y[i]) - 1.0) * h_inv
        end
    end

    # mean x
    sum_x = 0.0
    @inbounds  for i in 1:n
        sum_x += xh[i]
    end
    mean_x = sum_x * L.inv_n

    # mean y
    sum_y = 0.0
    @inbounds  for i in 1:m
        sum_y += yh[i]
    end
    mean_y = sum_y * L.inv_m

    # vaariance x
    ssx = 0.0
    @inbounds  for i in 1:n
        d = xh[i] - mean_x
        ssx += d * d
    end

    # variance y
    ssy = 0.0
    @inbounds  for i in 1:m
        d = yh[i] - mean_y
        ssy += d * d
    end

    ll = L.n_half * log(ssx * L.inv_n) +
         L.m_half * log(ssy * L.inv_m) +
         (h - 1.0) * L.sum_log_xy

    return -ll
end


function apply_box_cox(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    m = length(y)

    # Make all values positive if they are slightly negative or zero
    min_x = minimum(x)
    min_y = minimum(y)

    offset = (min(min_x, min_y) <= 0.0) ?
             (-min(min_x, min_y) + 5e-4) :
             0.0

    log_x = Vector{Float64}(undef, n)
    log_y = Vector{Float64}(undef, m)

    # Apply log with offset at the same time. This avoids allocations and branching
    @inbounds  for i in 1:n
        log_x[i] = log(x[i] + offset)
    end
    @inbounds  for i in 1:m
        log_y[i] = log(y[i] + offset)
    end

    sum_log_xy = sum(log_x) + sum(log_y)

    # To be filled inside likbox funtion
    xh = Vector{Float64}(undef, n)
    yh = Vector{Float64}(undef, m)

    lik = LikBox!(
        log_x, log_y, xh, yh,
        n, m,
        1.0 / n,
        1.0 / m,
        -0.5 * n,
        -0.5 * m,
        sum_log_xy
    )

    res = optimize(lik, -1.0, 1.0, Brent())
    h_optim = Optim.minimizer(res)

    if abs(h_optim) < 1e-5
        @inbounds  for i in 1:n
            xh[i] = log_x[i]
        end
        @inbounds  for i in 1:m
            yh[i] = log_y[i]
        end
    else
        h_inv = 1.0 / h_optim
        @inbounds  for i in 1:n
            xh[i] = (exp(h_optim * log_x[i]) - 1.0) * h_inv
        end
        @inbounds  for i in 1:m
            yh[i] = (exp(h_optim * log_y[i]) - 1.0) * h_inv
        end
    end

    # finally return the transformed data 
    return (
        transformed_x = copy(xh),
        transformed_y = copy(yh)
    )
end



function calculate_auc_normal(controls::Vector{Float64}, cases::Vector{Float64})
    n_ctrl = length(controls)
    n_case = length(cases)
    
    # Single-pass mean and variance calculation for controls
    sum_ctrl = 0.0
    sum_sq_ctrl = 0.0
    @inbounds for i in 1:n_ctrl
        val = controls[i]
        sum_ctrl += val
        sum_sq_ctrl += val * val
    end
    mean_ctrl = sum_ctrl / n_ctrl
    var_ctrl = (sum_sq_ctrl - n_ctrl * mean_ctrl * mean_ctrl) / (n_ctrl - 1)
    
    # Single-pass mean and variance calculation for cases
    sum_case = 0.0
    sum_sq_case = 0.0
    @inbounds @simd for i in 1:n_case
        val = cases[i]
        sum_case += val
        sum_sq_case += val * val
    end
    mean_case = sum_case / n_case
    var_case = (sum_sq_case - n_case * mean_case * mean_case) / (n_case - 1)
    
    # AUC calculation
    nAUC = cdf(Normal(0,1),((mean_case - mean_ctrl) / sqrt(var_ctrl + var_case)))
    
    return nAUC
end

function calculate_youden_normal(controls::Vector{Float64}, cases::Vector{Float64})
    n_ctrl = length(controls)
    n_case = length(cases)
    
    # Single-pass mean and variance calculation for controls
    sum_ctrl = 0.0
    sum_sq_ctrl = 0.0
    @inbounds @simd for i in 1:n_ctrl
        val = controls[i]
        sum_ctrl += val
        sum_sq_ctrl += val * val
    end
    mean_ctrl = sum_ctrl / n_ctrl
    var_ctrl = (sum_sq_ctrl - n_ctrl * mean_ctrl * mean_ctrl) / (n_ctrl - 1)
    sd_ctrl = sqrt(var_ctrl)
    
    # Single-pass mean and variance calculation for cases
    sum_case = 0.0
    sum_sq_case = 0.0
    @inbounds @simd for i in 1:n_case
        val = cases[i]
        sum_case += val
        sum_sq_case += val * val
    end
    mean_case = sum_case / n_case
    var_case = (sum_sq_case - n_case * mean_case * mean_case) / (n_case - 1)
    sd_case = sqrt(var_case)
    
    # Pre-compute repeated terms
    mean_diff = mean_ctrl - mean_case
    var_diff = var_ctrl - var_case
    
   if abs(var_diff) < 1e-16
        return 0.0
    end

    # Calculate cstar
    term1 = (mean_case * var_ctrl - mean_ctrl * var_case)
    term2 = sd_ctrl * sd_case * sqrt(mean_diff * mean_diff + var_diff * log(var_ctrl / var_case))
    cstar = (term1 - term2) / var_diff
    
    # Calculate Youden index
    youden = cdf(Normal(0,1), (cstar - mean_ctrl) / sd_ctrl) - cdf(Normal(0,1), (cstar - mean_case) / sd_case)
    
    return youden
end


@inline function standarization(eta::Float64)
    return log(eta + 1.0) / (1 + log(eta + 1.0))
end

function eta_from_roc_curves(
    roc::Vector{Float64},
    roc_dx::Vector{Float64},
    t0::Float64,
    mesh::Vector{Float64})
    n = length(mesh)
    
    # Find limit index (using searchsortedlast for O(log n) instead of O(n))
    lim = searchsortedlast(mesh, t0)
    
    # Early return if t0 is before first mesh point
    if lim == 0
        return standarization(0.0)
    end
    
    # Compute first eta value
    roc_dx_1 = roc_dx[1]
    diff_1_sq = (roc_dx_1 - 1.0) * (roc_dx_1 - 1.0)
    
    eta_sum = if roc_dx_1 <= 1.0
        diff_1_sq * mesh[1]
    else
        diff_1_sq / roc_dx_1 * roc[1]
    end
    
    # Accumulate remaining eta values up to lim (fused loop)
    @inbounds for i in 2:lim
        roc_dx_i = roc_dx[i]
        diff_i = roc_dx_i - 1.0
        diff_i_sq = diff_i * diff_i
        
        eta_i = if roc_dx_i <= 1.0
            diff_i_sq * (mesh[i] - mesh[i-1])
        else
            diff_i_sq / roc_dx_i * (roc[i] - roc[i-1])
        end
        
        eta_sum += eta_i
    end
    
    return standarization(eta_sum)
end

function parametric_eta(controls::Vector{Float64}, cases::Vector{Float64}, doBoxCox::Bool, t0::Float64, mesh::Vector{Float64})
    if doBoxCox
        transformed_controls, transformed_cases = apply_box_cox(controls, cases)
    else
        transformed_controls = controls
        transformed_cases = cases
    end

    # get statistics 
    n_ctrl = length(transformed_controls)
    n_case = length(transformed_cases)

    # Single-pass mean and variance calculation for controls
    sum_ctrl = 0.0
    sum_sq_ctrl = 0.0
    @inbounds for i in 1:n_ctrl
        val = transformed_controls[i]
        sum_ctrl += val
        sum_sq_ctrl += val * val
    end
    mean_ctrl = sum_ctrl / n_ctrl
    var_ctrl = (sum_sq_ctrl - n_ctrl * mean_ctrl * mean_ctrl) / (n_ctrl - 1)
    sd_controls = sqrt(var_ctrl)
    
    # Single-pass mean and variance calculation for cases
    sum_case = 0.0
    sum_sq_case = 0.0
    @inbounds for i in 1:n_case
        val = transformed_cases[i]
        sum_case += val
        sum_sq_case += val * val
    end
    mean_case = sum_case / n_case
    var_case = (sum_sq_case - n_case * mean_case * mean_case) / (n_case - 1)
    sd_cases = sqrt(var_case)

    rho = sd_controls / sd_cases
    delta = (mean_case - mean_ctrl) / sqrt(var_case)

    # Fill ROC and ROCPrima vectors

    length_mesh = length(mesh)
    ROC = Vector{Float64}(undef, length_mesh)
    ROCPrima = Vector{Float64}(undef, length_mesh)

    # Main computation loop
    @inbounds for i in 1:length_mesh
        p = mesh[i]

        qnorm_1mp = quantile(Normal(0,1), 1.0 - p)
        qnorm_p = quantile(Normal(0,1), p)


        ROC[i] = 1 - cdf(Normal(mean_ctrl + delta * sd_controls / rho , sd_controls / rho), quantile(Normal(mean_ctrl, sd_controls), 1 -p))
        ROCPrima[i] = (rho * exp(-0.5 * (delta + rho *quantile(Normal(0,1), p))^2)) / exp(-0.5 * quantile(Normal(0,1), p)^2)
    end
    return eta_from_roc_curves(ROC, ROCPrima, t0, mesh)
end

function kernel_density_estimation(
    data::Vector{Float64},
    mesh::Vector{Float64},
    bandwidth::Float64
)
    n_data = length(data)
    mesh_length = length(mesh)
    
    # Pre-compute constants
    inv_n_bw = 1.0 / (n_data * bandwidth)
    inv_bw = 1.0 / bandwidth
    
    # Pre-create standard normal distribution (reused for pdf calls)
    std_normal = Normal(0.0, 1.0)
    
    # Allocate output
    estimated_density = Vector{Float64}(undef, mesh_length)
    
    # Compute density at each point
    @inbounds for i in 1:mesh_length
        point = mesh[i]
        density_sum = 0.0
        
        # Sum kernel contributions from all data points
        @simd for j in 1:n_data
            z = (point - data[j]) * inv_bw
            density_sum += pdf(std_normal, z)
        end
        
        estimated_density[i] = density_sum * inv_n_bw
    end
    
    return estimated_density
end


function kernel_distribution_estimation(
    data::Vector{Float64},
    points::Vector{Float64},
    bandwidth::Float64
)
    n_data = length(data)
    n_points = length(points)
    
    # Pre-compute constants
    inv_n = 1.0 / n_data
    inv_bw = 1.0 / bandwidth
    
    # Pre-create standard normal distribution (reused for cdf calls)
    std_normal = Normal(0.0, 1.0)
    
    # Allocate output
    estimated_distribution = Vector{Float64}(undef, n_points)
    
    # Compute CDF at each point
    @inbounds for i in 1:n_points
        point = points[i]
        cdf_sum = 0.0
        
        # Sum CDF contributions from all data points
        @simd for j in 1:n_data
            z = (point - data[j]) * inv_bw
            cdf_sum += cdf(std_normal, z)
        end
        
        estimated_distribution[i] = cdf_sum * inv_n
    end
    
    return estimated_distribution
end


function kernel_distribution_estimation2(
    data::Vector{Float64},
    points::Vector{Float64},
    bandwidth::Float64
)
    n_data = length(data)
    n_points = length(points)
    
    # Pre-compute constants
    inv_n = 1.0 / n_data
    inv_bw = 1.0 / bandwidth
    
    # Pre-create standard normal distribution (reused for cdf calls)
    
    
    # Allocate output
    estimated_distribution = Vector{Float64}(undef, n_points)
    
    # Compute CDF at each point
    @inbounds for i in 1:n_points
        point = points[i]
        cdf_sum = 0.0
        
        # Sum CDF contributions from all data points
        @simd for j in 1:n_data
            z = (point - data[j]) * inv_bw
            cdf_sum += cdf(Normal(0, 1.0), z)
        end
        
        estimated_distribution[i] = cdf_sum * inv_n
    end
    
    return estimated_distribution
end



function evaluate_kernel_estimation(
    point::Float64,
    function_values::Vector{Float64},
    mesh::Vector{Float64}
)
    n_points = length(mesh)
    
    # Find position using binary search (O(log n) instead of O(n))
    position = searchsortedlast(mesh, point)
    
    # Determine interpolation indices
    if position == 0
        return function_values[1]
    elseif position >= n_points
        return function_values[n_points]
    else
        return 0.5 * (function_values[position] + function_values[position + 1])
    end
end

function inverse_kernel_estimation(
    point::Float64,
    function_values::Vector{Float64},
    mesh::Vector{Float64}
)
    n_points = length(mesh)
    
    # Find position where function_values crosses point
    position = searchsortedlast(function_values, point)
    
    # Determine interpolation indices
    if position == 0
        return mesh[1]
    elseif position >= n_points
        return mesh[n_points]
    else
        return 0.5 * (mesh[position] + mesh[position + 1])
    end
end



"""
    golden_section_search(f, a, b; tol=1e-5, max_iter=100)

Golden section search for univariate optimization.

# Arguments
- `f`: Objective function to minimize
- `a`: Lower bound
- `b`: Upper bound
- `tol`: Tolerance for convergence (default: 1e-5)
- `max_iter`: Maximum iterations (default: 100)

# Returns
- Optimal point x* that minimizes f(x) in [a, b]
"""
function golden_section_search(
    f::Function,
    a::Float64,
    b::Float64;
    tol::Float64=1e-5,
    max_iter::Int=100
)
    phi = 0.5 * (1.0 + sqrt(5.0))  # Golden ratio
    resphi = 2.0 - phi
    
    # Initial points
    x1 = a + resphi * (b - a)
    x2 = b - resphi * (b - a)
    f1 = f(x1)
    f2 = f(x2)
    
    iter = 0
    while abs(b - a) > tol && iter < max_iter
        if f1 < f2
            b = x2
            x2 = x1
            f2 = f1
            x1 = a + resphi * (b - a)
            f1 = f(x1)
        else
            a = x1
            x1 = x2
            f1 = f2
            x2 = b - resphi * (b - a)
            f2 = f(x2)
        end
        iter += 1
    end
    
    return 0.5 * (a + b)
end

"""
    hscv_optimized(data::Vector{Float64})

Optimized HSCV bandwidth selector using vectorized operations.
Significantly faster than the loop-based version for large datasets.
"""
function hscv_optimized(data::Vector{Float64})
    n = length(data)
    
    # Initial bandwidth using Silverman's rule
    sigma = std(data)
    h_init = 1.06 * sigma * n^(-0.2)
    
    # Define search range
    h_min = h_init * 0.1
    h_max = h_init * 3.0
    
    # Pre-compute pairwise differences (memory intensive for large n)
    # For n > 5000, consider using the loop version
    diff_matrix = data .- data'  # Broadcasting creates n×n matrix
    diff_sq = diff_matrix .^ 2
    
    # Constants
    inv_sqrt_2pi = 0.3989422804014327  # 1/sqrt(2π)
    inv_sqrt_4pi = 0.28209479177387814  # 1/sqrt(4π)
    
    function lscv_objective_optimized(h::Float64)
        if h <= 0.0
            return Inf
        end
        
        # Term 1: R(f̂) using convolution
        h_conv = h * 1.4142135623730951  # h * sqrt(2)
        kernel_conv = exp.(-0.5 .* diff_sq ./ (h_conv^2)) ./ (h_conv * 2.5066282746310002)
        term1 = sum(kernel_conv) / (n * n)
        
        # Term 2: Leave-one-out cross-validation
        # Set diagonal to zero (leave-one-out)
        kernel_loo = exp.(-0.5 .* diff_sq ./ (h^2)) ./ (h * 2.5066282746310002)
        kernel_loo[diagind(kernel_loo)] .= 0.0  # Exclude self-contribution
        term2 = -2.0 * sum(kernel_loo) / (n * (n - 1))
        
        return term1 + term2
    end
    
    # Optimize
    h_opt = golden_section_search(lscv_objective_optimized, h_min, h_max)
    
    return h_opt
end

"""
    hscv_hybrid(data::Vector{Float64}; chunk_size::Int=1000)

Hybrid HSCV implementation that balances speed and memory usage.
Uses chunked computation for large datasets.
"""
function hscv_hybrid(data::Vector{Float64}; chunk_size::Int=1000)
    n = length(data)
    
    # Use optimized version for small datasets
    if n <= chunk_size
        return hscv_optimized(data)
    end
    
    # For large datasets, use chunked approach
    sigma = std(data)
    h_init = 1.06 * sigma * n^(-0.2)
    
    h_min = h_init * 0.1
    h_max = h_init * 3.0
    
    inv_sqrt_2pi = 0.3989422804014327
    sqrt_2 = 1.4142135623730951
    
    function lscv_objective_chunked(h::Float64)
        if h <= 0.0
            return Inf
        end
        
        h_conv = h * sqrt_2
        inv_h_conv_sqrt_2pi = 1.0 / (h_conv * 2.5066282746310002)
        inv_h_sqrt_2pi = 1.0 / (h * 2.5066282746310002)
        
        term1 = 0.0
        term2 = 0.0
        
        # Process in chunks to control memory
        n_chunks = ceil(Int, n / chunk_size)
        
        for chunk_i in 1:n_chunks
            start_i = (chunk_i - 1) * chunk_size + 1
            end_i = min(chunk_i * chunk_size, n)
            
            @inbounds for i in start_i:end_i
                xi = data[i]
                
                # Term 1 contribution
                @simd for j in 1:n
                    diff = xi - data[j]
                    term1 += exp(-0.5 * (diff / h_conv)^2) * inv_h_conv_sqrt_2pi
                end
                
                # Term 2 contribution (leave-one-out)
                sum_loo = 0.0
                @simd for j in 1:n
                    if i != j
                        diff = xi - data[j]
                        sum_loo += exp(-0.5 * (diff / h)^2)
                    end
                end
                term2 += sum_loo * inv_h_sqrt_2pi
            end
        end
        
        term1 = term1 / (n * n)
        term2 = -2.0 * term2 / (n * (n - 1))
        
        return term1 + term2
    end
    
    h_opt = golden_section_search(lscv_objective_chunked, h_min, h_max)
    
    return h_opt
end

function kernel_eta(
    controls::Vector{Float64},
    cases::Vector{Float64};
    method::String = "optimal",
    t0::Float64 = 1.0,
    mesh_size_kernel::Int = 1000,
    box_cox::Bool = false
)
    # Apply Box-Cox transformation if specified
    if box_cox
        transformed = apply_box_cox(controls, cases)
        controls = transformed.transformed_x
        cases = transformed.transformed_y
    end
    
    n_ctrl = length(controls)
    n_case = length(cases)
    
    # Bandwidth selection
    bandwidth_controls, bandwidth_cases = if method == "optimal"
        # Silverman's rule of thumb
        bw_ctrl = 1.06 * std(controls) * n_ctrl^(-0.2)
        bw_case = 1.06 * std(cases) * n_case^(-0.2)
        (bw_ctrl, bw_case)
    elseif method == "iqr"
        # Robust rule based on IQR
        iqr_ctrl = iqr(controls)
        iqr_case = iqr(cases)
        bw_ctrl = 0.9 * min(std(controls), iqr_ctrl / 1.34) * n_ctrl^(-0.2)
        bw_case = 0.9 * min(std(cases), iqr_case / 1.34) * n_case^(-0.2)
        (bw_ctrl, bw_case)
    elseif method == "hscv"
        bw_ctrl = hscv_hybrid(controls)
        bw_case = hscv_hybrid(cases)
        (bw_ctrl, bw_case)
        else
        error("Unknown method: $method. Use 'optimal', 'iqr', or 'hscv'.")
    end
    
    # Sort samples (important for kernel estimation stability)
    sorted_controls = sort(controls)
    sorted_cases = sort(cases)
    
    # Create mesh for estimation
    min_val = min(minimum(controls), minimum(cases))
    max_val = max(maximum(controls), maximum(cases))
    mesh = collect(range(min_val, max_val, length=mesh_size_kernel))
    
    # Estimate distributions and densities
    estimated_dist_controls = kernel_distribution_estimation(sorted_controls, mesh, bandwidth_controls)
    estimated_dist_cases = kernel_distribution_estimation(sorted_cases, mesh, bandwidth_cases)
    estimated_density_controls = kernel_density_estimation(sorted_controls, mesh, bandwidth_controls)
    estimated_density_cases = kernel_density_estimation(sorted_cases, mesh, bandwidth_cases)
    
    # Create probability sequence
    p = collect(range(0.0001, 0.999, length=mesh_size_kernel))
    
    # Allocate ROC arrays
    roc = Vector{Float64}(undef, mesh_size_kernel)
    roc_prime = Vector{Float64}(undef, mesh_size_kernel)
    
    # Compute ROC and ROC' for each probability point
    @inbounds for i in 1:mesh_size_kernel
        p_opp = 1.0 - p[i]
        
        # Find threshold via inverse CDF of controls
        inv = inverse_kernel_estimation(p_opp, estimated_dist_controls, mesh)
        
        # Evaluate densities at threshold
        numerator = evaluate_kernel_estimation(inv, estimated_density_cases, mesh)
        denominator = evaluate_kernel_estimation(inv, estimated_density_controls, mesh)
        
        # Compute ROC and ROC'
        roc[i] = 1.0 - evaluate_kernel_estimation(inv, estimated_dist_cases, mesh)
        
        # Handle division by zero (numerical stability)
        roc_prime[i] = if denominator > 1e-12
            numerator / denominator
        else
            1.0  # Default value when denominator is too small
        end
    end
    
    # Calculate eta from ROC curves
    return eta_from_roc_curves(roc, roc_prime, t0, p)
end




#########################
# get_power translated
#########################

function get_power(file_name::String;
                   use_box_cox_in_parametric::Bool = false,
                   use_box_cox_in_kernel::Bool = false,
                   param_adjuster_function = sum,
                   case::String = "gaussian",
                   controls_params = Dict(:param1=>1.0, :param2=>1.0),
                   cases_params = Dict(:param1=>1.1, :param2=>1.1),
                   MC::Int = 1000,
                   BootstrapSize::Int = 500,
                   alpha::Float64 = 0.05)

    Random.seed!(1)

    AUCs = [0.6, 0.75, 0.9]
    ns = [20, 50, 100]
    t0s = [0.2, 0.4, 0.8, 1.0]

    # Sample generator
    sample_distribution = nothing
    if case == "gaussian"
        sample_distribution = (n, p1, p2) -> rand(Normal(p1, p2), n)
    elseif case == "lognormal"
        sample_distribution = (n, p1, p2) -> rand(LogNormal(p1, p2), n)
    elseif case == "gamma"
        sample_distribution = (n, p1, p2) -> rand(Gamma(p1, 1.0 / p2), n)
    else
        error("Parameter 'case' must be 'gaussian', 'lognormal' or 'gamma'")
    end

    results_summary = Dict{String,Any}("header" => "Simulation results for $file_name")

    mesh = collect(range(0.00001, 0.9999, length=10000))

    for (auc_idx, auc) in enumerate(AUCs)
        auc_list = Dict{String,Any}()

        # obtain missing param (either provided per-AUC or computed via adjuster)
        test_sample = sample_distribution(100_000, controls_params[:param1], controls_params[:param2])
        if haskey(cases_params, :param1) && isa(cases_params[:param1], AbstractVector)
            missing_param = cases_params[:param1][auc_idx]
        else
            error("no missing param")
        end

        for t0 in t0s
            t0_list = Dict{String,Any}()

            for n in ns

                parametric_base = Vector{Float64}(undef, MC)
                kernel_hscv_base = Vector{Float64}(undef, MC)
                kernel_opt_base = Vector{Float64}(undef, MC)
                kernel_iqr_base = Vector{Float64}(undef, MC)
                auc_base = Vector{Float64}(undef, MC)
                youden_base = Vector{Float64}(undef, MC)

                parametric_boot = Array{Float64}(undef, BootstrapSize, MC)
                kernel_hscv_boot = Array{Float64}(undef, BootstrapSize, MC)
                kernel_opt_boot = Array{Float64}(undef, BootstrapSize, MC)
                kernel_iqr_boot = Array{Float64}(undef, BootstrapSize, MC)
                auc_boot = Array{Float64}(undef, BootstrapSize, MC)
                youden_boot = Array{Float64}(undef, BootstrapSize, MC)

                parametric_p = Vector{Float64}(undef, MC)
                kernel_hscv_p = Vector{Float64}(undef, MC)
                kernel_opt_p = Vector{Float64}(undef, MC)
                kernel_iqr_p = Vector{Float64}(undef, MC)
                auc_p = Vector{Float64}(undef, MC)
                youden_p = Vector{Float64}(undef, MC)

                for mc_it in 1:MC
                    controls = sample_distribution(n, controls_params[:param1], controls_params[:param2])
                    cases = sample_distribution(n, missing_param, cases_params[:param2])

                    parametric_base[mc_it] = parametric_eta(controls, cases, use_box_cox_in_parametric, t0, mesh)
                    kernel_hscv_base[mc_it] = kernel_eta(controls, cases; method="hscv", t0=t0, box_cox=use_box_cox_in_kernel)
                    kernel_opt_base[mc_it] = kernel_eta(controls, cases; method="optimal", t0=t0, box_cox=use_box_cox_in_kernel)
                    kernel_iqr_base[mc_it] = kernel_eta(controls, cases; method="iqr", t0=t0, box_cox=use_box_cox_in_kernel)
                    auc_base[mc_it] = max(calculate_auc_normal(cases, controls), calculate_auc_normal(controls, cases))
                    youden_base[mc_it] = max(calculate_youden_normal(cases, controls), calculate_youden_normal(controls, cases))

                    # bootstrap
                    for bc_it in 1:BootstrapSize
                        combined_boot = sample(vcat(controls, cases), 2n; replace=true)
                        controls_b = combined_boot[1:n]
                        cases_b = combined_boot[n+1:2n]

                        parametric_boot[bc_it, mc_it] = parametric_eta(controls_b, cases_b, use_box_cox_in_parametric, t0, mesh)
                        kernel_hscv_boot[bc_it, mc_it] = kernel_eta(controls_b, cases_b; method="hscv", t0=t0, box_cox=use_box_cox_in_kernel)
                        kernel_opt_boot[bc_it, mc_it] = kernel_eta(controls_b, cases_b; method="optimal", t0=t0, box_cox=use_box_cox_in_kernel)
                        kernel_iqr_boot[bc_it, mc_it] = kernel_eta(controls_b, cases_b; method="iqr", t0=t0, box_cox=use_box_cox_in_kernel)
                        auc_boot[bc_it, mc_it] = max(calculate_auc_normal(cases_b, controls_b), calculate_auc_normal(controls_b, cases_b))
                        youden_boot[bc_it, mc_it] = max(calculate_youden_normal(cases_b, controls_b), calculate_youden_normal(controls_b, cases_b))
                    end

                    # p-values from bootstrap
                    parametric_p[mc_it] = mean(parametric_boot[:, mc_it] .>= parametric_base[mc_it])
                    kernel_hscv_p[mc_it] = mean(kernel_hscv_boot[:, mc_it] .>= kernel_hscv_base[mc_it])
                    kernel_opt_p[mc_it] = mean(kernel_opt_boot[:, mc_it] .>= kernel_opt_base[mc_it])
                    kernel_iqr_p[mc_it] = mean(kernel_iqr_boot[:, mc_it] .>= kernel_iqr_base[mc_it])
                    auc_p[mc_it] = mean(auc_boot[:, mc_it] .>= auc_base[mc_it])
                    youden_p[mc_it] = mean(youden_boot[:, mc_it] .>= youden_base[mc_it])
                end

                # summary metrics (allow mixed value types)
                summary_result = Dict{String,Any}(
                    "parametric" => Dict("power" => mean(parametric_p .< alpha)),
                    "kernel_hscv" => Dict("power" => mean(kernel_hscv_p .< alpha)),
                    "kernel_opt" => Dict("power" => mean(kernel_opt_p .< alpha)),
                    "kernel_iqr" => Dict("power" => mean(kernel_iqr_p .< alpha)),
                    "auc" => Dict("power" => mean(auc_p .< alpha)),
                    "youden" => Dict("power" => mean(youden_p .< alpha))
                )

                # save detailed RDS-equivalent as JSON (compressed not applied)
                dir_path = joinpath(pwd(), "results_powers_sim_julia")
                isdir(dir_path) || mkpath(dir_path)
                rds_path = joinpath(dir_path, string(file_name, "_n", n, "_t0", replace(string(t0), "."=>"_"), "_auc", replace(string(auc), "."=>"_"), ".json"))

                # Convert bootstrap matrices to arrays of vectors for JSON compatibility
                parametric_boot_json = [vec(parametric_boot[:, i]) for i in 1:size(parametric_boot, 2)]
                kernel_hscv_boot_json = [vec(kernel_hscv_boot[:, i]) for i in 1:size(kernel_hscv_boot, 2)]
                kernel_opt_boot_json = [vec(kernel_opt_boot[:, i]) for i in 1:size(kernel_opt_boot, 2)]
                kernel_iqr_boot_json = [vec(kernel_iqr_boot[:, i]) for i in 1:size(kernel_iqr_boot, 2)]
                auc_boot_json = [vec(auc_boot[:, i]) for i in 1:size(auc_boot, 2)]
                youden_boot_json = [vec(youden_boot[:, i]) for i in 1:size(youden_boot, 2)]

                detailed = Dict(
                    "meta" => Dict("file"=>file_name, "n"=>n, "t0"=>t0, "auc"=>auc, "case"=>case),
                    "base" => Dict("parametric"=>collect(parametric_base), "kernel_hscv"=>collect(kernel_hscv_base), "kernel_opt"=>collect(kernel_opt_base), "kernel_iqr"=>collect(kernel_iqr_base), "auc"=>collect(auc_base), "youden"=>collect(youden_base)),
                    "bootstrap" => Dict("parametric"=>parametric_boot_json, "kernel_hscv"=>kernel_hscv_boot_json, "kernel_opt"=>kernel_opt_boot_json, "kernel_iqr"=>kernel_iqr_boot_json, "auc"=>auc_boot_json, "youden"=>youden_boot_json),
                    "p_values" => Dict("parametric"=>collect(parametric_p), "kernel_hscv"=>collect(kernel_hscv_p), "kernel_opt"=>collect(kernel_opt_p), "kernel_iqr"=>collect(kernel_iqr_p), "auc"=>collect(auc_p), "youden"=>collect(youden_p))
                )

                # Try to serialize; on failure write a types map for debugging
                safe_detailed = stringify_keys(detailed)
                try
                    open(rds_path, "w") do io
                        JSON.print(io, safe_detailed)
                    end
                catch e
                    types_map = json_type_map(detailed)
                    types_path = rds_path * ".types.json"
                    open(types_path, "w") do io
                        JSON.print(io, types_map)
                    end
                    println("JSON serialization failed for ", rds_path, ". Types map written to ", types_path)
                    rethrow(e)
                end

                summary_result["rds"] = String(split(rds_path, "/") |> last)
                t0_list[string("size:", n)] = summary_result
            end
            auc_list[string("t0:", t0)] = t0_list
        end
        results_summary[string("AUC:", auc)] = auc_list
    end

    # write summary JSON
    json_path = joinpath(pwd(), "results_powers_sim_julia", string(file_name, "_summary.json"))
    open(json_path, "w") do io
        JSON.print(io, results_summary)
    end

    return (summary = json_path, detailed_dir = joinpath(pwd(), "results_powers_sim_julia"))
end

# Scenario definitions (only gaussian as in original)
gaussian_configs = Dict(
    :normal_1 => Dict(:file_name => "normal_1_box_cox_parametric_and_kernel", :case => "gaussian", :use_box_cox_in_parametric => true, :use_box_cox_in_kernel => true, :param_adjuster_function => sum, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[0.3515625, 0.94970703125, 1.81640625], :param2=>1.0)),
    :normal_2 => Dict(:file_name => "normal_2_box_cox_parametric_and_kernel", :case => "gaussian", :use_box_cox_in_parametric => true, :use_box_cox_in_kernel => true, :param_adjuster_function => sum, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[0.4296875, 1.15234375, 2.20703125],:param2=>1.4)),
    :normal_3 => Dict(:file_name => "normal_3_box_cox_parametric_and_kernel", :case => "gaussian", :use_box_cox_in_parametric => true, :use_box_cox_in_kernel => true, :param_adjuster_function => sum, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[0.80078125, 2.13623046875, 4.0625],:param2=>3.0))
)

lognormal_configs = Dict(
    :lognormal_1 => Dict(:file_name => "lognormal_1_box_cox_parametric", :case => "lognormal", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[1.4013671875, 2.0703125, 3.671875],:param2=>0.5)),
    :lognormal_2 => Dict(:file_name => "lognormal_2_box_cox_parametric", :case => "lognormal", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[1.806640625, 2.705078125, 4.3359375],:param2=>3/2)),
    :lognormal_3 => Dict(:file_name => "lognormal_3_box_cox_parametric", :case => "lognormal", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[1.30859375, 1.982421875, 3.6328125],:param2=>0.2)),
    :lognormal_4 => Dict(:file_name => "lognormal_4_box_cox_parametric", :case => "lognormal", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[1.982421875, 3.046875, 4.8828125],:param2=>2.0)),
    :lognormal_1_bc => Dict(:file_name => "lognormal_1_box_cox_parametric_and_kernel", :case => "lognormal", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :use_box_cox_in_kernel => true, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[1.4013671875, 2.0703125, 3.671875],:param2=>0.5)),
    :lognormal_2_bc => Dict(:file_name => "lognormal_2_box_cox_parametric_and_kernel", :case => "lognormal", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :use_box_cox_in_kernel => true, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[1.806640625, 2.705078125, 4.3359375],:param2=>3/2)),
    :lognormal_3_bc => Dict(:file_name => "lognormal_3_box_cox_parametric_and_kernel", :case => "lognormal", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :use_box_cox_in_kernel => true, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[1.30859375, 1.982421875, 3.6328125],:param2=>0.2)),
    :lognormal_4_bc => Dict(:file_name => "lognormal_4_box_cox_parametric_and_kernel", :case => "lognormal", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :use_box_cox_in_kernel => true, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[1.982421875, 3.046875, 4.8828125],:param2=>2.0)),
)

gamma_configs = Dict(
    :gamma_1 => Dict(:file_name => "gamma_1_box_cox_parametric", :case => "gamma", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :controls_params => Dict(:param1=>0.5,:param2=>0.5), :cases_params => Dict(:param1=>[1.083984375, 1.845703125, 3.57421875],:param2=>1.0)),
    :gamma_2 => Dict(:file_name => "gamma_2_box_cox_parametric", :case => "gamma", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :controls_params => Dict(:param1=>0.5,:param2=>0.5), :cases_params => Dict(:param1=>[3.251953125, 5.80078125, 11.8359375],:param2=>4.0)),
    :gamma_3 => Dict(:file_name => "gamma_3_box_cox_parametric", :case => "gamma", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :controls_params => Dict(:param1=>0.5,:param2=>0.5), :cases_params => Dict(:param1=>[0.36346435546875, 0.574951171875, 1.025390625],:param2=>1/8)),
    :gamma_1_bc => Dict(:file_name => "gamma_1_box_cox_parametric_and_kernel", :case => "gamma", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :use_box_cox_in_kernel => true, :controls_params => Dict(:param1=>0.5,:param2=>0.5), :cases_params => Dict(:param1=>[1.083984375, 1.845703125, 3.57421875],:param2=>1.0)),
    :gamma_2_bc => Dict(:file_name => "gamma_2_box_cox_parametric_and_kernel", :case => "gamma", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :use_box_cox_in_kernel => true, :controls_params => Dict(:param1=>0.5,:param2=>0.5), :cases_params => Dict(:param1=>[3.251953125, 5.80078125, 11.8359375],:param2=>4.0)),
    :gamma_3_bc => Dict(:file_name => "gamma_3_box_cox_parametric_and_kernel", :case => "gamma", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :use_box_cox_in_kernel => true, :controls_params => Dict(:param1=>0.5,:param2=>0.5), :cases_params => Dict(:param1=>[0.36346435546875, 0.574951171875, 1.025390625],:param2=>1/8)),
)

all_configs = Dict(:gaussian => gaussian_configs, :lognormal => lognormal_configs, :gamma => gamma_configs)



#########################
# get_power translated
#########################

function get_power2(file_name::String;
                   use_box_cox_in_parametric::Bool = false,
                   use_box_cox_in_kernel::Bool = false,
                   param_adjuster_function = sum,
                   case::String = "gaussian",
                   controls_params = Dict(:param1=>1.0, :param2=>1.0),
                   cases_params = Dict(:param1=>1.1, :param2=>1.1),
                   MC::Int = 1000,
                   BootstrapSize::Int = 500,
                   alpha::Float64 = 0.05)

    Random.seed!(1)

    AUCs = [0.6, 0.75, 0.9]
    ns = [20, 50, 100]
    t0s = [0.2, 0.4, 0.8, 1.0]

    # Sample generator
    sample_distribution = nothing
    if case == "gaussian"
        sample_distribution = (n, p1, p2) -> rand(Normal(p1, p2), n)
    elseif case == "lognormal"
        sample_distribution = (n, p1, p2) -> rand(LogNormal(p1, p2), n)
    elseif case == "gamma"
        sample_distribution = (n, p1, p2) -> rand(Gamma(p1, 1.0 / p2), n)
    else
        error("Parameter 'case' must be 'gaussian', 'lognormal' or 'gamma'")
    end

    results_summary = Dict{String,Any}("header" => "Simulation results for $file_name")

    mesh = collect(range(0.00001, 0.9999, length=10000))

    for (auc_idx, auc) in enumerate(AUCs)
        auc_list = Dict{String,Any}()

        # obtain missing param (either provided per-AUC or computed via adjuster)
        test_sample = sample_distribution(100_000, controls_params[:param1], controls_params[:param2])
        if haskey(cases_params, :param1) && isa(cases_params[:param1], AbstractVector)
            missing_param = cases_params[:param1][auc_idx]
        else
            error("no missing param")
        end

        for t0 in t0s
            t0_list = Dict{String,Any}()

            for n in ns

                parametric_base = Vector{Float64}(undef, MC)
                kernel_hscv_base = Vector{Float64}(undef, MC)
                kernel_opt_base = Vector{Float64}(undef, MC)
                kernel_iqr_base = Vector{Float64}(undef, MC)
                auc_base = Vector{Float64}(undef, MC)
                youden_base = Vector{Float64}(undef, MC)

                parametric_boot = Array{Float64}(undef, BootstrapSize, MC)
                kernel_hscv_boot = Array{Float64}(undef, BootstrapSize, MC)
                kernel_opt_boot = Array{Float64}(undef, BootstrapSize, MC)
                kernel_iqr_boot = Array{Float64}(undef, BootstrapSize, MC)
                auc_boot = Array{Float64}(undef, BootstrapSize, MC)
                youden_boot = Array{Float64}(undef, BootstrapSize, MC)

                parametric_p = Vector{Float64}(undef, MC)
                kernel_hscv_p = Vector{Float64}(undef, MC)
                kernel_opt_p = Vector{Float64}(undef, MC)
                kernel_iqr_p = Vector{Float64}(undef, MC)
                auc_p = Vector{Float64}(undef, MC)
                youden_p = Vector{Float64}(undef, MC)

                for mc_it in 1:MC
                    controls = sample_distribution(n, controls_params[:param1], controls_params[:param2])
                    cases = sample_distribution(n, missing_param, cases_params[:param2])

                    parametric_base[mc_it] = parametric_eta(controls, cases, use_box_cox_in_parametric, t0, mesh)
                    kernel_hscv_base[mc_it] = kernel_eta(controls, cases; method="hscv", t0=t0, box_cox=use_box_cox_in_kernel)
                    kernel_opt_base[mc_it] = kernel_eta(controls, cases; method="optimal", t0=t0, box_cox=use_box_cox_in_kernel)
                    kernel_iqr_base[mc_it] = kernel_eta(controls, cases; method="iqr", t0=t0, box_cox=use_box_cox_in_kernel)
                    auc_base[mc_it] = max(calculate_auc_normal(cases, controls), calculate_auc_normal(controls, cases))
                    youden_base[mc_it] = max(calculate_youden_normal(cases, controls), calculate_youden_normal(controls, cases))

                    # bootstrap
                    for bc_it in 1:BootstrapSize
                        combined_boot = sample(vcat(controls, cases), 2n; replace=true)
                        controls_b = combined_boot[1:n]
                        cases_b = combined_boot[n+1:2n]

                        parametric_boot[bc_it, mc_it] = parametric_eta(controls_b, cases_b, use_box_cox_in_parametric, t0, mesh)
                        kernel_hscv_boot[bc_it, mc_it] = kernel_eta(controls_b, cases_b; method="hscv", t0=t0, box_cox=use_box_cox_in_kernel)
                        kernel_opt_boot[bc_it, mc_it] = kernel_eta(controls_b, cases_b; method="optimal", t0=t0, box_cox=use_box_cox_in_kernel)
                        kernel_iqr_boot[bc_it, mc_it] = kernel_eta(controls_b, cases_b; method="iqr", t0=t0, box_cox=use_box_cox_in_kernel)
                        auc_boot[bc_it, mc_it] = max(calculate_auc_normal(cases_b, controls_b), calculate_auc_normal(controls_b, cases_b))
                        youden_boot[bc_it, mc_it] = max(calculate_youden_normal(cases_b, controls_b), calculate_youden_normal(controls_b, cases_b))
                    end

                    # p-values from bootstrap
                    parametric_p[mc_it] = mean(parametric_boot[:, mc_it] .>= parametric_base[mc_it])
                    kernel_hscv_p[mc_it] = mean(kernel_hscv_boot[:, mc_it] .>= kernel_hscv_base[mc_it])
                    kernel_opt_p[mc_it] = mean(kernel_opt_boot[:, mc_it] .>= kernel_opt_base[mc_it])
                    kernel_iqr_p[mc_it] = mean(kernel_iqr_boot[:, mc_it] .>= kernel_iqr_base[mc_it])
                    auc_p[mc_it] = mean(auc_boot[:, mc_it] .>= auc_base[mc_it])
                    youden_p[mc_it] = mean(youden_boot[:, mc_it] .>= youden_base[mc_it])
                end

                # summary metrics (allow mixed value types)
                summary_result = Dict{String,Any}(
                    "parametric" => Dict("power" => mean(parametric_p .< alpha)),
                    "kernel_hscv" => Dict("power" => mean(kernel_hscv_p .< alpha)),
                    "kernel_opt" => Dict("power" => mean(kernel_opt_p .< alpha)),
                    "kernel_iqr" => Dict("power" => mean(kernel_iqr_p .< alpha)),
                    "auc" => Dict("power" => mean(auc_p .< alpha)),
                    "youden" => Dict("power" => mean(youden_p .< alpha))
                )

                # save detailed RDS-equivalent as JSON (compressed not applied)
                dir_path = joinpath(pwd(), "results_powers_sim_julia")
                isdir(dir_path) || mkpath(dir_path)
                rds_path = joinpath(dir_path, string(file_name, "_n", n, "_t0", replace(string(t0), "."=>"_"), "_auc", replace(string(auc), "."=>"_"), ".json"))

                # Convert bootstrap matrices to arrays of vectors for JSON compatibility
                parametric_boot_json = [vec(parametric_boot[:, i]) for i in 1:size(parametric_boot, 2)]
                kernel_hscv_boot_json = [vec(kernel_hscv_boot[:, i]) for i in 1:size(kernel_hscv_boot, 2)]
                kernel_opt_boot_json = [vec(kernel_opt_boot[:, i]) for i in 1:size(kernel_opt_boot, 2)]
                kernel_iqr_boot_json = [vec(kernel_iqr_boot[:, i]) for i in 1:size(kernel_iqr_boot, 2)]
                auc_boot_json = [vec(auc_boot[:, i]) for i in 1:size(auc_boot, 2)]
                youden_boot_json = [vec(youden_boot[:, i]) for i in 1:size(youden_boot, 2)]

                detailed = Dict(
                    "meta" => Dict("file"=>file_name, "n"=>n, "t0"=>t0, "auc"=>auc, "case"=>case, "true_eta"=>true_eta),
                    "base" => Dict("parametric"=>collect(parametric_base), "kernel_hscv"=>collect(kernel_hscv_base), "kernel_opt"=>collect(kernel_opt_base), "kernel_iqr"=>collect(kernel_iqr_base), "auc"=>collect(auc_base), "youden"=>collect(youden_base)),
                    "bootstrap" => Dict("parametric"=>parametric_boot_json, "kernel_hscv"=>kernel_hscv_boot_json, "kernel_opt"=>kernel_opt_boot_json, "kernel_iqr"=>kernel_iqr_boot_json, "auc"=>auc_boot_json, "youden"=>youden_boot_json),
                    "p_values" => Dict("parametric"=>collect(parametric_p), "kernel_hscv"=>collect(kernel_hscv_p), "kernel_opt"=>collect(kernel_opt_p), "kernel_iqr"=>collect(kernel_iqr_p), "auc"=>collect(auc_p), "youden"=>collect(youden_p))
                )

                # Try to serialize; on failure write a types map for debugging
                safe_detailed = stringify_keys(detailed)
                try
                    open(rds_path, "w") do io
                        JSON.print(io, safe_detailed)
                    end
                catch e
                    types_map = json_type_map(detailed)
                    types_path = rds_path * ".types.json"
                    open(types_path, "w") do io
                        JSON.print(io, types_map)
                    end
                    println("JSON serialization failed for ", rds_path, ". Types map written to ", types_path)
                    rethrow(e)
                end

                summary_result["rds"] = String(split(rds_path, "/") |> last)
                t0_list[string("size:", n)] = summary_result
            end
            auc_list[string("t0:", t0)] = t0_list
        end
        results_summary[string("AUC:", auc)] = auc_list
    end

    # write summary JSON
    json_path = joinpath(pwd(), "results_powers_sim_julia", string(file_name, "_summary.json"))
    open(json_path, "w") do io
        JSON.print(io, results_summary)
    end

    return (summary = json_path, detailed_dir = joinpath(pwd(), "results_powers_sim_julia"))
end

# Scenario definitions (only gaussian as in original)
gaussian_configs = Dict(
    :normal_1 => Dict(:file_name => "normal_1_box_cox_parametric_and_kernel", :case => "gaussian", :use_box_cox_in_parametric => true, :use_box_cox_in_kernel => true, :param_adjuster_function => sum, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[0.3515625, 0.94970703125, 1.81640625], :param2=>1.0)),
    :normal_2 => Dict(:file_name => "normal_2_box_cox_parametric_and_kernel", :case => "gaussian", :use_box_cox_in_parametric => true, :use_box_cox_in_kernel => true, :param_adjuster_function => sum, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[0.4296875, 1.15234375, 2.20703125],:param2=>1.4)),
    :normal_3 => Dict(:file_name => "normal_3_box_cox_parametric_and_kernel", :case => "gaussian", :use_box_cox_in_parametric => true, :use_box_cox_in_kernel => true, :param_adjuster_function => sum, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[0.80078125, 2.13623046875, 4.0625],:param2=>3.0))
)

lognormal_configs = Dict(
    :lognormal_1 => Dict(:file_name => "lognormal_1_box_cox_parametric", :case => "lognormal", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[1.4013671875, 2.0703125, 3.671875],:param2=>0.5)),
    :lognormal_2 => Dict(:file_name => "lognormal_2_box_cox_parametric", :case => "lognormal", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[1.806640625, 2.705078125, 4.3359375],:param2=>3/2)),
    :lognormal_3 => Dict(:file_name => "lognormal_3_box_cox_parametric", :case => "lognormal", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[1.30859375, 1.982421875, 3.6328125],:param2=>0.2)),
    :lognormal_4 => Dict(:file_name => "lognormal_4_box_cox_parametric", :case => "lognormal", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[1.982421875, 3.046875, 4.8828125],:param2=>2.0)),
    :lognormal_1_bc => Dict(:file_name => "lognormal_1_box_cox_parametric_and_kernel", :case => "lognormal", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :use_box_cox_in_kernel => true, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[1.4013671875, 2.0703125, 3.671875],:param2=>0.5)),
    :lognormal_2_bc => Dict(:file_name => "lognormal_2_box_cox_parametric_and_kernel", :case => "lognormal", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :use_box_cox_in_kernel => true, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[1.806640625, 2.705078125, 4.3359375],:param2=>3/2)),
    :lognormal_3_bc => Dict(:file_name => "lognormal_3_box_cox_parametric_and_kernel", :case => "lognormal", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :use_box_cox_in_kernel => true, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[1.30859375, 1.982421875, 3.6328125],:param2=>0.2)),
    :lognormal_4_bc => Dict(:file_name => "lognormal_4_box_cox_parametric_and_kernel", :case => "lognormal", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :use_box_cox_in_kernel => true, :controls_params => Dict(:param1=>0.0,:param2=>1.0), :cases_params => Dict(:param1=>[1.982421875, 3.046875, 4.8828125],:param2=>2.0)),
)

gamma_configs = Dict(
    :gamma_1 => Dict(:file_name => "gamma_1_box_cox_parametric", :case => "gamma", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :controls_params => Dict(:param1=>0.5,:param2=>0.5), :cases_params => Dict(:param1=>[1.083984375, 1.845703125, 3.57421875],:param2=>1.0)),
    :gamma_2 => Dict(:file_name => "gamma_2_box_cox_parametric", :case => "gamma", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :controls_params => Dict(:param1=>0.5,:param2=>0.5), :cases_params => Dict(:param1=>[3.251953125, 5.80078125, 11.8359375],:param2=>4.0)),
    :gamma_3 => Dict(:file_name => "gamma_3_box_cox_parametric", :case => "gamma", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :controls_params => Dict(:param1=>0.5,:param2=>0.5), :cases_params => Dict(:param1=>[0.36346435546875, 0.574951171875, 1.025390625],:param2=>1/8)),
    :gamma_1_bc => Dict(:file_name => "gamma_1_box_cox_parametric_and_kernel", :case => "gamma", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :use_box_cox_in_kernel => true, :controls_params => Dict(:param1=>0.5,:param2=>0.5), :cases_params => Dict(:param1=>[1.083984375, 1.845703125, 3.57421875],:param2=>1.0)),
    :gamma_2_bc => Dict(:file_name => "gamma_2_box_cox_parametric_and_kernel", :case => "gamma", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :use_box_cox_in_kernel => true, :controls_params => Dict(:param1=>0.5,:param2=>0.5), :cases_params => Dict(:param1=>[3.251953125, 5.80078125, 11.8359375],:param2=>4.0)),
    :gamma_3_bc => Dict(:file_name => "gamma_3_box_cox_parametric_and_kernel", :case => "gamma", :param_adjuster_function => sum, :use_box_cox_in_parametric => true, :use_box_cox_in_kernel => true, :controls_params => Dict(:param1=>0.5,:param2=>0.5), :cases_params => Dict(:param1=>[0.36346435546875, 0.574951171875, 1.025390625],:param2=>1/8)),
)

all_configs = Dict(:gaussian => gaussian_configs, :lognormal => lognormal_configs, :gamma => gamma_configs)

gaussians = Dict(:gaussian => gaussian_configs)

function simulate(category_configs::Dict; MC::Int=100, BootstrapSize::Int=50, alpha::Float64=0.05, verbose::Bool=true)
    results = Dict{Any,Any}()
    for (category_name, configs) in category_configs
        if verbose
            println("Running scenarios for: ", category_name)
        end
        cat_results = Dict{Any,Any}()
        for (cfg_name, cfg) in configs
            try
                if verbose
                    println("  Running: ", cfg[:file_name])
                end
                res = get_power(cfg[:file_name];
                                use_box_cox_in_parametric = get(cfg, :use_box_cox_in_parametric, false),
                                use_box_cox_in_kernel = get(cfg, :use_box_cox_in_kernel, false),
                                param_adjuster_function = get(cfg, :param_adjuster_function, sum),
                                case = get(cfg, :case, "gaussian"),
                                controls_params = get(cfg, :controls_params, Dict(:param1=>1.0,:param2=>1.0)),
                                cases_params = get(cfg, :cases_params, Dict(:param1=>1.1,:param2=>1.1)),
                                MC = MC,
                                BootstrapSize = BootstrapSize,
                                alpha = alpha)
                cat_results[cfg_name] = res
                if verbose
                    println("    Finished: ", cfg[:file_name], " -> ", res.summary)
                end
            catch e
                cat_results[cfg_name] = Dict("error" => sprint(showerror, e))
                if verbose
                    println("    Error in: ", cfg[:file_name], " -> ", e)
                end
            end
        end
        results[category_name] = cat_results
    end
    return results
end


function simulate2(category_configs::Dict; MC::Int=100, BootstrapSize::Int=50, alpha::Float64=0.05, verbose::Bool=true)
    results = Dict{Any,Any}()
    for (category_name, configs) in category_configs
        if verbose
            println("Running scenarios for: ", category_name)
        end
        cat_results = Dict{Any,Any}()
        for (cfg_name, cfg) in configs
            try
                if verbose
                    println("  Running: ", cfg[:file_name])
                end
                res = get_power2(cfg[:file_name];
                                use_box_cox_in_parametric = get(cfg, :use_box_cox_in_parametric, false),
                                use_box_cox_in_kernel = get(cfg, :use_box_cox_in_kernel, false),
                                param_adjuster_function = get(cfg, :param_adjuster_function, sum),
                                case = get(cfg, :case, "gaussian"),
                                controls_params = get(cfg, :controls_params, Dict(:param1=>1.0,:param2=>1.0)),
                                cases_params = get(cfg, :cases_params, Dict(:param1=>1.1,:param2=>1.1)),
                                MC = MC,
                                BootstrapSize = BootstrapSize,
                                alpha = alpha)
                cat_results[cfg_name] = res
                if verbose
                    println("    Finished: ", cfg[:file_name], " -> ", res.summary)
                end
            catch e
                cat_results[cfg_name] = Dict("error" => sprint(showerror, e))
                if verbose
                    println("    Error in: ", cfg[:file_name], " -> ", e)
                end
            end
        end
        results[category_name] = cat_results
    end
    return results
end


