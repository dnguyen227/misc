using DisjunctiveProgramming, InfiniteOpt, Gurobi, Plots

# Set parameters
Tsp = 1.38
num_grid_pts = 111
# Add the steady-state diffusion constraints
D = 0.05

xs = collect(LinRange(-1, 1, num_grid_pts))
num_heaters = 6 # choose such that (num_grid_pts - num_heaters) / (num_heaters + 1) is an integer

# Determine heater placements
spacing = Int((num_grid_pts - num_heaters) / (num_heaters + 1)) + 1
heater_positions = xs[spacing:spacing:num_grid_pts] 

# Initialize the model
model = InfiniteGDPModel(Gurobi.Optimizer)
@infinite_parameter(model, x ∈ [-1, 1], independent = true, num_supports = num_grid_pts, derivative_method = FiniteDifference(Central()))

# Create the variables
@variable(model, 0 ≤ T, Infinite(x))
@variable(model, u == 0, Infinite(x)) # heater input

# Make it so heater only heats at certain points
for xp in heater_positions
        u_pt = u(xp)
        unfix(u_pt)
        set_lower_bound(u_pt, 0.0)
        set_upper_bound(u_pt, 50)
        set_start_value(u_pt, 25)
end

# Set up global constraints
@constraint(model, D * ∂(T, x, x) + u - 0.3 == 0)
@constraint(model, T(-1) == 0)
@constraint(model, T(1) == 0)

# Set up disjunctions
@variable(model, W[1:2], InfiniteLogical(x))
@variable(model, Q, InfiniteLogical(x))
@constraint(model, T ≤ 1.9, Disjunct(W[1]))
@constraint(model, T ≤ 1.4, Disjunct(W[2]))
@disjunction(model, W)

# Add logical proposition
@constraint(model, W[1] ⟹ Q := true)

# Set the objective function
@objective(model, Min, ∫((T - Tsp)^2 + 1E-02*binary_variable(Q), x))

# Choose a reformulation method and solve
set_silent(model)
method = "bigM"
if method == "indicator"
    optimize!(model, gdp_method = Indicator())
elseif method == "bigM"
    optimize!(model, gdp_method = BigM())
elseif method == "hull"
    optimize!(model, gdp_method = Hull())
else
    error("Please specify a GDP method.")
end

# Check & extract the result values for plotting
if has_values(model)
    x_vals = value.(x)
    Tsp_opt = ones(length(x_vals))*Tsp
    u_opt = value.(u)
    T_opt = value.(T)
    Y1_opt = value.(Y[1])
    Y2_opt = value.(Y[2])
    obj_opt = objective_value(model)
end

# Plot the results
mode_combined = [findfirst(x -> x == 1, [Y1_opt[i], Y2_opt[i]]) for i in 1:length(Y1_opt)]
Tmax_mode = [Tmax == 1 ? 1.9 : 1.4 for Tmax in mode_combined]
l = @layout [a ; b]
p1 = bar(x_vals, u_opt, label = "", xlabel = "Spatial position (x)", ylabel = "Heater input (u)", color="green")
p2 = plot(x_vals, T_opt, label = "", xlabel = "Spatial position (x)", ylabel = "Temperature (T)")
plot!(x_vals, Tsp_opt, label = "Tsp", linestyle=:dot)
plot!(x_vals, Tmax_mode, label = "Tmax", xlabel = "Spatial position (x)", colour="violet")
display(plot(p2, p1, layout = l, size = (700, 600)))
savefig("optimal-strip.png")