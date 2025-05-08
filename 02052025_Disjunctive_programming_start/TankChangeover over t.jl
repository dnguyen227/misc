using DisjunctiveProgramming, InfiniteOpt, SCIP

m = InfiniteGDPModel(SCIP.Optimizer)
n = 6
I = 1:3
J = 1:n
D = [1,2,3,4,5,6]

@infinite_parameter(m, t ∈ [0,6] , num_supports = 360)
@variable(m, -4<= z <= 4)
@variable(m, W[1:3,J], InfiniteLogical(t))


for j in J
    @constraint(m, [W[1,j],W[2,j],W[3,j]] in Exactly(1))
end

@variable(m, -5 <= y[J] <= 5, Infinite(t))

@constraint(m, y[1](0) == 1)
@constraint(m, [j ∈ 1:5], y[j+1](D[j]) == y[j](D[j]))
# @constraint(m, [j ∈ 1:5], y(j) == y(j+1))
@constraint(m, [j ∈ J], deriv(y[j],t) == -2t + 0.3z -20*y[j], Disjunct(W[1,j]))
@constraint(m, [j ∈ J],deriv(y[j],t) == -2z + 0.4t -4, Disjunct(W[2,j]))
@constraint(m, [j ∈ J],deriv(y[j],t) == 2z + 4(t - y[j] - 1), Disjunct(W[3,j]))
# @disjunction(m, [j ∈ J], W[:, j])

@constraint(m, W[1,1] ∧ W[3,2] ∧ W[2,3] ∧ W[2,4] ∧ W[1,5] ∧ W[2,6] := true)

# @disjuntion(m, W)

# @objective(m, Min, integral(y^2,t))
@objective(m, Min, sum(integral(y[j]^2, t) for j in J))

# @objective(m, Min, sum(integral(y(j)^2, t, bounds = (j,j+1)) for j in 1:5))

optimize!(m, gdp_method = BigM())

# print(value.(W))


if termination_status(m) == MOI.OPTIMAL
    println("The solution is optimal.")
elseif termination_status(m) == MOI.INFEASIBLE
    println("The problem is infeasible.")
elseif termination_status(m) == MOI.UNBOUNDED
    println("The problem is unbounded.")
else
    println("The solution status is: ", termination_status(m))
end
println("Objective value: ", objective_value(m))

using Plots
# Extract the solution for y over time
t_values = collect(range(0, 6, length=360))  # Compute t as τ + j - 1 for each index j
y_values = vcat(value.(m[:y]))  # Extract y values for each j
y_values = reduce(vcat, y_values)  # Flatten the array
# println(y_values)
# y_values = value(y)  # Extract y values for each j

plot(t_values, y_values, label="y(t)", xlabel="Time (t)", ylabel="y", title="y vs t", legend=:topright)

# Extract the solution for i (disjunction index) over time
# i_values = value()
# # Plot i over time
# plot!(t_values, i_values, label="i(t)", xlabel="Time (t)", ylabel="i", title="y and i over time", legend=:topright)