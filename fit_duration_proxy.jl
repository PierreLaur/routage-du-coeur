using DataFrames
using GLM
using XLSX
using StatsBase

distance_matrix = XLSX.readtable("data/euclidean_matrix.xlsx", 1) |> DataFrame
duration_matrix = XLSX.readtable("data/duration_matrix_w_traffic.xlsx", 1) |> DataFrame

distances = reshape(Matrix(distance_matrix), :) ./ 1000 # distances in km
durations = reshape(Matrix(duration_matrix), :) # durations in minutes

df = DataFrame(Distance=distances, Duration=durations) .|> Float64

using Plots
scatter(df.Distance, df.Duration)

model_1 = lm(@formula(Duration ~ 1 + Distance), df)
r²(model_1)
model_2 = lm(@formula(Duration ~ 1 + Distance + Distance^2), df)
r²(model_2)
model_3 = lm(@formula(Duration ~ 1 + Distance + Distance^2 + Distance^3), df)
r²(model_3)
model_4 = lm(@formula(Duration ~ 1 + Distance + Distance^2 + Distance^3 + Distance^4), df)
r²(model_4)

adjr²(model_1)
adjr²(model_2)
adjr²(model_3)
adjr²(model_4)

error = mean(sqrt.((predict(model_1, df) .- durations) .^ 2))
error = mean(sqrt.((predict(model_2, df) .- durations) .^ 2))
error = mean(sqrt.((predict(model_3, df) .- durations) .^ 2))
error = mean(sqrt.((predict(model_4, df) .- durations) .^ 2))

indices = sort(eachindex(df.Distance), by=x -> df.Distance[x])

plot!(df.Distance[indices], predict(model_1, df)[indices])
plot!(df.Distance[indices], predict(model_2, df)[indices])
plot!(df.Distance[indices], predict(model_3, df)[indices])
plot!(df.Distance[indices], predict(model_4, df)[indices])

a, b = coef(model_1)

a
b / 1000