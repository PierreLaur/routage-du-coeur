using DataFrames
using GLM
using XLSX
using StatsBase
using Statistics

distance_matrix = XLSX.readtable("data/euclidean_matrix.xlsx", 1) |> DataFrame
duration_matrix = XLSX.readtable("data/duration_matrix.xlsx", 1) |> DataFrame

distances = reshape(Matrix(distance_matrix), :)
durations = reshape(Matrix(duration_matrix), :)

df = DataFrame(Distance=distances, Duration=durations) .|> Float64
model = lm(@formula(Duration ~ 1 + Distance), df)

error = mean(sqrt.((predict(model, df) .- durations) .^ 2))

a, b = coef(model_1)

a
b