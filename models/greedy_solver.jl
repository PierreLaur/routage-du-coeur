using CSV, DataFrames

function process_data()

	ignore_centres = ["Fronton"]

	centres = DataFrame(CSV.File("data/centres.csv", header = 1))
	points_de_ramasse = DataFrame(CSV.File("data/points_de_ramasse.csv", header = 1))
	vehicles = DataFrame(CSV.File("data/vehicules.csv", header = 1))
	matrix = DataFrame(CSV.File("data/euclidean_matrix.csv", header = 1))


	n = size(centres, 1)
	n_pdr = size(points_de_ramasse, 1)
	n_days = 5  # Single week scheduling
	m = size(vehicles, 1)

	capacities = vehicles[:, "Capacité (kg)"]
	sizes = vehicles[:, "Taille(Palettes)"]
	max_palette_capacity = 800

	demands = Dict(
		"a" => centres[:, "Tonnage Ambiant (kg)"],
		"f" => centres[:, "Tonnage Frais (kg)"],
		"s" => centres[:, "Tonnage Surgelé (kg)"],
	)

	# Ignore some centres (deliver them the other week)
	indexes = [findfirst(c -> c == cname, centres[:, "Nom"]) for cname in ignore_centres]
	for i in indexes
		demands["a"][i] = 0
		demands["f"][i] = 0
		demands["s"][i] = 0
	end

	# Add 15% for robustness
	for c in 2:n
		demands["a"][c] = trunc(Int, demands["a"][c] * 1.15)
		demands["f"][c] = trunc(Int, demands["f"][c] * 1.15)
		demands["s"][c] = trunc(Int, demands["s"][c] * 1.15)
	end

	for d in values(demands)
		d[1] = 0
	end

	freqs_pdr = points_de_ramasse[:, "Fréquence de Ramasse(/w)"]
	frais = [f == "Oui" for f in vehicles[:, "Frais"]]

	return matrix, n, n_pdr, n_days, m, capacities, sizes, max_palette_capacity, demands, freqs_pdr, frais

end

function solve()
	matrix, n, n_pdr, n_days, m, capacities, sizes, max_palette_capacity, demands, freqs_pdr, frais = process_data()

	# TODO: Implement solver
end
