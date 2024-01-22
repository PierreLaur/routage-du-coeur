use csv;
use json;
use io;
use localsolver;
use random;
use datetime;

function input() {

    usage = "Usage : localsolver ls_solver.lsp <init_solution or 'nil'> <desired output file path> <week number (1 or 2)>";
    if (outfile == nil) throw usage;
    if (week == nil) throw usage;

    params = json.parse("data/params.json") ;
    
    centres = csv.parse("data/centres.csv");
    pdr = csv.parse("data/points_de_ramasse.csv");
    matrix = csv.parse("data/euclidean_matrix.csv");
    vehicles = csv.parse("data/vehicules.csv");

    n_centres = centres.nbRows ;

    // Drop some unused centres (delivered the other week)
    // Create a map to be able to fetch info for centres in the datafiles
    index_to_centre = {} ;
    i = 0 ;
    for [c in 0...n_centres - 1] {
        semaine_livraison = centres.rows[c+1][7] ;
        if(semaine_livraison == 0 || semaine_livraison == week) {
            index_to_centre[i] = c ;
            i += 1 ;
        }
    }

    // The number of centres to serve this week
    n = index_to_centre.count();
    m = vehicles.nbRows;
    n_pdr = pdr.nbRows;
    n_nodes = n + n_pdr;
    n_days = 5 ;

    for [node in 0...n_nodes] {
        node_index = node < n ? index_to_centre[node] : node - n + n_centres - 1 ;

        dist_from_depot[node] = matrix.rows[0][node_index+2] ;
        dist_to_depot[node] = matrix.rows[node_index+1][1] ;

        for [node2 in 0...n_nodes] {
            node_index2 = node2 < n ? index_to_centre[node2] : node2 - n + n_centres - 1 ;
            
            dist_matrix[node][node2] = matrix.rows[node_index+1][node_index2+2];
        }
    }

    // Demands
    for [i in 0...n] {
        centre = index_to_centre[i] ;
        demands["a"][i] = ceil(centres.rows[centre+1][4] * 1.15); 
        demands["f"][i] = ceil(centres.rows[centre+1][5] * 1.15); // Add 15% for robustness
        demands["s"][i] = ceil(centres.rows[centre+1][6] * 1.15);

        /* 
        Note : if the demand is normal, we can enforce demand constraints satisfied 95% of the time by
            setting the demand to mu + Ksigma
            with K = Q0.975 = 1.645
            if sigma = X% of mu, this is mu(1 + KX)
        */

    }

    // Pickup weights (per pickup)
    for [p in 0...n_pdr] {
        weight[p] = ceil(pdr.rows[p][5] * 1.15);
    }

    for [i in 0...m] {
        row = vehicles.rows[i];
        // Capacities in kg
        capacities[i] = row[1];
        
        // Sizes in number of pallets
        sizes[i] = row[2];

        // Fuel consumption (L/100km)
        fuel[i] = row[3];

        frais[i] = row[4];
    }

    jours_map = {"Lundi": 0, "Mardi": 1, "Mercredi": 2, "Jeudi": 3, "Vendredi": 4};
    // When are we allowed to deliver each centre ?
    for [c in 0...n] {
        index = index_to_centre[c]+1 ;
        jours = centres.rows[index][10].split(", ") ;
        j_de_livraison[c] = {};
        for [k in 0...jours.count()] {
            j_de_livraison[c].add(jours_map[jours[k]]);
        }
    }

    // When do we have to pickup each site ?
    for [i in 0...n_pdr] {
        jours = pdr.rows[i][4].split(", ") ;
        j_de_ramasse[i] = {};
        for [k in 0...jours.count()] {
            j_de_ramasse[i].add(jours_map[jours[k]]);
        }
    }
}


function add_specific_requirements() {
    /* 
    Adds a few "hard-coded" constraints 
    These are not written in the input files but have been specified by the client 
    */

    // The PL must visit Carrefour Centrale on Wednesday
    constraint visits_r[0][2][0][5] ;
    // He must only do one pickup
    constraint n_ramasses[0][2][0] == 1;

    // The first camion Frigo must visit Carrefour Centrale on Wednesday and Friday
    constraint visits_r[0][2][2][5] ;
    constraint visits_r[0][4][2][5] ;
    // To simplify we say it only does one pickup (in reality he delivers a centre afterwards)
    constraint n_ramasses[0][2][2] == 1;
    constraint n_ramasses[0][4][2] == 1;

}

function model() {

    livraisons[h in 0...2][d in 0...n_days][v in 0...m] <- list(n);
    ramasses[h in 0...2][d in 0...n_days][v in 0...m] <- list(n_pdr);

    constraint cover[h in 0...2][d in 0...n_days][v in 0...m](livraisons[h][d][v]);
    constraint cover[h in 0...2][d in 0...n_days][v in 0...m](ramasses[h][d][v]);

    for [h in 0...2] {
        for [d in 0...n_days] {
            for [v in 0...m] {

                // Deliver quantity for each product type (ambiant/frais/surgelé)
                if (params["vehicle_allowed"][v]) {
                    serve_a[h][d][v][c in 0...n] <- int(0, demands["a"][c]) ;
                    serve_f[h][d][v][c in 0...n] <- frais[v] == "Oui" ? int(0, demands["f"][c]) : 0 ;
                    serve_s[h][d][v][c in 0...n] <- int(0, demands["s"][c]) ;
                    norvegiennes[h][d][v][c in 0...n] <- int(0, params["n_norvegiennes"]);
                } else {
                    serve_a[h][d][v][c in 0...n] <- 0 ;
                    serve_f[h][d][v][c in 0...n] <- 0 ;
                    serve_s[h][d][v][c in 0...n] <- 0 ;
                    norvegiennes[h][d][v][c in 0...n] <- 0 ;
                }

                n_livraisons[h][d][v] <- count(livraisons[h][d][v]) ;
                n_ramasses[h][d][v] <- count(ramasses[h][d][v]) ;

                local liv <- n_livraisons[h][d][v] > 0 ;
                local ram <- n_ramasses[h][d][v] > 0 ;

                visits_l[h][d][v][c in 0...n] <- contains(livraisons[h][d][v], c) ;
                visits_r[h][d][v][p in 0...n_pdr] <- contains(ramasses[h][d][v], p);
                
                for [c in 0...n] {
                    load[h][d][v][c] <- visits_l[h][d][v][c] * (serve_a[h][d][v][c] + serve_f[h][d][v][c] + serve_s[h][d][v][c]) ;
                    palettes[h][d][v][c] <- ceil(visits_l[h][d][v][c] * serve_a[h][d][v][c] / params["max_palette_capacity"]) ;
                    demi_palettes[h][d][v][c] <- ceil(visits_l[h][d][v][c] * serve_f[h][d][v][c] / params["demi_palette_capacity"]);

                    demi_palettes_s[h][d][v][c] <- frais[v] == "Oui" ? 
                                        ceil(
                                            (visits_l[h][d][v][c] * serve_s[h][d][v][c] - 
                                            norvegiennes[h][d][v][c] * params["norvegienne_capacity"]) / params["demi_palette_capacity"]
                                            ) : 
                                        0;

                    constraint visits_l[h][d][v][c] * serve_s[h][d][v][c] <= 
                                    norvegiennes[h][d][v][c] * params["norvegienne_capacity"] +
                                    demi_palettes_s[h][d][v][c] * params["demi_palette_capacity"] ;
                    constraint demi_palettes_s[h][d][v][c] >= 0 ;
                } 


                // Delivery capacity constraints
                constraint sum[c in 0...n](palettes[h][d][v][c] + 0.5 * (demi_palettes[h][d][v][c] + demi_palettes_s[h][d][v][c])) <= sizes[v] ;
                constraint sum[c in 0...n](load[h][d][v][c]) <= capacities[v] ;

                // Pickup capacity constraints
                constraint sum[p in 0...n_pdr](weight[p] * visits_r[h][d][v][p]) <= capacities[v];
                // We assume 2 palettes are used for each pickup
                constraint sum[p in 0...n_pdr](2 * visits_r[h][d][v][p]) <= sizes[v] ;


                // Deliver first, then pickup
                route_dist[h][d][v] <- 
                    + (liv ?  dist_from_depot[livraisons[h][d][v][0]] 
                            + sum(1...n_livraisons[h][d][v], i => dist_matrix[livraisons[h][d][v][i - 1]][livraisons[h][d][v][i]])
                            + (ram ? 
                                dist_matrix[livraisons[h][d][v][n_livraisons[h][d][v]-1]][ramasses[h][d][v][0] + n]
                                : dist_to_depot[livraisons[h][d][v][n_livraisons[h][d][v]-1]])

                            : (ram ? dist_from_depot[ramasses[h][d][v][0] + n] : 0)) 
                    + (ram ? dist_to_depot[ramasses[h][d][v][n_ramasses[h][d][v]-1] + n] : 0) 
                            + sum(1...n_ramasses[h][d][v], i => dist_matrix[ramasses[h][d][v][i - 1] + n][ramasses[h][d][v][i] + n])
                    ;

                // Constrain the estimated tour duration
                // Based on a linear proxy from the tour length
                // See fit_duration_proxy.ipynb
                tour_duration[h][d][v] <- 
                            params["duration_coefficients"][0] * 
                                    (route_dist[h][d][v] - (ram ? 
                                        dist_to_depot[ramasses[h][d][v][n_ramasses[h][d][v]-1] + n] : 
                                        (liv ? dist_to_depot[livraisons[h][d][v][n_livraisons[h][d][v]-1]] : 0))) + 
                            params["duration_coefficients"][1] + 
                            params["wait_at_centres"] * 60 * n_livraisons[h][d][v] +
                        params["wait_at_pdrs"] * 60 * n_ramasses[h][d][v]
                            ;

                constraint tour_duration[h][d][v] <= params["max_tour_duration"] * 60 ;
                constraint n_livraisons[h][d][v] + n_ramasses[h][d][v] <= params["max_stops"]; // This is basically redundant
            }

            constraint sum[c in 0...n][v in 0...m](norvegiennes[h][d][v][c]) <= params["n_norvegiennes"] ;

        }
    }

    // Time window constraints
    for [c in 0...n] {
        for [d in 0...n_days] {
            allowed[d][c] = false ;
            for [a in j_de_livraison[c]] {
                if (a==d) {
                    allowed[d][c] = true ;
                }
            }
            if (!allowed[d][c]) {
                // println("Disallowing ", centres.rows[index_to_centre[c]+1][0]," on day ", d) ;
                constraint sum[v in 0...m][h in 0...2](visits_l[h][d][v][c]) == 0 ;
            }
        }
    }

    // Pickup each site enough times each week
    for [p in 0...n_pdr] {
        for [d in 0...n_days] {

            // How many times do we have to pickup p on day d ?
            n_ram = 0 ;
            for [j in j_de_ramasse[p]] {
                if (j==d) {
                    n_ram += 1 ;
                }
            }
            constraint sum[v in 1...m : frais[v] == "Oui"][h in 0...2](visits_r[h][d][v][p]) == n_ram ; 
        }
    }

    // Meet demands for each product type
    for [c in 0...n] {
        constraint sum[h in 0...2][d in 0...n_days][v in 0...m](visits_l[h][d][v][c] * serve_a[h][d][v][c]) >= demands["a"][c];
        constraint sum[h in 0...2][d in 0...n_days][v in 0...m](visits_l[h][d][v][c] * serve_f[h][d][v][c]) >= demands["f"][c];
        constraint sum[h in 0...2][d in 0...n_days][v in 0...m](visits_l[h][d][v][c] * serve_s[h][d][v][c]) >= demands["s"][c];
    }

    add_specific_requirements() ;

    // total_distance <- sum[h in 0...2][d in 0...n_days][v in 0...m](route_dist[h][d][v]);
    fuel_consumption <- sum[h in 0...2][d in 0...n_days][v in 0...m](route_dist[h][d][v] * fuel[v] / 100000);

    used[v in 0...m] <- sum[h in 0...2][d in 0...n_days](n_livraisons[h][d][v] + n_ramasses[h][d][v]) > 0 ;
    minimize sum[v in 0...m](used[v]);
    // constraint sum[v in 0...m](used[v]) <= 3;

    // for [v in 3...9] {
    //     constraint used[v] <= used[v-1] ;
    // }
    

    minimize fuel_consumption ;

}

function main(args) {

    initfile=args[0] == "nil" ? nil : args[0];
    outfile=args[1];
    week=args[2].toInt();
    timelimit= args.count() >= 4 ? args[3].toInt() : nil ;

    lastBestValue = inf - 1 ;
    lastBestRunningTime = 0;
    lastSolutionWritten = false;

    input();

    with(ls = localsolver.create()) {

        model();

        ls.model.close();

        if (timelimit != nil) {
            ls.param.timeLimit = timelimit ;
        }

        ls.solve();
    }
}