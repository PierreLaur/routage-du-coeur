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
    
    centres = csv.parse("data/centres_keep.csv");
    pdr = csv.parse("data/points_de_ramasse.csv");
    matrix = csv.parse("data/euclidean_matrix.csv");
    matrix2 = csv.parse("data/duration_matrix_w_traffic.csv");
    vehicles = csv.parse("data/vehicules.csv");

    // The number of centres to serve this week
    n = centres.nbRows - 1 ;
    m = vehicles.nbRows;
    n_pdr = pdr.nbRows;
    n_cp = n + n_pdr;
    n_days = 5 ;

    // Demands
    for [i in 0...n] {
        demand_kg = ceil(centres.rows[i+1][4] * params["robustness_factor"]) ;
        demands[i]["a"] = ceil(demand_kg / params["demi_palette_capacity"]) ;

        demand_kg = ceil(centres.rows[i+1][5] * params["robustness_factor"]) ;
        demands[i]["f"] = ceil(demand_kg / params["demi_palette_capacity"]) ;
        
        demand_kg = ceil(centres.rows[i+1][6] * params["robustness_factor"]) ;
        demands[i]["s"] = ceil(demand_kg / params["demi_palette_capacity"]) ;
    }

    // Pickup weights (per pickup)
    for [p in 0...n_pdr] {
        demands[n + p]["a"] = 0;
        demands[n + p]["f"] = 2;
        demands[n + p]["s"] = 0;
    }

    realindex = {} ;
    dupnodes = {} ;

    // Duplicate the nodes for each palette and store a map
    i = 0 ;
    for [cp in 0...n_cp] {
        dupnodes[cp] = {} ;
        for [type in {"a", "f", "s"}] {
            for [j in 0...demands[cp][type]] {
                dupnodes[cp].add(i) ;
                realindex[i] = cp ;
                i += 1 ;
            }
        }
        if (cp == n - 1) {
            n_centre_nodes = i ;
        }
    }
    n_nodes = i ;
    n_pdr_nodes = n_nodes - n_centre_nodes ;

    // Write down the distance matrix for nodes
    for [cp in 0...n_cp] {

        for [node in dupnodes[cp]] {
            dist_from_depot[node] = matrix.rows[0][cp+2] ;
            dist_to_depot[node] = matrix.rows[cp+1][1] ;
            time_from_depot[node] = matrix2.rows[0][cp+2] ;
            time_to_depot[node] = matrix2.rows[cp+1][1] ;

            for [node2 in 0...n_nodes] {
                dist_matrix[node][node2] = matrix.rows[cp+1][realindex[node2]+2];
                time_matrix[node][node2] = matrix2.rows[cp+1][realindex[node2]+2];
            }
        }
    }

    for [i in 0...m] {
        row = vehicles.rows[i];
        // Capacities in kg
        capacities[i] = row[1];
        
        // Sizes in number of pallets
        sizes[i] = row[2];

        // Fuel consumption (L/100km)
        fuel[i] = row[3];

        // Whether the vehicle can transport products that need to be refrigerated
        frais[i] = row[4];
    }

    jours_map = {"Lundi": 0, "Mardi": 1, "Mercredi": 2, "Jeudi": 3, "Vendredi": 4};
    // When are we allowed to deliver each centre ?
    for [c in 0...n] {
        index = c+1 ;
        jours = centres.rows[index][10].replace(" ", "").split(",") ;
        j_de_livraison[c] = {};
        for [k in 0...jours.count()] {
            j_de_livraison[c].add(jours_map[jours[k]]);
        }
    }
    for [c in 0...n] {
        for [d in 0...n_days] {
            delivery_allowed[d][c] = false ;
            for [a in j_de_livraison[c]] {
                if (a==d) {
                    delivery_allowed[d][c] = true ;
                    break ;
                }
            }
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
    constraint visits_r[2][0][5] ;
    // He must only do one pickup
    constraint n_ramasses[2][0] == 1;

    // The first camion Frigo must visit Carrefour Centrale on Wednesday and Friday
    constraint visits_r[2][2][5] ;
    constraint visits_r[4][2][5] ;
    // To simplify we say it only does one pickup (in reality he delivers a centre afterwards)
    constraint n_ramasses[2][2] == 1;
    constraint n_ramasses[4][2] == 1;

    // The PL cannot deliver Toulouse/Grande-Bretagne
    index_gde_bretagne = 24 - 1 ;
    for [c in 0...n] {
        if (index_to_centre[c] == index_gde_bretagne) {
            constraint and[d in 0...n_days](!visits_l[d][0][c]);
            break ;
        }
    }

    // A camion frigo must deliver Revel on Tuesdays and do nothing else
    index_revel = 20 - 1 ;
    for [c in 0...n : index_to_centre[c] == index_revel] {
        constraint visits_l[1][2][c] ;
        constraint n_livraisons[1][2] == 1 ;
        constraint n_ramasses[1][2] == 0 ;
        constraint and[d in 0...n_days][v in 0...m : params["vehicle_allowed"][v] && (d!=1 || v!=2)]
                    (!visits_l[d][v][c]) ;
    }

    // A Camion Frigo must deliver les arènes on Friday then pickup Leclerc Blagnac
    index_lc_blagnac = 3 ;
    index_arenes = 26 - 1 ;
    for [c in 0...n : index_to_centre[c] == index_arenes] {
        constraint visits_l[4][3][c] ;
        constraint visits_r[4][3][index_lc_blagnac];
    }
}

function model() {

    livraisons[d in 0...n_days][v in 0...m : params["vehicle_allowed"][v]] <- list(n_centre_nodes) ;
    ramasses[d in 0...n_days][v in 0...m : params["vehicle_allowed"][v]] <- list(n_pdr_nodes);

    // Every node has to be delivered/picked up exactly once
    constraint partition[d in 0...n_days][v in 0...m : params["vehicle_allowed"][v]](livraisons[d][v]);
    constraint partition[d in 0...n_days][v in 0...m : params["vehicle_allowed"][v]](ramasses[d][v]);

    for [d in 0...n_days] {
        for [v in 0...m : params["vehicle_allowed"][v]] {
            n_livraisons[d][v] <- count(livraisons[d][v]) ;
            n_ramasses[d][v] <- count(ramasses[d][v]) ;

            liv[d][v] <- n_livraisons[d][v] > 0 ;
            ram[d][v] <- n_ramasses[d][v] > 0 ;

            visits_l[d][v][c in 0...n_centre_nodes] <- contains(livraisons[d][v], c) ;
            visits_r[d][v][p in 0...n_pdr_nodes] <- contains(ramasses[d][v], p);

            // Delivery capacity constraints
            constraint sum[c in 0...n_centre_nodes](visits_l[d][v][c]) <= 2 * sizes[v] ;
            constraint sum[c in 0...n_centre_nodes](visits_l[d][v][c] * params["demi_palette_capacity"]) <= capacities[v] ;

            // Pickup capacity constraints
            constraint sum[p in 0...n_pdr](params["max_palette_capacity"] * visits_r[d][v][p]) <= capacities[v];
            constraint sum[p in 0...n_pdr_nodes](visits_r[d][v][p]) <= sizes[v] ;

            // Deliver first, then pickup
            route_dist[d][v] <- (liv[d][v] ? dist_from_depot[livraisons[d][v][0]]
                        + sum(1...n_livraisons[d][v], i => dist_matrix[livraisons[d][v][i - 1]][livraisons[d][v][i]])
                        + (ram[d][v] ? 
                            dist_matrix[livraisons[d][v][n_livraisons[d][v]-1]][ramasses[d][v][0] + n_centre_nodes]
                            : dist_to_depot[livraisons[d][v][n_livraisons[d][v]-1]])

                        : (ram[d][v] ? dist_from_depot[ramasses[d][v][0] + n_centre_nodes] : 0)) 
                    + (ram[d][v] ? dist_to_depot[ramasses[d][v][n_ramasses[d][v]-1] + n_centre_nodes]  
                        + sum(1...n_ramasses[d][v], i => dist_matrix[ramasses[d][v][i - 1] + n_centre_nodes][ramasses[d][v][i] + n_centre_nodes]): 0)
                    ;

            // // Constrain the tour duration
            // local end_of_delivery_time <- (liv[d][v] ? 
            //                                     time_from_depot[livraisons[d][v][n_livraisons[d][v]-1]] 
            //                                     + sum(1...n_livraisons[d][v], 
            //                                         i => time_matrix[livraisons[d][v][i - 1]][livraisons[d][v][i]] 
            //                                         + (realindex[i-1]!=realindex[i] ? params["wait_at_centres"] * 60 : 0)
            //                                         )
            //                                     + (ram[d][v] ? time_matrix[livraisons[d][v][n_livraisons[d][v]-1]][ramasses[d][v][0] + n_centre_nodes] : 0)
            //             : (ram[d][v] ? time_from_depot[ramasses[d][v][0] + n_centre_nodes] : 0)) 
            //         ;


            // tour_duration[d][v] <- end_of_delivery_time
            //         + (ram[d][v] ? sum(1...n_ramasses[d][v], 
            //                     i => time_matrix[ramasses[d][v][i - 1] + n][ramasses[d][v][i] + n] +
            //                     (realindex[i-1]!=realindex[i] ? params["wait_at_pdrs"] * 60 : 0))
            //                     : 0)
            //         ;

            // tour_duration[d][v] <- end_of_delivery_time
            //         + (ram[d][v] ? sum(1...n_ramasses[d][v], i => time_matrix[ramasses[d][v][i - 1] + n][ramasses[d][v][i] + n]) + params["wait_at_pdrs"] * 60 * n_ramasses[d][v] : 0)
            //         ;



            local seq <- livraisons[d][v] ;
            local c <- count(livraisons[d][v]) ;
           
            // Minimal segfault example
            // (c > 0 ? dist_from_depot[seq[1]] : 0) works
            // (c > 0 ? dist_from_depot[seq[c-2]] : 0) works
            // (c > 0 ? dist_from_depot[seq[c-1]] : 0) produces segfault
            // (c > 0 ? dist_from_depot[seq[0]] : 0) produces segfault here but is used line 208 without any problems
            tour_duration[d][v] <- (c > 0 ? dist_from_depot[seq[0]] : 0) ;

            constraint tour_duration[d][v] <= params["max_tour_duration"] * 60 ;

            // constraint n_livraisons[d][v] + n_ramasses[d][v] <= params["max_stops"] ;
        }

        // constraint sum[c in 0...n : delivery_allowed[d][c]]
        //             [v in 0...m : params["vehicle_allowed"][v]]
        //             (norvegiennes[d][v][c]) <= params["n_norvegiennes"] ;
    }

    // add_specific_requirements() ;

    fuel_consumption <- sum[d in 0...n_days][v in 0...m : params["vehicle_allowed"][v]](route_dist[d][v] * fuel[v] / 100000);
    used[v in 0...m : params["vehicle_allowed"][v]] <- or[d in 0...n_days](liv[d][v] || ram[d][v]) ;

    n_used <- sum[v in 0...m : params["vehicle_allowed"][v]](used[v]) ;
    total_cost <- fuel_consumption * params["fuel_cost"] + params["weekly_fixed_cost"] * n_used ;

    minimize total_cost ;

}

function main(args) {

    initfile=args[0] == "nil" ? nil : args[0];
    outfile=args[1];
    week=args[2].toInt();


    input();

    with(ls = localsolver.create()) {
        ls.param.verbosity = 2 ;

        model();

        ls.model.close();

        if (initfile != nil) {
            set_initial_solution(false, nil);
        }

        ls.solve();
    }
}