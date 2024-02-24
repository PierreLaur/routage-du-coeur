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

    n_centres = centres.nbRows ;

    delivered_centres = {} ;
    delivered_centres_and_pdrs = {} ;
    semihebdo_index = 0 ;
    for [c in 1...n_centres] {
        semaine_livraison = centres.rows[c][7] ;

        if (semaine_livraison != 0) {
            semaine_livraison = params["week_assignments"][semihebdo_index] ;
            semihebdo_index += 1 ;
        }

        if(semaine_livraison == 0 || semaine_livraison == week) {
            delivered_centres.add(c) ;
            delivered_centres_and_pdrs.add(c) ;
        }
    }

    // The number of centres to deliver this week, not including the depot
    n = delivered_centres.count();

    m = vehicles.nbRows;
    n_pdr = pdr.nbRows;
    n_cp = n + n_pdr;
    n_days = 5 ;

    for [p in n_centres...n_centres+n_pdr] {
        delivered_centres_and_pdrs.add(p) ;
    }

    // Demands
    for [i in delivered_centres] {
        demands["A"][i] = ceil(centres.rows[i][4] * params["robustness_factor"]) ;
        pallet_demands[i]["A"] = ceil(demands["A"][i] / params["max_palette_capacity"]) ;

        demands["F"][i] = ceil(centres.rows[i][5] * params["robustness_factor"]) ;
        pallet_demands[i]["F"] = ceil(demands["F"][i] / params["demi_palette_capacity"]) ;
        
        demands["S"][i] = ceil(centres.rows[i][6] * params["robustness_factor"]) ;
        pallet_demands[i]["S"] = ceil(demands["S"][i] / params["demi_palette_capacity"]) ;
    }

    // Pickup weights (per pickup)
    for [p in n_centres...n_centres + n_pdr] {
        weight[p] = ceil(pdr.rows[p - n_centres][5]);
    }

    realindex = {} ;
    dupnodes = {} ;
    nodetype = {} ;

    jours_map = {"Lundi": 0, "Mardi": 1, "Mercredi": 2, "Jeudi": 3, "Vendredi": 4};
    // When do we have to pickup each site ?
    for [i in 0...n_pdr] {
        jours = pdr.rows[i][4].split(", ") ;
        j_de_ramasse[i] = {};
        for [k in 0...jours.count()] {
            j_de_ramasse[i].add(jours_map[jours[k]]);
        }
    }

    // Duplicate the nodes for each palette and store a map
    i = 0 ;
    for [cp in delivered_centres_and_pdrs] {
        dupnodes[cp] = {} ;

        if (cp < n_centres) {
            for [type in {"A", "F", "S"}] {
                for [j in 0...pallet_demands[cp][type]] {
                    dupnodes[cp].add(i) ;
                    realindex[i] = cp ;
                    nodetype[i] = type ;
                    i += 1 ;
                }
            }
        } else {
            for [j in j_de_ramasse[cp - n_centres]] {
                dupnodes[cp].add(i) ;
                realindex[i] = cp ;
                nodetype[i] = pdr.rows[cp - n_centres][6] ;
                i += 1 ;
            }
        }

        // n_centre_nodes will be nil if we don't deliver the last centre
        if (cp == n_centres - 1) {
            n_centre_nodes = i ;
        }
    }

    n_nodes = i ;
    n_pdr_nodes = n_nodes - n_centre_nodes ;

    // Write down the distance matrix for nodes
    for [cp in delivered_centres_and_pdrs] {

        for [node in dupnodes[cp]] {
            dist_from_depot[node] = matrix.rows[0][cp+1] ;
            dist_to_depot[node] = matrix.rows[cp][1] ;
            time_from_depot[node] = matrix2.rows[0][cp+1] ;
            time_to_depot[node] = matrix2.rows[cp][1] ;

            for [node2 in 0...n_nodes] {
                dist_matrix[node][node2] = matrix.rows[cp][realindex[node2]+1];
                time_matrix[node][node2] = matrix2.rows[cp][realindex[node2]+1];
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

    // When are we allowed to deliver each centre ?
    for [c in delivered_centres] {
        jours = centres.rows[c][10].replace(" ", "").split(",") ;
        j_de_livraison[c] = {};
        for [k in 0...jours.count()] {
            j_de_livraison[c].add(jours_map[jours[k]]);
        }
    }
    for [c in delivered_centres] {
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

    for [p in 0...n_pdr] {
        i = 0 ;
        for [d in 0...n_days] {
            for [node in dupnodes[p + n_centres]] {
                pickup_d_p[d][node] = false ;
            }

            // Do we have to pickup p on day d ?
            for [j in j_de_ramasse[p]] {
                if (j==d) {
                    pickup_d_p[d][dupnodes[p + n_centres][i]] = true ;
                    i += 1 ;
                    break ;
                }
            }
        }
    }
}


function add_specific_requirements() {
    /* 
    Adds a few "hard-coded" constraints 
    These are not written in the input files but have been specified by the client 
    */

    cca_nodes = dupnodes[n_centres + 6] ;
    cc_nodes = dupnodes[n_centres + 5] ;

    // The PL must visit Carrefour Centrale on Wednesday
    constraint visits_r[2][0][cca_nodes[0] - n_centre_nodes] ;
    // He must only do one pickup
    constraint n_ramasses[2][0] == 1;

    // The first camion Frigo must visit Carrefour Centrale on Wednesday and Friday
    constraint visits_r[2][2][cc_nodes[0] - n_centre_nodes] ;
    constraint visits_r[4][2][cc_nodes[1] - n_centre_nodes] ;

    // // To simplify we say it only does one pickup (in reality he delivers a centre afterwards)
    constraint n_stops_r[2][2] == 1; 
    constraint n_stops_r[4][2] == 1;

    // The PL cannot deliver Toulouse/Grande-Bretagne
    index_gde_bretagne = 24 ;
    for [node in dupnodes[index_gde_bretagne]]{
        constraint and[d in 0...n_days](!visits_l[d][0][node]);
    }

    // A camion frigo must deliver Revel on Tuesdays and do nothing else
    index_revel = 20 ;
    for [node in dupnodes[index_revel]] {
        constraint visits_l[1][2][node] ;
        constraint n_stops_l[1][2] == 1 ; 
        constraint n_stops_r[1][2] == 0 ;
        constraint and[d in 0...n_days][v in 0...m : params["vehicle_allowed"][v] && (d!=1 || v!=2)]
                    (!visits_l[d][v][node]) ;
    }

    // A Camion Frigo must deliver les arènes on Friday then pickup Leclerc Blagnac
    index_arenes = 26 ;
    constraint or[node in dupnodes[index_arenes]](visits_l[4][3][node]) ;

    index_lc_blagnac = 3 + n_centres;
    nodes_blagnac = dupnodes[index_lc_blagnac] ;
    constraint visits_r[4][3][nodes_blagnac[nodes_blagnac.count()-1] - n_centre_nodes];
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

            if (frais[v] != "Oui") {
                for [c in 0...n_centre_nodes : nodetype[c] != "A"] {
                    constraint !visits_l[d][v][c] ;
                }
                for [p in 0...n_pdr_nodes : nodetype[p + n_centre_nodes] == "F"] {
                    constraint !visits_r[d][v][p] ;
                }
            }

            for [c in 0...n_centre_nodes : !delivery_allowed[d][realindex[c]]] {
                constraint !visits_l[d][v][c] ;
            }

            for [p in 0...n_pdr_nodes : !pickup_d_p[d][p+n_centre_nodes]] {
                constraint !visits_r[d][v][p] ;
            }

            quantity[d][v][c in 0...n_centre_nodes] <- int(0, min(
                    nodetype[c] == "A" ? params["max_palette_capacity"] : params["demi_palette_capacity"],
                    demands[nodetype[c]][realindex[c]]
                    )
                ) ;

            demi_palettes_s[d][v][c in 0...n_centre_nodes] <- int(0, 1) ;
            norvegiennes[d][v][c in 0...n_centre_nodes] <- nodetype[c] == "S" ? max(0, ceil(
                                            (visits_l[d][v][c] * quantity[d][v][c] - demi_palettes_s[d][v][c] * params["demi_palette_capacity"]) 
                                            / params["norvegienne_capacity"] 
                                            )) : 0;

            for [c in 0...n_centre_nodes : nodetype[c] == "S"] {
                constraint visits_l[d][v][c] * quantity[d][v][c] <= 
                            norvegiennes[d][v][c] * params["norvegienne_capacity"] +
                            demi_palettes_s[d][v][c] * params["demi_palette_capacity"] ;
            }


            // Delivery capacity constraints
            constraint sum[c in 0...n_centre_nodes](visits_l[d][v][c] * quantity[d][v][c]) <= capacities[v] ;
            constraint sum[c in 0...n_centre_nodes](
                    2 * visits_l[d][v][c] * (nodetype[c] == "A") + 
                    visits_l[d][v][c] * (nodetype[c] == "F") + 
                    demi_palettes_s[d][v][c] * (nodetype[c] == "S")) <= 2 * sizes[v] ;

            // Pickup capacity constraints
            constraint sum[p in 0...n_pdr_nodes](weight[realindex[p+n_centre_nodes]] * visits_r[d][v][p]) <= capacities[v];
            constraint sum[p in 0...n_pdr_nodes](2 * visits_r[d][v][p]) <= sizes[v] ; 

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

            n_stops_l[d][v] <- liv[d][v] ? 
                            count(distinct(0...n_livraisons[d][v], i => realindex[livraisons[d][v][i]]))
                            : 0 ;
            n_stops_r[d][v] <- ram[d][v] ? 
                            count(distinct(0...n_ramasses[d][v], i => realindex[ramasses[d][v][i] + n_centre_nodes]))
                            : 0 ;

            // Constrain the tour duration
            local end_of_delivery_time <- (liv[d][v] ? time_from_depot[livraisons[d][v][0]] 
                                        + sum(1...n_livraisons[d][v], i => time_matrix[livraisons[d][v][i - 1]][livraisons[d][v][i]]) 
                                            + n_stops_l[d][v] * params["wait_at_centres"] * 60
                                            + (ram[d][v] ? time_matrix[livraisons[d][v][n_livraisons[d][v]-1]][ramasses[d][v][0] + n_centre_nodes] : 0)
                                : (ram[d][v] ? time_from_depot[ramasses[d][v][0] + n_centre_nodes] : 0))
                    ;


            tour_duration[d][v] <- end_of_delivery_time
                    + (ram[d][v] ? 
                            sum(1...n_ramasses[d][v], i => time_matrix[ramasses[d][v][i - 1] + n_centre_nodes][ramasses[d][v][i] + n_centre_nodes]) 
                            + n_stops_r[d][v] * params["wait_at_pdrs"] * 60
                        : 0)
                    ;

            constraint tour_duration[d][v] <= params["max_tour_duration"] * 60 ;


            constraint n_stops_l[d][v] + n_stops_r[d][v] <= params["max_stops"] ;
        }

        // constraint sum[c in 0...n : delivery_allowed[d][c]]
        //             [v in 0...m : params["vehicle_allowed"][v]]
        //             (norvegiennes[d][v][c]) <= params["n_norvegiennes"] ;
    }

    add_specific_requirements() ;

    // Meet demands for each product type
    for [c in delivered_centres] {
        for [type in {"A", "F", "S"}] {
            constraint sum[d in 0...n_days : delivery_allowed[d][c]]
                            [v in 0...m : params["vehicle_allowed"][v]]
                            [node in dupnodes[c] : nodetype[node]==type]
                            (visits_l[d][v][node] * quantity[d][v][node]) >= demands[type][c];
        }
    }

    fuel_consumption <- sum[d in 0...n_days][v in 0...m : params["vehicle_allowed"][v]](route_dist[d][v] * fuel[v] / 100000);
    used[v in 0...m : params["vehicle_allowed"][v]] <- or[d in 0...n_days](liv[d][v] || ram[d][v]) ;

    n_used <- sum[v in 0...m : params["vehicle_allowed"][v]](used[v]) ;
    total_cost <- fuel_consumption * params["fuel_cost"] + 15 * n_used ;

    minimize total_cost ;

}

function write_solution() {
    if (outfile == nil) return;
    local outf = io.openWrite(outfile);

    sol = {};
    sol["total_distance"] = sum[d in 0...n_days][v in 0...m : params["vehicle_allowed"][v]](route_dist[d][v].value);
    sol["fuel_consumption"] = fuel_consumption.value;
    sol["total_cost"] = total_cost.value;
    sol["vehicles_used"] = {} ;
    for [v in 0...m] {
        if (params["vehicle_allowed"][v]) {
            sol["vehicles_used"].add(used[v].value) ;
        } else {
            sol["vehicles_used"].add(0) ;
        }
    }
    sol["tour_durations"] = {} ;

    sol["tours"] = {};

    for [d in 0...n_days] {
        for [v in 0...m : params["vehicle_allowed"][v]] {
            if (livraisons[d][v].value.count() == 0 && ramasses[d][v].value.count() == 0) continue;

            key = d + ", " + v ;

            sol["tour_durations"][key] = tour_duration[d][v].value ;
            if (n_ramasses[d][v].value > 0) {
                sol["tour_durations"][key] += time_to_depot[ramasses[d][v].value[n_ramasses[d][v].value-1]] ;
            }

            sol["tours"][key] = {};

            i = 0 ;
            pals = {0, 0, 0} ;
            deliv = {0, 0, 0} ;
            norv = 0 ;
            typemap = {"A":0, "F":1, "S":2} ;
            while (i < n_livraisons[d][v].value) {
                node = livraisons[d][v].value[i] ;
                type = typemap[nodetype[node]] ;
                if (i==0 || realindex[node] == realindex[prev_node]) {
                } else {
                    sol["tours"][key].add({
                        "index":realindex[prev_node],
                        "type":"livraison",
                        "name":matrix.rows[realindex[prev_node]][0],
                        "delivery":deliv,
                        "palettes":pals,
                        "norvegiennes":norv
                    });
                    pals = {0, 0, 0} ;
                    deliv = {0, 0, 0} ;
                    norv = 0 ;
                }
                if (type == 2) {
                    pals[type] += demi_palettes_s[d][v][node].value/2 ;
                } else if (type == 1) {
                    pals[type] += 0.5 ;
                } else {
                    pals[type] += 1 ;
                }
                deliv[type] += quantity[d][v][node].value ;
                norv += norvegiennes[d][v][node].value ;
                prev_node = node ;
                i += 1 ;
            }
            if (n_livraisons[d][v].value > 0) {
                sol["tours"][key].add({
                        "index":realindex[prev_node],
                        "type":"livraison",
                        "name":matrix.rows[realindex[prev_node]][0],
                        "delivery":deliv,
                        "palettes":pals,
                        "norvegiennes":norv
                    });
            }

            i = 0 ;
            while (i < n_ramasses[d][v].value) {
                node = ramasses[d][v].value[i] + n_centre_nodes ;
                if (i==0 || realindex[node] == realindex[prev_node]) {
                } else {
                    sol["tours"][key].add({
                        "index":realindex[prev_node],
                        "type":"ramasse",
                        "name":matrix.rows[realindex[prev_node]][0]
                    });
                }
                prev_node = node ;
                i += 1 ;
            }
            if (n_ramasses[d][v].value > 0) {
                sol["tours"][key].add({
                        "index":realindex[prev_node],
                        "type":"ramasse",
                        "name":matrix.rows[realindex[prev_node]][0]
                    });
            }
        }
    }
    json.dump(sol, outf);
}

function callback(ls, cbType) {
    local stats = ls.statistics;
    local obj <- ls.model.objectives[0];
    if (ls.solution.status == "FEASIBLE" && obj.value < lastBestValue) {
        lastBestRunningTime = stats.runningTime;
        lastBestValue = obj.value;
        lastSolutionWritten = false ;

        println(
            "[   ", stats.runningTime, "s] : ", total_cost.value, "E    ", fuel_consumption.value, "L   ", n_used.value, " vehicles"
        );
    }

    if (stats.runningTime - lastBestRunningTime > 2 && !lastSolutionWritten) {
        println(">>>>>>> No improvement during 5 seconds: writing current solution to file");
        write_solution();
        lastSolutionWritten = true ;
    }
}

function main(args) {

    initfile=args[0] == "nil" ? nil : args[0];
    outfile=args[1];
    week=args[2].toInt();

    lastBestValue = inf - 1 ;
    lastBestRunningTime = 0;
    lastSolutionWritten = false;

    input();

    with(ls = localsolver.create()) {
        ls.addCallback("ITERATION_TICKED", callback);
        ls.param.verbosity = 0 ;

        model();

        ls.param.setAdvancedParam("clarkeAndWrightEnabled", false);
        // ls.param.timeLimit = 10;
        ls.model.close();

        if (initfile != nil) {
            set_initial_solution(false, nil);
        }

        ls.solve();
    }
}