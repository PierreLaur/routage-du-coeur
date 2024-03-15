use csv;
use json;
use io;
use localsolver;
use random;
use datetime;

function input() {

    problem = json.parse(problem_file) ;

    params = problem["params"] ;

    centres = problem["centres"];
    pdrs = problem["pdrs"];
    matrix = problem["distance_matrix"];
    matrix2 = problem["duration_matrix"];
    vehicles = problem["vehicles"];

    n_centres = problem["n_centres"];

    delivered_centres = {} ;
    delivered_centres_and_pdrs = {} ;
    semihebdo_index = 0 ;
    for [c in 1...n_centres] {
        semaine_livraison = centres[c]["delivery_week"] ;
        if(semaine_livraison == "ANY" || semaine_livraison == params["week"]) {
            delivered_centres.add(c) ;
            delivered_centres_and_pdrs.add(c) ;
        }
    }

    // The number of centres to deliver this week, not including the depot
    n = delivered_centres.count();

    m = vehicles.count();
    n_pdr = pdrs.count();
    n_cp = n + n_pdr;
    n_days = 5 ;

    for [p in n_centres...n_centres+n_pdr] {
        delivered_centres_and_pdrs.add(p) ;
    }

    product_types = {"A", "F", "S"};

    // Demands
    for [i in delivered_centres] {
        demands["A"][i] = centres[i]["demands"]["A"] ;
        pallet_demands[i]["A"] = ceil(demands["A"][i] / params["max_palette_capacity"]) ;

        demands["F"][i] = centres[i]["demands"]["A"] ;
        pallet_demands[i]["F"] = ceil(demands["F"][i] / params["demi_palette_capacity"]) ;
        
        demands["S"][i] = centres[i]["demands"]["A"] ;
        pallet_demands[i]["S"] = ceil(demands["S"][i] / params["demi_palette_capacity"]) ;
    }

    realindex = {} ;
    dupnodes = {} ;
    nodetype = {} ;

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
            for [j in pdrs[cp - n_centres]["required_days"]] {
                dupnodes[cp].add(i) ;
                realindex[i] = cp ;
                nodetype[i] = pdrs[cp - n_centres]["product_type"] ;
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

            dist_from_depot[node] = matrix[0][cp] ;
            time_from_depot[node] = matrix2[0][cp] ;
            dist_to_depot[node] = matrix[cp][0] ;
            time_to_depot[node] = matrix2[cp][0] ;

            for [node2 in 0...n_nodes] {
                dist_matrix[node][node2] = matrix[cp][realindex[node2]];
                time_matrix[node][node2] = matrix2[cp][realindex[node2]];
            }
        }
    }

    allowed_vehicles = {};
    for [v in 0...m] {

        if (vehicles[v]["allowed"]) {
            allowed_vehicles.add(v);
        }

        // Whether the vehicle can transport each type of product
        for [t in product_types] {
            can_carry[v][t] = false;
            for [t2 in vehicles[v]["can_carry"]] {
                if (t == t2) {
                    can_carry[v][t] = true ;
                    break;
                }
            }
        }
    }

    // Can we deliver c on day d ?
    for [c in delivered_centres] {
        for [d in 0...n_days] {
            delivery_allowed[d][c] = false ;
            for [a in centres[c]["allowed_days"]] {
                if (a==d) {
                    delivery_allowed[d][c] = true ;
                    break ;
                }
            }
        }
    }

    // Do we have to pickup p on day d ?
    for [p in 0...n_pdr] {
        for [d in 0...n_days] {
            pickup_d_p[d][p] = false ;
            for [j in pdrs[p]["required_days"]] {
                if (j==d) {
                    pickup_d_p[d][p] = true ;
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
        constraint and[d in 0...n_days][v in allowed_vehicles : (d!=1 || v!=2)]
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

    livraisons[d in 0...n_days][v in allowed_vehicles] <- list(n_centre_nodes) ;
    ramasses[d in 0...n_days][v in allowed_vehicles] <- list(n_pdr_nodes);

    // Every node has to be delivered/picked up exactly once
    constraint partition[d in 0...n_days][v in allowed_vehicles](livraisons[d][v]);
    constraint partition[d in 0...n_days][v in allowed_vehicles](ramasses[d][v]);

    for [d in 0...n_days] {
        for [v in allowed_vehicles] {
            n_livraisons[d][v] <- count(livraisons[d][v]) ;
            n_ramasses[d][v] <- count(ramasses[d][v]) ;

            liv[d][v] <- n_livraisons[d][v] > 0 ;
            ram[d][v] <- n_ramasses[d][v] > 0 ;

            visits_l[d][v][c in 0...n_centre_nodes] <- contains(livraisons[d][v], c) ;
            visits_r[d][v][p in 0...n_pdr_nodes] <- contains(ramasses[d][v], p);

            for [t in product_types] {
                if (can_carry[v][t]) continue ;
                for [c in 0...n_centre_nodes : nodetype[c] == t] {
                    constraint !visits_l[d][v][c] ;
                }
                for [p in 0...n_pdr_nodes : nodetype[p + n_centre_nodes] == t] {
                    constraint !visits_r[d][v][p] ;
                }
            }

            for [c in 0...n_centre_nodes : !delivery_allowed[d][realindex[c]]] {
                constraint !visits_l[d][v][c] ;
            }

            for [p in 0...n_pdr_nodes : !pickup_d_p[d][realindex[p + n_centre_nodes] - n_centres]] {
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
            constraint sum[c in 0...n_centre_nodes](visits_l[d][v][c] * quantity[d][v][c]) <= vehicles[v]["capacity"] ;
            constraint sum[c in 0...n_centre_nodes](
                    2 * visits_l[d][v][c] * (nodetype[c] == "A") + 
                    visits_l[d][v][c] * (nodetype[c] == "F") + 
                    demi_palettes_s[d][v][c] * (nodetype[c] == "S")) <= 2 * vehicles[v]["size"] ;

            // Pickup capacity constraints
            constraint sum[p in 0...n_pdr_nodes](pdrs[realindex[p+n_centre_nodes] - n_centres]["weight"] * visits_r[d][v][p]) <= vehicles[v]["capacity"];
            constraint sum[p in 0...n_pdr_nodes](2 * visits_r[d][v][p]) <= vehicles[v]["size"] ; 

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
        //             [v in allowed_vehicles]
        //             (norvegiennes[d][v][c]) <= params["n_norvegiennes"] ;
    }

    add_specific_requirements() ;

    // Meet demands for each product type
    for [c in delivered_centres] {
        for [type in product_types] {
            constraint sum[d in 0...n_days : delivery_allowed[d][c]]
                            [v in allowed_vehicles]
                            [node in dupnodes[c] : nodetype[node]==type]
                            (visits_l[d][v][node] * quantity[d][v][node]) >= demands[type][c];
        }
    }

    used[v in allowed_vehicles] <- or[d in 0...n_days](liv[d][v] || ram[d][v]) ;
    fixed_costs <- sum[v in allowed_vehicles](vehicles[v]["fixed_cost"] * used[v]);
    n_used <- sum[v in allowed_vehicles](used[v]) ;

    variable_costs <- sum[d in 0...n_days][v in allowed_vehicles](vehicles[v]["cost_per_km"] * route_dist[d][v] / 1000);

    total_costs <- variable_costs + fixed_costs ;
    total_distance <- sum[v in allowed_vehicles][d in 0...n_days](route_dist[d][v]);

    // constraint total_costs <= 455 ;

    minimize total_costs ;

}

function write_solution() {
    if (outfile == nil) return;
    local outf = io.openWrite(outfile);

    sol = {};
    sol["total_distance"] = sum[d in 0...n_days][v in allowed_vehicles](route_dist[d][v].value);
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
        for [v in allowed_vehicles] {
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

    usage = "Usage: ls_solver_split problem_file";
    if (args.count() < 1) {
        throw usage;
    }

    problem_file = args[0];

    input();

    lastBestValue = inf - 1 ;
    lastBestRunningTime = 0;
    lastSolutionWritten = false;

    with(ls = localsolver.create()) {
        ls.param.verbosity = 2 ;

        model();

        ls.model.close();

        ls.solve();
    }
}