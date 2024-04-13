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
    livraisons_de_ramasses = problem["livraisons_de_ramasses"];

    n_centres = problem["n_centres"];

    // list variables require 0-based dense indices 
    // drop the depot & any centres that are not delivered this week to achieve that
    // create a map to be able to link new indices to indices in the problem data
    index_to_centre = {} ;
    i = 0 ;
    for [c in 1...n_centres] {
        semaine_livraison = centres[c]["delivery_week"] ;
        if(semaine_livraison == "ANY" || semaine_livraison == params["week"]) {
            index_to_centre[i] = c ;
            i += 1 ;
        }
    }

    // The number of centres to serve this week
    n = index_to_centre.count();

    m = vehicles.count();
    n_pdr = pdrs.count();
    n_nodes = n + n_pdr;
    n_days = 5;
    n_trips = params["max_trips"];

    product_types = {"A", "F", "S"};

    // Adjust the distance & time matrices to the new indices
    for [node in 0...n_nodes] {
        node_index = node < n ? index_to_centre[node] : node - n + n_centres ;

        dist_from_depot[node] = matrix[0][node_index] ;
        time_from_depot[node] = matrix2[0][node_index] ;
        dist_to_depot[node] = matrix[node_index][0] ;
        time_to_depot[node] = matrix2[node_index][0] ;

        for [node2 in 0...n_nodes] {
            node_index2 = node2 < n ? index_to_centre[node2] : node2 - n + n_centres ;
            
            dist_matrix[node][node2] = matrix[node_index][node_index2];
            time_matrix[node][node2] = matrix2[node_index][node_index2];
        }
    }

    // Demands
    n_scenarios = centres[0]["demands"].count() ;
    for [c in 0...n] {
        centre = index_to_centre[c] ;

        for [s in 0...n_scenarios] {
            for [t in product_types] {
                demands[s][c][t] = centres[centre]["demands"][s][t];
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
    for [c in 0...n] {
        for [d in 0...n_days] {
            delivery_allowed[d][c] = false ;
            for [a in centres[index_to_centre[c]]["allowed_days"]] {
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

    for [d in 0...n_days] {
        for [v in 0...m] {
            for[trip in 0...n_trips] {
                for [c in 0...n] {
                    is_ldr[d][v][trip][c] = false ;
                }
            }
        }
    }

    for [dvtc in livraisons_de_ramasses] {
        d = dvtc[0] ;
        v = dvtc[1] ;
        trip = dvtc[2] ;
        c = dvtc[3] ;
        for [c2 in 0...n : index_to_centre[c2] == c] {
            is_ldr[d][v][trip][c2] = true ;
        }
    }


}


function add_specific_requirements() {
    /* 
    Adds a few "hard-coded" constraints 
    These are not written in the input files but have been specified by the client 
    */

    // The PL must visit Carrefour Centrale on Wednesday
    constraint visits_r[2][0][0][6] ;
    // He must only do one pickup
    constraint n_ramasses[2][0][0] == 1;

    // The first camion Frigo must visit Carrefour Centrale on Wednesday and Friday
    constraint visits_r[2][6][0][5] ;
    constraint visits_r[4][4][0][5] ;
    // To simplify we say it only does one pickup (in reality he delivers a centre afterwards)
    constraint n_ramasses[2][6][0] == 1;
    constraint n_ramasses[4][4][0] == 1;

    // The PL cannot deliver Toulouse/Grande-Bretagne
    index_gde_bretagne = 24 ;
    for [c in 0...n : index_to_centre[c] == index_gde_bretagne ] {
        constraint and[d in 0...n_days]
            [trip in 0...n_trips](!visits_l[d][0][trip][c]);
    }

    // Bessières cannot be delivered with camions frigos
    index_bessieres = 3 ;
    for [c in 0...n : index_to_centre[c] == index_bessieres] {
        for [v in allowed_vehicles : can_carry[v]["F"]] {
            constraint and[d in 0...n_days]
                [trip in 0...n_trips](!visits_l[d][v][trip][c]);
        }
    }

    // St Orens has to be picked up with a specific camion frigo
    index_st_orens = 1;
    for [v in allowed_vehicles : v != 2] {
        constraint and[d in 0...n_days](!visits_r[d][v][0][index_st_orens]);
    }

    // "Livraisons de Ramasses" these deliver the food gathered in previous pickups. Only a single delivery can be done each time.
    for [dvtc in livraisons_de_ramasses] {
        d = dvtc[0] ;
        v = dvtc[1] ;
        trip = dvtc[2] ;
        c = dvtc[3] ;
        for [c2 in 0...n] {
            if (index_to_centre[c2] == c) {
                constraint visits_l[d][v][trip][c2] == 1;
            } else {
                constraint visits_l[d][v][trip][c2] == 0;
            }
        }
    }
}

function add_time_windows() {
    for [d in 0...n_days] {
        for [v in allowed_vehicles] {
            for [trip in 0...n_trips] {
                for [c in 0...n : !is_ldr[d][v][trip][c] && !delivery_allowed[d][c]] {
                    constraint !visits_l[d][v][trip][c] ;
                }
            }
            for [trip in 1...n_trips] {
                constraint count(ramasses[d][v][trip]) == 0 ;
            }
            for [p in 0...n_pdr] {
                if (!pickup_d_p[d][p] || !can_carry[v][pdrs[p]["product_type"]]) {
                    constraint !visits_r[d][v][trip][p] ;
                }
            }
        }
    }
}

function add_capacity_constraints() {

    for [d in 0...n_days] {
        for [v in allowed_vehicles] {
            for [trip in 0...n_trips] {
                // Delivery capacity & size constraints
                constraint sum[c in 0...n : delivery_allowed[d][c]](palettes[d][v][trip][c] + 0.5 * (demi_palettes[d][v][trip][c]["F"] + demi_palettes[d][v][trip][c]["S"])) <= vehicles[v]["size"] ;
                constraint sum[c in 0...n : delivery_allowed[d][c]](load[d][v][trip][c]) <= vehicles[v]["capacity"] ;

                // Pickup capacity & size constraints
                constraint sum[p in 0...n_pdr : pickup_d_p[d][p]](pdrs[p]["weight"] * visits_r[d][v][trip][p]) <= vehicles[v]["capacity"];
                // We assume 2 palettes are used for each pickup
                constraint sum[p in 0...n_pdr : pickup_d_p[d][p]](2 * visits_r[d][v][trip][p]) <= vehicles[v]["size"] ;
            }
        }
        constraint sum[c in 0...n : delivery_allowed[d][c]][v in allowed_vehicles]
            [trip in 0...n_trips](norvegiennes[d][v][trip][c]) <= params["n_norvegiennes"] ;
    }
}

function add_duration_constraints() {
    for [d in 0...n_days] {
        for [v in allowed_vehicles] {
            for [trip in 0...n_trips] {
                // Deliver first, then pickup
                route_dist[d][v][trip] <- (liv[d][v][trip] ?  dist_from_depot[livraisons[d][v][trip][0]] 
                            + sum(1...n_livraisons[d][v][trip], i => dist_matrix[livraisons[d][v][trip][i - 1]][livraisons[d][v][trip][i]])
                            + (ram[d][v][trip] ? 
                                dist_matrix[livraisons[d][v][trip][n_livraisons[d][v][trip]-1]][ramasses[d][v][trip][0] + n]
                                : dist_to_depot[livraisons[d][v][trip][n_livraisons[d][v][trip]-1]])

                            : (ram[d][v][trip] ? dist_from_depot[ramasses[d][v][trip][0] + n] : 0)) 
                        + (ram[d][v][trip] ? dist_to_depot[ramasses[d][v][trip][n_ramasses[d][v][trip]-1] + n]  
                            + sum(1...n_ramasses[d][v][trip], i => dist_matrix[ramasses[d][v][trip][i - 1] + n][ramasses[d][v][trip][i] + n]): 0)
                        ;

                // TODO : refactor this garbage

                trip_duration[d][v][trip] <- (liv[d][v][trip] ? 
                                                time_from_depot[livraisons[d][v][trip][0]] 
                                                    + sum(1...n_livraisons[d][v][trip], i => time_matrix[livraisons[d][v][trip][i - 1]][livraisons[d][v][trip][i]])
                                                    + params["wait_at_centres"] * 60 * n_livraisons[d][v][trip]
                                                    + (ram[d][v][trip] ? 
                                                        time_matrix[livraisons[d][v][trip][n_livraisons[d][v][trip]-1]][ramasses[d][v][trip][0] + n] 
                                                        + time_to_depot[ramasses[d][v][trip][n_ramasses[d][v][trip]-1] + n]
                                                        : time_to_depot[livraisons[d][v][trip][n_livraisons[d][v][trip]-1]])
                                                : (ram[d][v][trip] ? time_from_depot[ramasses[d][v][trip][0] + n] + time_to_depot[ramasses[d][v][trip][n_ramasses[d][v][trip]-1] + n] : 0))
                                                
                                            + (ram[d][v][trip] ? 
                                                sum(1...n_ramasses[d][v][trip], i => time_matrix[ramasses[d][v][trip][i - 1] + n][ramasses[d][v][trip][i] + n]) + params["wait_at_pdrs"] * 60 * n_ramasses[d][v][trip] : 0)
                                            ;
                constraint n_livraisons[d][v][trip] + n_ramasses[d][v][trip] <= params["max_stops"] ;
            }
            constraint ram[d][v][0] ? (trip_duration[d][v][0] <= params["max_tour_duration_with_pickup"] * 60) : true; 
            tour_duration[d][v] <- sum[trip in 0...n_trips](trip_duration[d][v][trip]);
            constraint tour_duration[d][v] <= params["max_tour_duration"] * 60;
        }
    }
}

function create_variables() {

    livraisons[d in 0...n_days][v in allowed_vehicles][trip in 0...n_trips] <- list(n);
    ramasses[d in 0...n_days][v in allowed_vehicles][trip in 0...n_trips] <- list(n_pdr);

    for [d in 0...n_days] {
        for [v in allowed_vehicles] {
            for [trip in 0...n_trips] {

                // Deliver quantity for each product type (ambiant/frais/surgelé)
                for [c in 0...n] {
                    for [t in product_types] {
                        max_demands = 0;
                        for [s in 0...n_scenarios] {
                            if (demands[s][c][t] > max_demands) {
                                max_demands = demands[s][c][t];
                            }
                        }

                        no_delivery = is_ldr[d][v][trip][c] || !delivery_allowed[d][c] || !can_carry[v][t];

                        quantity[d][v][trip][c][t] <-
                                no_delivery ? 
                                0 :int(0, max_demands) ;
                    }
                }

                norvegiennes[d][v][trip][c in 0...n] <- !is_ldr[d][v][trip][c] && delivery_allowed[d][c] && can_carry[v]["S"] ? int(0, params["n_norvegiennes"]) : 0;

                n_livraisons[d][v][trip] <- count(livraisons[d][v][trip]) ;
                n_ramasses[d][v][trip] <- count(ramasses[d][v][trip]) ;

                liv[d][v][trip] <- n_livraisons[d][v][trip] > 0 ;
                ram[d][v][trip] <- n_ramasses[d][v][trip] > 0 ;

                visits_l[d][v][trip][c in 0...n] <- contains(livraisons[d][v][trip], c) ;
                visits_r[d][v][trip][p in 0...n_pdr] <- contains(ramasses[d][v][trip], p);

                
                for [c in 0...n] {
                    load[d][v][trip][c] <- visits_l[d][v][trip][c] * sum[t in product_types](quantity[d][v][trip][c][t]);

                    palettes[d][v][trip][c] <- ceil(visits_l[d][v][trip][c] * quantity[d][v][trip][c]["A"] / params["max_palette_capacity"]) ;
                    demi_palettes[d][v][trip][c]["F"] <- 2 * ceil(visits_l[d][v][trip][c] * quantity[d][v][trip][c]["F"] / params["demi_palette_capacity"]);
                    demi_palettes[d][v][trip][c]["S"] <- delivery_allowed[d][c] && can_carry[v]["S"] && vehicles[v]["allows_isotherm_cover"] ?
                                            max(0, ceil(
                                                (visits_l[d][v][trip][c] * quantity[d][v][trip][c]["S"] - 
                                                norvegiennes[d][v][trip][c] * params["norvegienne_capacity"]) / params["demi_palette_capacity"]
                                                )) : 
                                            0;

                    constraint visits_l[d][v][trip][c] * quantity[d][v][trip][c]["S"] <= 
                                    norvegiennes[d][v][trip][c] * params["norvegienne_capacity"] +
                                    demi_palettes[d][v][trip][c]["S"] * params["demi_palette_capacity"] ;
                } 
            }
        }
    }
}

function add_cover_constraints() {
    // Every centre & pdr has to be delivered/picked up at least once
    constraint cover[d in 0...n_days][v in allowed_vehicles][trip in 0...n_trips](livraisons[d][v][trip]);
    constraint cover[d in 0...n_days][v in allowed_vehicles][trip in 0...n_trips](ramasses[d][v][trip]);

    // Pickup each site enough times each week
    for [p in 0...n_pdr] {
        for [d in 0...n_days] {
            if (pickup_d_p[d][p]) {
                constraint xor[v in allowed_vehicles][trip in 0...n_trips](visits_r[d][v][trip][p]) ;
                // constraint or[v in allowed_vehicles : can_carry[v][pdrs[p]["product_type"]]](visits_r[d][v][trip][p]) ;
            } else {
                // constraint and[v in allowed_vehicles](visits_r[d][v][trip][p] == 0) ;
            }
        }
    }

    // Redundant
    for [d in 0...n_days] {
        // constraint disjoint[v in allowed_vehicles](ramasses[d][v][trip]) ;
        for [p in 0...n_pdr : pickup_d_p[d][p]] {
            constraint sum[v in allowed_vehicles : can_carry[v][pdrs[p]["product_type"]]][trip in 0...n_trips](visits_r[d][v][trip][p]) == 1 ;
        }
    }

    // Don't visit the same centre more than 4 times
    for [c in 0...n] {
        constraint sum[d in 0...n_days : delivery_allowed[d][c]]
                    [v in allowed_vehicles : !is_ldr[d][v][trip][c]]
                    [trip in 0...n_trips]
            (visits_l[d][v][trip][c]) <= 4 ;
        
        // // redundant
        // if (sum[t in product_types](demands[c][t]) > 1200) {
        //     constraint sum[d in 0...n_days : delivery_allowed[d][c]][v in allowed_vehicles : !is_ldr[d][v][trip][c]]
        //         (visits_l[d][v][trip][c]) >= 2 ;
        // }
    }

}

function add_demand_constraints() {
    // Meet demands for each product type in each scenario
    // demands_met[s in 0...n_scenarios] <- bool();
    for [s in 0...n_scenarios] {
        for [c in 0...n] {
            for [t in product_types] {
                if (n_scenarios == 1) {
                    constraint sum[d in 0...n_days : delivery_allowed[d][c]]
                            [v in allowed_vehicles : !is_ldr[d][v][trip][c]]
                            [trip in 0...n_trips]
                            (visits_l[d][v][trip][c] * quantity[d][v][trip][c][t]) == demands[s][c][t];
                } else {
                    constraint sum[d in 0...n_days : delivery_allowed[d][c]]
                            [v in allowed_vehicles : !is_ldr[d][v][trip][c]]
                            [trip in 0...n_trips]
                            (visits_l[d][v][trip][c] * quantity[d][v][trip][c][t]) >= demands[s][c][t];
                }
            }
        }
    }
}

function add_objectives() {
    used[v in allowed_vehicles] <- or[d in 0...n_days][trip in 0...n_trips](liv[d][v][trip] || ram[d][v][trip]) ;
    fixed_costs <- sum[v in allowed_vehicles](vehicles[v]["fixed_cost"] * used[v]);
    n_used <- sum[v in allowed_vehicles](used[v]) ;

    variable_costs <- sum[d in 0...n_days][v in allowed_vehicles][trip in 0...n_trips](vehicles[v]["cost_per_km"] * route_dist[d][v][trip] / 1000);
    total_costs <- variable_costs + fixed_costs ;
    total_distance <- sum[d in 0...n_days][v in allowed_vehicles][trip in 0...n_trips](route_dist[d][v][trip]);

    // constraint total_costs <= 900 ;
    minimize total_costs ;
}

function model() {
    create_variables();
    add_time_windows();
    add_capacity_constraints();
    add_duration_constraints();
    add_cover_constraints();
    add_specific_requirements() ;
    add_demand_constraints();
    add_objectives();
}


function set_initial_solution(force, specific_vd) {
    /* Set an initial solution.
    If force == true, add constraints instead to match the initial solution instead.
    Only fix tours of vehicle-days in specific_vd */

    if (initfile == nil) return;

    if (specific_vd == nil) {
        fix_vd = {} ;
        for [vd in 0...n_days * m]
            fix_vd.add(vd);
    }
    else 
        fix_vd = specific_vd ;

    tours_init = json.parse(initfile)["tours"];

    for [d in 0...n_days] {
        for [v in allowed_vehicles] {

            vd = v + d * m;

            if (fix_vd[vd] == nil)
                continue;

            key = d + ", " + v;

            if (tours_init[key] == nil) {
                if (force) {
                    constraint n_livraisons[d][v] == 0 ;
                    constraint n_ramasses[d][v] == 0 ;
                }
                continue;
            }

            livraison_init = {};
            ramasse_init = {};
            if (!force) {
                livraisons[d][v].value.clear();
                ramasses[d][v].value.clear();
            }
            i = 0;
            j = 0 ;
            for [place in tours_init[key]] {
                if (place["type"] == "Livraison") {
                    index = -1 ;
                    for [c in 0...n] {
                        if (index_to_centre[c] == place["index"]) {
                            index = c ;
                            break ;
                        }
                    }
                    if (index == -1) {
                        println("Warning : skipping ", place["name"], " ", place["index"], " in input solution") ;
                        break ;
                    }

                    livraison_init.add(index);
                    if (force) {
                        constraint livraisons[d][v][i] == index;
                        i += 1 ;
                        constraint visits_l[d][v][index] == 1 ;
                        constraint norvegiennes[d][v][index] == place["norvegiennes"];
                        constraint quantity[d][v][index]["A"] == place["delivery"][0] ;
                        constraint quantity[d][v][index]["F"] == place["delivery"][1] ;
                        constraint quantity[d][v][index]["S"] == place["delivery"][2] ;
                    } else {
                        livraisons[d][v].value.add(index);

                        try {
                            if (can_carry[v]["A"]) {
                                quantity[d][v][index]["A"].value = place["delivery"][0];
                            }
                            if (can_carry[v]["F"]) {
                                quantity[d][v][index]["F"].value = place["delivery"][1];
                            }
                            if (can_carry[v]["S"]) {
                                quantity[d][v][index]["S"].value = place["delivery"][2];
                                norvegiennes[d][v][index].value = place["norvegiennes"];
                            }
                        } catch (error) {
                        }
                    }

                } else {
                    ramasse_init.add(place["index"] - n_centres);
                    if (force){
                        constraint ramasses[d][v][j] == place["index"] - n_centres;
                        j += 1;
                    } else {
                        ramasses[d][v].value.add(place["index"] - n_centres);
                    }
                }
            }
            if (force) {
                constraint count(livraisons[d][v]) == livraison_init.count();
                constraint count(ramasses[d][v]) == ramasse_init.count();
            }
        }
    }
    
}

function fix(percentage) {
    // Fix some percentage of the tours to locally improve the solution (basic LNS procedure)
    vd_to_fix = {};
    n_fix = round(percentage * n_days * m);

    vd = 0 ;
    randgen = random.create();
    while (vd_to_fix.count() < n_fix) {
        vd_to_fix.add(randgen.next(0, n_days * m)) ;
    }

    set_initial_solution(true, vd_to_fix);
}

function set_current_tours() {
    /* Constrains the solution to follow the current tours */
    weeknum = params["week"] == "ODD" ? "1" : "2" ;
    current_tours = json.parse("data/current/tours_tournees_actuelles_w" + weeknum + ".json");
    for [d in 0...n_days] {
        for [v in allowed_vehicles] {
            key = "(" + d + ", " + v + ")" ;
            if (current_tours[key] == nil) {
                constraint n_livraisons[d][v] == 0 ;
                constraint n_ramasses[d][v] == 0 ;
                continue;
            }
            init_liv[d][v] = {} ;
            init_ram = {} ;
            i = 0 ;
            j = 0 ;
            for [k in 1...current_tours[key].count()-1] {
                for [c in 0...n] {
                    if (index_to_centre[c] == current_tours[key][k]) {
                        constraint livraisons[d][v][i] == c;
                        i += 1 ;
                        // init_liv.add(c);
                        break;
                    }
                }
                for [p in 0...n_pdr] {
                    if (p + n_centres == current_tours[key][k]) {
                        constraint ramasses[d][v][j] == p;
                        j += 1 ;
                        // init_ram.add(p);
                        break;
                    }
                }
            }
            constraint n_livraisons[d][v] == i;
            constraint n_ramasses[d][v] == j;
        }
    }
}

function check_solution_robustness(n_runs, sigma) {
    /*
    Checks the robustness of the tours provided as input

    Simulates many demands from the given distribution, then for each one
    constrains the solution to be almost identical (we can only add products to existing tours, no new stops or tours)
    Prints how many feasible solutions it finds
    */

    // TODO : make sure this still works if needed

    if (initfile == nil) {
        println("Please enter an input file to check") ;
        return ;
    }

    // this is hardcoded for 1.15 robustness factor
    for [p in product_types]{
        for [c in 0...n] {
            original_demands[c][p] = demands[c][p] * (1/1.15) ;
        }
    }
    rng = random.create() ;

    feasible = 0 ;
    st = datetime.now() ;
    for [run in 0...n_runs] {
        println("Simulating run ", run, "/",n_runs,"... time=",datetime.now()-st,"s feasible=",feasible,"/",run) ;
        for [p in product_types] {
            for [c in 0...n] {
                rand = rng.nextNormal(1, sigma) ;
                rand = max(1 - 1.5*sigma, rand) ;
                rand = min(1 + 1.5*sigma, rand) ;
                demands[c][p] = max(1, ceil(original_demands[c][p] * rand)) ;
            }
        }

        with(ls = localsolver.create()) {
            // ls.addCallback("ITERATION_TICKED", callback);

            model();
            set_initial_solution(true, nil); // Force the tours to be identical to the given solution
            ls.model.close();
            set_initial_solution(false, nil); // Give the given solution as a hint
            ls.param.timeLimit = 2 ;
            ls.param.verbosity = 0 ;
            ls.solve();

            if (ls.solution.status == "FEASIBLE") {
                feasible += 1 ;
            }
        }
    }

    perc_feasible = round(100 * feasible / n_runs) ;
    println("Percentage of feasible runs : ", perc_feasible, "%") ;
}

function aggregate_tours() {
    for [d in 0...n_days] {
        for [v in allowed_vehicles] {
            tours[d][v] = {} ;
            for [trip in 0...n_trips] {
                for [node in livraisons[d][v][trip].value] {
                    tours[d][v].add(node) ;
                }
                for [node in ramasses[d][v][trip].value] {
                    tours[d][v].add(node+n) ;
                }
                if (trip < n_trips - 1 && (liv[d][v][trip+1].value || ram[d][v][trip+1].value)) {
                    tours[d][v].add(-1) ;
                }
            }
        }
    }
    return tours;
}

function write_solution() {
    if (outfile == nil) return;
    local outf = io.openWrite(outfile);

    sol = {};
    sol["week"] = params["week"] ;
    sol["total_distance"] = total_distance.value;
    sol["total_costs"] = total_costs.value;
    sol["variable_costs"] = variable_costs.value;
    sol["fixed_costs"] = fixed_costs.value;
    sol["vehicles_used"] = {} ;
    for [v in 0...m] {
        if (vehicles[v]["allowed"]) {
            sol["vehicles_used"].add(used[v].value == 1 ? true : false) ;
        } else {
            sol["vehicles_used"].add(false) ;
        }
    }
    sol["tour_durations"] = {} ;

    sol["tours"] = {};

    tours = aggregate_tours();

    for [d in 0...n_days] {
        for [v in allowed_vehicles] {
            if (tours[d][v].count() == 0) continue;

            key = d + ", " + v ;

            sol["tour_durations"][key] = tour_duration[d][v].value ;
            sol["tours"][key] = {};

            trip = 0;
            for [node in tours[d][v]] {
                place = {};

                if (node == -1) {
                    place["index"] = 0 ;
                    place["name"] = centres[0]["name"];
                    place["type"] = "Livraison";
                    place["delivery"] = {0,0,0};
                    place["palettes"] = {0,0,0};
                    place["norvegiennes"] = 0;
                    trip += 1;
                } else if (node < n) {
                    place["index"] = index_to_centre[node];
                    place["name"] = centres[index_to_centre[node]]["name"];
                    place["type"] = "Livraison";
                    place["delivery"] = {
                        quantity[d][v][trip][node]["A"].value, 
                        quantity[d][v][trip][node]["F"].value, 
                        quantity[d][v][trip][node]["S"].value
                    };

                    place["palettes"] = {
                        palettes[d][v][trip][node].value,
                        demi_palettes[d][v][trip][node]["F"].value/2,
                        demi_palettes[d][v][trip][node]["S"].value/2
                    };

                    // ignore any superfluous norvégienne
                    place["norvegiennes"] = quantity[d][v][trip][node]["S"].value > demi_palettes[d][v][trip][node]["S"].value * params["demi_palette_capacity"] ? 
                                                norvegiennes[d][v][trip][node].value : 0;
                } else {
                    place["index"] = node - n + n_centres;
                    place["name"] = pdrs[node - n]["name"];
                    place["type"] = "Ramasse";
                }

                sol["tours"][key].add(place);
            }
        }
    }
    json.dump(sol, outf);
}

function output(){
    write_solution();
}

function callback(ls, cbType) {
    local stats = ls.statistics;
    local obj <- ls.model.objectives[0];
    if (ls.solution.status == "FEASIBLE" && obj.value < lastBestValue) {
        lastBestRunningTime = stats.runningTime;
        lastBestValue = obj.value;
        lastSolutionWritten = false ;

        perc_fixed = fixed_costs.value/total_costs.value ;
        perc_var = variable_costs.value/total_costs.value ;

        println(
            "[   ", stats.runningTime, "s] : ", total_costs.value, "E   (", 
            round(perc_fixed*100), "%F ", round(perc_var*100), "%V)   ", 
            total_distance.value/1000, "km   ", n_used.value, " vehicles"
        );
    }
    if (outfile != nil && stats.runningTime - lastBestRunningTime > 5 && !lastSolutionWritten) {
        println(">>>>>>> No improvement during 5 seconds: writing current solution to file");
        write_solution();
        lastSolutionWritten = true ;
    }
}


function main(args) {

    usage = "Usage: ls_solver problem_file (initsol | 'nil') (outfile | 'nil') [timelimit]";
    if (args.count() < 1) {
        throw usage;
    }

    problem_file = args[0];
    initfile = args[1] == "nil" || args.count() < 2 ? nil : args[1];
    outfile = args[2] == "nil"  || args.count() < 3 ? nil : args[2];
    timelimit = args.count() < 4 ? nil : args[3].toInt();

    input();

    lastBestValue = inf - 1 ;
    lastBestRunningTime = 0;
    lastSolutionWritten = false;

    with(ls = localsolver.create()) {
        ls.addCallback("ITERATION_TICKED", callback);
        ls.param.verbosity = 0 ;
        // ls.param.verbosity = 2 ;

        model();

        // set_current_tours();

        // // Force the solution to match the input solution
        // set_initial_solution(true, nil);

        ls.model.close();

        if (timelimit != nil) {
            ls.param.timeLimit = timelimit ;
        }

        if (initfile != nil) {
            set_initial_solution(false, nil);
        }

        ls.solve();
        output();
    }
}