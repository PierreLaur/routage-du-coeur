use csv;
use json;
use io;
use localsolver;
use random;
use datetime;

function make_index_map(centres, n_centres, week) {
    // list variables require 0-based dense indices 
    // drop the depot & any centres that are not delivered this week to achieve that
    // create a map to be able to link new indices to indices in the problem data
    index_to_centre = {} ;
    i = 0 ;
    for [c in 1...n_centres] {
        semaine_livraison = centres[c]["delivery_week"] ;
        if(semaine_livraison == "ANY" || semaine_livraison == "ODD" && week == 1 || semaine_livraison == "EVEN" && week == 2) {
            index_to_centre[i] = c ;
            i += 1 ;
        }
    }

    return index_to_centre;
}

function read_demands(demand_dict, n, index_to_centre) {
    // where is this function :(
    int_to_string = {
        0 : "0",
        1 : "1",
        2 : "2",
        3 : "3",
        4 : "4",
        5 : "5",
        6 : "6",
        7 : "7",
        8 : "8",
        9 : "9",
        10 : "10",
        11 : "11",
        12 : "12",
        13 : "13",
        14 : "14",
        15 : "15",
        16 : "16",
        17 : "17",
        18 : "18",
        19 : "19",
        20 : "20",
        21 : "21",
        22 : "22",
        23 : "23",
        24 : "24",
        25 : "25",
        26 : "26",
        27 : "27",
        28 : "28",
        29 : "29",
        30 : "30"
    };

    // Demands
    for [c in 0...n] {
        centre = index_to_centre[c] ;
        demands[c] = demand_dict[int_to_string[centre]];
    }

    return demands;
}

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

    index_to_centre = make_index_map(centres, n_centres, week);

    // The number of centres to serve this week
    n = index_to_centre.count();

    m = vehicles.count();
    n_pdr = pdrs.count();
    n_nodes = n + n_pdr;
    n_days = 5;
    n_trips = params["max_trips"];

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

    demands = read_demands(problem["demands"], n, index_to_centre);

    product_types = {"A", "F", "S"};

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

    // Is d, v, trip, c supposed to be a livraison de ramasse ?
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



function create_variables() {

    livraisons[d in 0...n_days][v in allowed_vehicles][trip in 0...n_trips] <- list(n);
    ramasses[d in 0...n_days][v in allowed_vehicles][trip in 0...n_trips] <- list(n_pdr);

    for [d in 0...n_days] {
        for [v in allowed_vehicles] {
            for [trip in 0...n_trips] {
                // whether each demand is delivered by this day-vehicle-trip
                for [c in 0...n] {
                    for [i in 0...demands[c].count()] {
                        demand = demands[c][i];

                        can_deliver = true;
                        can_deliver = can_deliver && can_carry[v][demand["product_type"]];
                        if (
                            demand["product_type"] == "S"
                            && demand["palettes"] > 0
                            && !vehicles[v]["allows_isotherm_cover"]
                        ) {
                            can_deliver = false;
                        }

                        if (!delivery_allowed[d][c]) {
                            if (!is_ldr[d][v][trip][c]) {
                                can_deliver = false;
                            }
                        }

                        if (is_ldr[d][v][trip][c]) {
                            if (!delivery_allowed[d][c] || demand["product_type"] != "S") {
                                can_deliver = false;
                            }
                        }

                        if (can_deliver) {
                            delivers[d][v][trip][c][i] <- bool();
                        } else {
                            delivers[d][v][trip][c][i] <- false;
                        }

                    }
                }

                n_livraisons[d][v][trip] <- count(livraisons[d][v][trip]) ;
                n_ramasses[d][v][trip] <- count(ramasses[d][v][trip]) ;

                liv[d][v][trip] <- n_livraisons[d][v][trip] > 0 ;
                ram[d][v][trip] <- n_ramasses[d][v][trip] > 0 ;

                visits_l[d][v][trip][c in 0...n] <- contains(livraisons[d][v][trip], c) ;
                visits_r[d][v][trip][p in 0...n_pdr] <- contains(ramasses[d][v][trip], p);
                
                for [c in 0...n] {
                    if (!is_ldr[d][v][trip][c]) {
                        // constraint and[i in 0...demands[c].count()](!delivers[d][v][trip][c][i]) || visits_l[d][v][trip][c];
                        // for [i in 0...demands[c].count()] {
                        //     constraint delivers[d][v][trip][c][i] <= visits_l[d][v][trip][c];
                        // }
                    }

                    load[d][v][trip][c] <- !visits_l[d][v][trip][c] ? 0 : sum[i in 0...demands[c].count()]
                                                (delivers[d][v][trip][c][i] * demands[c][i]["weight"]);
                    palettes[d][v][trip][c] <-  visits_l[d][v][trip][c] * sum[i in 0...demands[c].count()]
                                                (delivers[d][v][trip][c][i] * demands[c][i]["palettes"]);
                    norvegiennes[d][v][trip][c] <- visits_l[d][v][trip][c] * sum[i in 0...demands[c].count()]
                                                (delivers[d][v][trip][c][i] * demands[c][i]["norvegiennes"]);
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

    // both PLs cannot deliver Escalquens
    index_escalquens = 8;
    for [c in 0...n : index_to_centre[c] == index_escalquens ] {
        for [v in {0, 1}] {
            constraint and[d in 0...n_days]
                [trip in 0...n_trips](!visits_l[d][v][trip][c]);
        }
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

    // Fenouillet not in first (too early)
    index_fenouillet = 9;
    for [c in 0...n : index_to_centre[c] == index_fenouillet] {
        // constraint and[d in 0...n_days][v in 0...m](indexOf(livraisons[d][v][0], c) != 0 || n_livraisons[d][v][0] <= 2);
        for [d in 0...n_days][v in 0...m] {
            constraint and[trip in 0...n_trips](!visits_l[d][v][trip][c]) 
                        || tour_duration[d][v] <= (params["max_tour_duration"] - 45) * 60;
        }
    }

    // Fronton in first only (they distribute food in the morning)
    index_fronton = 11;
    for [c in 0...n : index_to_centre[c] == index_fronton] {
        constraint and[d in 0...n_days][v in 0...m](indexOf(livraisons[d][v][0], c) <= 0);
        constraint and[d in 0...n_days][v in 0...m][trip in 1...n_trips](!visits_l[d][v][trip][c]);
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

    // Leclerc Rouffiac with camions frigos except on Tuesdays
    index_lc_rouffiac = 4;
    for [d in 0...n_days : d != 1] {
        constraint and[v in allowed_vehicles : !can_carry[v]["F"]]
            (!visits_r[d][v][0][index_lc_rouffiac]);
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
                constraint !ram[d][v][trip] ;
                constraint (liv[d][v][trip] + ram[d][v][trip] == 0) || (liv[d][v][trip-1] + ram[d][v][trip-1] >= 1) ;
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
                constraint sum[c in 0...n : delivery_allowed[d][c]](palettes[d][v][trip][c]) <= vehicles[v]["size"] ;
                constraint sum[c in 0...n : delivery_allowed[d][c]](load[d][v][trip][c]) <= vehicles[v]["capacity"] ;

                // Pickup capacity & size constraints
                constraint sum[p in 0...n_pdr : pickup_d_p[d][p]](pdrs[p]["weight"] * visits_r[d][v][trip][p]) <= vehicles[v]["capacity"];
                
                // We assume 2 palettes are used for each pickup
                constraint sum[p in 0...n_pdr : pickup_d_p[d][p]](pdrs[p]["palettes"] * visits_r[d][v][trip][p]) <= vehicles[v]["size"] ;
            }
        }
        constraint sum[c in 0...n : delivery_allowed[d][c]][v in allowed_vehicles]
            [trip in 0...n_trips](norvegiennes[d][v][trip][c]) <= params["n_norvegiennes"] ;
    }
}

function add_duration_constraints() {
    index_luchon = 2;
    index_montrejeau = 15;
    for [c in 0...n] {
        if (index_to_centre[c] == index_luchon) {
            node_luchon = c;
        }
        if (index_to_centre[c] == index_montrejeau) {
            node_montrejeau = c;
        }
    }
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

                liv_time <- liv[d][v][trip] ? 
                                time_from_depot[livraisons[d][v][trip][0]] 
                                + sum(1...n_livraisons[d][v][trip], i => time_matrix[livraisons[d][v][trip][i - 1]][livraisons[d][v][trip][i]])
                                + params["wait_at_centres"] * 60 * n_livraisons[d][v][trip]
                                : 0;
                ram_time <- ram[d][v][trip] ?
                                sum(1...n_ramasses[d][v][trip], i => time_matrix[ramasses[d][v][trip][i - 1] + n][ramasses[d][v][trip][i] + n])
                                + params["wait_at_pdrs"] * 60 * n_ramasses[d][v][trip] 
                                + time_to_depot[ramasses[d][v][trip][n_ramasses[d][v][trip]-1] + n]
                                : 0;
                
                other_time <- ram[d][v][trip] && liv[d][v][trip] ? 
                                    time_matrix[livraisons[d][v][trip][n_livraisons[d][v][trip]-1]][ramasses[d][v][trip][0] + n] 
                                    : 0
                            + ram[d][v][trip] && !liv[d][v][trip] ? time_from_depot[ramasses[d][v][trip][0] + n] 
                                    : 0
                            + liv[d][v][trip] && !ram[d][v][trip] ? time_to_depot[livraisons[d][v][trip][n_livraisons[d][v][trip]-1]] 
                                    : 0;

                trip_duration[d][v][trip] <- liv_time + ram_time + other_time;
            }
            constraint sum[trip in 0...n_trips](n_livraisons[d][v][trip] + n_ramasses[d][v][trip]) <= params["max_stops"] ;
            constraint ram[d][v][0] ? (trip_duration[d][v][0] <= params["max_tour_duration_with_pickup"] * 60) : true; 
            tour_duration[d][v] <- sum[trip in 0...n_trips](trip_duration[d][v][trip])
                                + sum[trip in 1...n_trips](liv[d][v][trip]) * params["wait_between_trips"] * 60;

            leeway <- 90 * 60 * or[trip in 0...n_trips]
                        ((node_luchon != nil) ? visits_l[d][v][trip][node_luchon] : false || 
                        (node_montrejeau != nil) ? visits_l[d][v][trip][node_montrejeau] : false);

            constraint tour_duration[d][v] <= params["max_tour_duration"] * 60 + leeway;
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
                constraint xor[v in allowed_vehicles : can_carry[v][pdrs[p]["product_type"]]][trip in 0...n_trips](visits_r[d][v][trip][p] == 1) ;
                // constraint or[v in allowed_vehicles : can_carry[v][pdrs[p]["product_type"]]][trip in 0...n_trips](visits_r[d][v][trip][p]) ;
            } else {
                // constraint and[v in allowed_vehicles][trip in 0...n_trips](visits_r[d][v][trip][p] == 0) ;
            }
        }
    }

    // Redundant
    for [d in 0...n_days] {
        constraint disjoint[v in allowed_vehicles][trip in 0...n_trips](ramasses[d][v][trip]) ;
        // for [p in 0...n_pdr : pickup_d_p[d][p]] {
        //     constraint sum[v in allowed_vehicles : can_carry[v][pdrs[p]["product_type"]]][trip in 0...n_trips](visits_r[d][v][trip][p]) == 1 ;
        // }
    }

    fourgon_capacity = 1200 ;
    for [c in 0...n] {

        sum_demands = sum[dem in demands[c]](dem["weight"]);
        if (sum_demands > fourgon_capacity) {
            min_visits = 2;
        } else {
            min_visits = 1;
        }

        // more than 3 visits is unreasonable in practice
        max_visits = 3;

        n_visits <- sum[d in 0...n_days : delivery_allowed[d][c]]
                    [v in allowed_vehicles : !is_ldr[d][v][trip][c]]
                    [trip in 0...n_trips]
            (visits_l[d][v][trip][c]) ;
        // constraint n_visits <= max_visits;
        // constraint n_visits >= min_visits;
    }
}

function add_demand_constraints() {
    for [c in 0...n] {
        for [i in 0...demands[c].count()] {
            // constraint xor[d in 0...n_days : delivery_allowed[d][c]][v in allowed_vehicles][trip in 0...n_trips]
            //         (visits_l[d][v][trip][c] && delivers[d][v][trip][c][i]);
            constraint sum[d in 0...n_days : delivery_allowed[d][c]][v in allowed_vehicles][trip in 0...n_trips]
                (visits_l[d][v][trip][c] * delivers[d][v][trip][c][i]) == 1;
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
    // constraint n_used <= 9 ;
    minimize total_costs ;

    // add_robustness_cuts();
}

function add_robustness_cuts() {
    constraint n_used >= 9;
    constraint used[1];
}

function model() {
    create_variables();
    add_time_windows();
    add_capacity_constraints();
    add_duration_constraints();
    add_cover_constraints();
    add_specific_requirements() ;
    add_demand_constraints();
}


function set_initial_solution(solution_file, force, days) {
    /* 
    Set an initial solution.
    If force == true, add constraints instead to match the initial solution instead.
    */

    fix_d = {};
    if (days == nil) {
        for [d in 0...n_days]
            fix_d[d] = true;
    } else {
        for [d in days] {
            fix_d[d] = true;
        }
    }

    tours_init = json.parse(solution_file)["tours"];

    for [d in 0...n_days] {
        if (fix_d[d] == nil)
            continue;

        for [v in allowed_vehicles] {

            key = d + ", " + v;

            if (!force) {
                livraisons[d][v][trip].value.clear();
                ramasses[d][v][trip].value.clear();
            }

            if (tours_init[key] == nil) {
                for [trip in 0...n_trips] {
                    if (force) {
                        constraint n_livraisons[d][v][trip] == 0 ;
                        constraint n_ramasses[d][v][trip] == 0 ;
                    } else {
                        livraisons[d][v][trip].value.clear();
                        ramasses[d][v][trip].value.clear();
                    }
                }
                continue;
            }
            livraison_init = {};
            ramasse_init = {};
            trip = 0;
            i = 0;
            j = 0 ;
            for [place in tours_init[key]] {
                if (place["index"] == 0) {
                    trip += 1 ;
                    i = 0 ;
                    livraison_init = {};
                    ramasse_init = {};
                    continue;
                }

                if (place["stop_type"] == "Livraison" || place["stop_type"] == "Liv_Ramasse") {
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
                        constraint livraisons[d][v][trip][i] == index;
                        i += 1 ;
                        constraint visits_l[d][v][trip][index] == 1 ;
                        constraint norvegiennes[d][v][trip][index] == place["norvegiennes"];
                        constraint palettes[d][v][trip][index] == sum[p in 0...3](place["palettes"][p]);
                    } else {
                        livraisons[d][v][trip].value.add(index);
                    }

                } else {
                    ramasse_init.add(place["index"] - n_centres);
                    if (force){
                        constraint ramasses[d][v][trip][j] == place["index"] - n_centres;
                        j += 1;
                    } else {
                        ramasses[d][v][trip].value.add(place["index"] - n_centres);
                    }
                }
            }
            if (force) {
                constraint count(livraisons[d][v][trip]) == livraison_init.count();
                constraint count(ramasses[d][v][trip]) == ramasse_init.count();
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

    
    if (initfile == nil) return;

    set_initial_solution(initfile, true, vd_to_fix);
}

function set_current_tours() {
    /* Constrains the solution to follow the current tours */
    weeknum = week == 1 ? "1" : "2" ;
    current_tours_file = "data/current/2602.json";
    current_tours = json.parse(current_tours_file);
    for [d in 0...n_days] {
        for [v in allowed_vehicles] {
            key = "(" + d + ", " + v + ")" ;
            if (current_tours[key] == nil) {
                for [trip in 0...n_trips] {
                    constraint n_livraisons[d][v][trip] == 0 ;
                    constraint n_ramasses[d][v][trip] == 0 ;
                }
                continue;
            }
            i = 0 ;
            j = 0 ;
            trip = 0 ;
            for [k in 1...current_tours[key].count()-1] {
                if (current_tours[key][k] == 0) {
                    constraint n_livraisons[d][v][trip] == i;
                    constraint n_ramasses[d][v][trip] == j;
                    trip += 1;
                    i = 0;
                    j = 0;
                    continue;
                }

                for [c in 0...n] {
                    if (index_to_centre[c] == current_tours[key][k]) {
                        constraint livraisons[d][v][trip][i] == c;
                        i += 1 ;
                        break;
                    }
                }
                for [p in 0...n_pdr] {
                    if (p + n_centres == current_tours[key][k]) {
                        constraint ramasses[d][v][trip][j] == p;
                        j += 1 ;
                        break;
                    }
                }
            }

            if (trip < n_trips) {
                constraint n_livraisons[d][v][trip] == i;
                constraint n_ramasses[d][v][trip] == j;
            }

            for [later_trips in trip+1...n_trips] {
                constraint n_livraisons[d][v][later_trips] == 0;
                constraint n_ramasses[d][v][later_trips] == 0;
            }
        }
    }
}

function set_same_visits(solution_file) {

    solution = json.parse(solution_file);
    tours_init = solution["tours"];

    recourses = {};
    for [d in 0...n_days] {
        for [v in allowed_vehicles] {

            key = d + ", " + v;

            if (tours_init[key] == nil) {
                for [trip in 0...n_trips] {
                    constraint n_ramasses[d][v][trip] == 0;
                }
                recourses.add(sum[trip in 0...n_trips](n_livraisons[d][v][trip]));
                continue;
            }

            livraison_init = {};
            ramasse_init = {};
            trip = 0;
            i = 0;
            j = 0 ;
            for [place in tours_init[key]] {
                if (place["index"] == 0) {
                    liv_diff <- n_livraisons[d][v][trip] - livraison_init.count();
                    recourses.add(liv_diff > 0 ? liv_diff : 0);
                    trip += 1 ;
                    i = 0 ;
                    livraison_init = {};
                    ramasse_init = {};
                    continue;
                }

                if (place["stop_type"] == "Livraison" || place["stop_type"] == "Liv_Ramasse") {
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
                    constraint or[trip in 0...n_trips](visits_l[d][v][trip][index]);
                    i += 1 ;

                } else {
                    ramasse_init.add(place["index"] - n_centres);
                    constraint ramasses[d][v][trip][j] == place["index"] - n_centres;
                    constraint visits_r[d][v][trip][place["index"] - n_centres] ;
                    j += 1;
                }
            }
            liv_diff <- n_livraisons[d][v][trip] - livraison_init.count();
            recourses.add(liv_diff > 0 ? liv_diff : 0);
            constraint n_ramasses[d][v][trip] == ramasse_init.count();

            for [later_trip in trip+1...n_trips] {
                recourses.add(n_livraisons[d][v][later_trip]);
                constraint n_ramasses[d][v][later_trip] == 0;
            }
        }
    }

    return recourses;
    
}



function aggregate_tours() {
    tours = {};
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

function write_solution(outfile) {
    if (outfile == nil) return;
    local outf = io.openWrite(outfile);

    sol = {};
    sol["week"] = week ;
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
                    place["stop_type"] = "Livraison";
                    place["delivery"] = {0,0,0};
                    place["palettes"] = {0,0,0};
                    place["norvegiennes"] = 0;
                    trip += 1;
                } else if (node < n) {
                    place["index"] = index_to_centre[node];
                    place["name"] = centres[index_to_centre[node]]["name"];
                    if (is_ldr[d][v][trip][node]) {
                        place["stop_type"] = "Liv_Ramasse";
                    } else {
                        place["stop_type"] = "Livraison";
                    }
                    place["delivery"] = {
                        sum[i in 0...demands[node].count() : demands[node][i]["product_type"] == "A"](demands[node][i]["weight"] * visits_l[d][v][trip][node].value * delivers[d][v][trip][node][i].value),
                        sum[i in 0...demands[node].count() : demands[node][i]["product_type"] == "F"](demands[node][i]["weight"] * visits_l[d][v][trip][node].value * delivers[d][v][trip][node][i].value),
                        sum[i in 0...demands[node].count() : demands[node][i]["product_type"] == "S"](demands[node][i]["weight"] * visits_l[d][v][trip][node].value * delivers[d][v][trip][node][i].value)
                    };

                    place["palettes"] = {
                        sum[i in 0...demands[node].count() : demands[node][i]["product_type"] == "A"](demands[node][i]["palettes"] * visits_l[d][v][trip][node].value * delivers[d][v][trip][node][i].value),
                        sum[i in 0...demands[node].count() : demands[node][i]["product_type"] == "F"](demands[node][i]["palettes"] * visits_l[d][v][trip][node].value * delivers[d][v][trip][node][i].value),
                        sum[i in 0...demands[node].count() : demands[node][i]["product_type"] == "S"](demands[node][i]["palettes"] * visits_l[d][v][trip][node].value * delivers[d][v][trip][node][i].value)
                    };

                    place["norvegiennes"] = sum[i in 0...demands[node].count()](demands[node][i]["norvegiennes"] * visits_l[d][v][trip][node].value * delivers[d][v][trip][node][i].value);
                } else {
                    place["index"] = node - n + n_centres;
                    place["name"] = pdrs[node - n]["name"];
                    place["stop_type"] = "Ramasse";
                }

                sol["tours"][key].add(place);
            }
        }
    }
    json.dump(sol, outf);
}

function output(){
    write_solution(outfile);
}



function callback(ls, cbType) {
    local stats = ls.statistics;
    local obj <- ls.model.objectives[0];
    if (ls.solution.status == "FEASIBLE" && obj.value != lastObjectiveValue) {

        if (obj.value < lastBestValue) {
            lastBestRunningTime = stats.runningTime;
            lastImprovementRunningTime = stats.runningTime;
            lastBestValue = obj.value;
            lastSolutionWritten = false ;

            perc_fixed = fixed_costs.value/total_costs.value ;
            perc_var = variable_costs.value/total_costs.value ;

            println(
                "[   ", stats.runningTime, "s] : ", total_costs.value, "E   (", 
                round(perc_fixed*100), "%F ", round(perc_var*100), "%V)   ", 
                total_distance.value/1000, "km   ", n_used.value, " vehicles"
            );
        } else {
            println(
                "[   ", stats.runningTime, "s] : ", lastBestValue, "E   ",
                "current = ", obj.value
            );
            lastImprovementRunningTime = stats.runningTime;
        }
        lastObjectiveValue = obj.value;
    }
    if (outfile != nil && stats.runningTime - lastBestRunningTime > 1 && !lastSolutionWritten) {
        write_solution(outfile);
        println(">>>>>>> wrote solution to ", outfile);
        lastSolutionWritten = true ;
    }

    if (improvement_limit != nil && stats.runningTime - lastImprovementRunningTime > improvement_limit) {
        ls.stop();
    }
}

function solve_multi(num_tries, time_limit) {
    improvement_limit = time_limit ;
    for [t in 0...num_tries] {
        with(ls = localsolver.create()) {
            ls.addCallback("ITERATION_TICKED", callback);
            ls.param.verbosity = 0 ;
            ls.param.seed = random.create().next(0, 10000);
            model();
            add_objectives();
            ls.model.close();
            ls.solve();
            ls.model.open();
        }
    }

    output();
}

function solve() {
    with(ls = localsolver.create()) {
        ls.addCallback("ITERATION_TICKED", callback);
        ls.param.verbosity = 2 ;
        ls.param.seed = random.create().next(0, 10000);

        model();
        add_objectives();
        
        // set_current_tours();

        // Force the solution to match the input solution
        // set_initial_solution(initfile, true, nil);

        ls.model.close();

        if (timelimit != nil) {
            ls.param.timeLimit = timelimit ;
        }

        if (initfile != nil) {
            set_initial_solution(initfile, false, nil);
        }

        ls.solve();
        output();
    }
}

function evaluate_flexibility(solution_file, silent) {
    test_file = "problems/demands/test_demands.json";
    scenarios = json.parse(test_file);

    summary = {};

    with(ls = localsolver.create()) {
        ls.param.verbosity = 0 ;
        model();
        add_objectives();
        ls.model.close();
        set_initial_solution(solution_file, false, nil);
        tours = aggregate_tours();
        summary["original_tours"] = tours;
    }

    sc_index = 0;
    n_rec = {};
    for [scenario in scenarios.keys()] {
        sc_index += 1;
        if (!silent)
            print("[EVAL] ", sc_index, "/", scenarios.count(), " ");
        demands = read_demands(scenarios[scenario], n, index_to_centre);

        with(ls = localsolver.create()) {
            ls.param.verbosity = 0 ;
            ls.param.timeLimit = 2 ;

            model();

            recourses = set_same_visits(solution_file);

            n_recourses <- sum[r in recourses](r) ; 
            minimize n_recourses;

            add_objectives();
            for [v in allowed_vehicles] {
                constraint used[v] == solution["vehicles_used"][v];
            }
            constraint n_recourses <= 5;

            ls.model.close();
            set_initial_solution(solution_file, false, nil);

            ls.solve();
            if (!silent)
                print("       ", ls.solution.status);

            summary[scenario]["status"] = ls.solution.status;

            if (ls.solution.status == "OPTIMAL" || ls.solution.status == "FEASIBLE") {
                if (!silent)
                    print("       cost=", round(total_costs.value));
                summary[scenario]["total_costs"] = total_costs.value;
                if (!silent)
                    print("       recourses=", n_recourses.value);
                summary[scenario]["recourses"] = n_recourses.value;
                n_rec.add(n_recourses.value);
                tours = aggregate_tours();
                summary[scenario]["tours"] = tours;
            }
            if (!silent)
                println();
        }
    }

    n_feasible = sum[scenario in scenarios.keys() : summary[scenario]["status"] == "FEASIBLE" || summary[scenario]["status"] == "OPTIMAL"](1);
    average_cost = sum[scenario in scenarios.keys() : summary[scenario]["status"] == "FEASIBLE" || summary[scenario]["status"] == "OPTIMAL"]
                            (summary[scenario]["total_costs"]) 
                    / n_feasible;

    if (!silent) {
        println("Feasible: ", n_feasible, "/", scenarios.count());
        println("Average cost: ", round(average_cost), "E");
        println("Recourses: ", n_rec);
    }

    json.dump(summary, "solutions/flexibility/summary.json");
    return sum[scenario in scenarios.keys()]
        (summary[scenario]["status"] == "INFEASIBLE" ? 2000 : summary[scenario]["total_costs"]) / scenarios.count();
}

function optimize_bilevel(solution_file) {

    original_score = 1340;
    best_score = original_score;
    best_solution = solution_file;
    for [d in 0...n_days] {
        print("[OPTIM] ", d, "/", n_days, " ");
        for [iteration in 0...10] {
            print("Current best score: ", round(best_score), "  |  improv=", round((original_score - best_score) / original_score * 100), "%");
            println("   |   file: ", best_solution);
            with(ls = localsolver.create()) {
                ls.param.verbosity = 0 ;
                ls.param.seed = random.create().next(0, 10000);
                demands = read_demands(problem["demands"], n, index_to_centre);
                model();
                add_objectives();
                days_to_fix = {};
                for [d2 in 0...n_days] {
                    if (d != d2) {
                        days_to_fix.add(d2);
                    }
                }
                set_initial_solution(best_solution, true, days_to_fix);
                ls.model.close();
                ls.param.timeLimit = 5 ;
                // set_initial_solution(best_solution, false, {d});
                ls.solve();

                if (ls.solution.status == "INFEASIBLE") {
                    continue;
                }

                tmp_file = "solutions/tmp" + iteration + ".json";
                write_solution(tmp_file);
            }

            score = evaluate_flexibility(tmp_file, true);

            if (score < best_score) {
                best_score = score;
                best_solution = tmp_file;
            }
        }
    }
    println("Best score: ", best_score);
    println("Best solution: ", best_solution);
}


function lns(solution_file) {
    original_score = 1340;
    best_score = original_score;
    best_solution = solution_file;
    for [d in 0...n_days] {
        print("[LNS] ", d, "/", n_days, " ");
        for [iteration in 0...10] {
            print("Current best score: ", round(best_score), "  |  improv=", round((original_score - best_score) / original_score * 100), "%");
            println("   |   file: ", best_solution);
            with(ls = localsolver.create()) {
                ls.param.verbosity = 0 ;
                ls.param.seed = random.create().next(0, 10000);
                demands = read_demands(problem["demands"], n, index_to_centre);
                model();
                add_objectives();
                days_to_fix = {};
                for [d2 in 0...n_days] {
                    if (d != d2) {
                        days_to_fix.add(d2);
                    }
                }
                set_initial_solution(best_solution, true, days_to_fix);
                ls.model.close();
                ls.param.timeLimit = 5 ;
                if (iteration == 0) {
                    set_initial_solution(best_solution, false, {d});
                }
                ls.solve();

                if (ls.solution.status == "INFEASIBLE") {
                    continue;
                }

                tmp_file = "solutions/tmp" + iteration + ".json";
                write_solution(tmp_file);
                score = ls.model.objectives[0].value;
            }


            if (score < best_score) {
                best_score = score;
                best_solution = tmp_file;
            }
        }
    }
    println("Best score: ", best_score);
    println("Best solution: ", best_solution);
}



function parse_args(args) {
    usage = "Usage: ls_solver problem_file -w week -i initsol -o outfile -t timelimit -e eval_file -b file_to_optimize";
    n_args = args.count();
    if (n_args == 0) {
        throw usage;
    }

    try {
        problem_file = args[0];

        i = 1;
        while (i < n_args) {
            arg = args[i];
            if (arg == "-w") {
                i += 1;
                week = args[i].toInt();
            } else if (arg == "-i") {
                i += 1;
                initfile = args[i];
            } else if (arg == "-o") {
                i += 1;
                outfile = args[i];
            } else if (arg == "-t") {
                i += 1;
                timelimit = args[i].toInt();
            } else if (arg == "-e") {
                i += 1;
                eval_file = args[i];
            } else if (arg == "-b") {
                i += 1;
                file_to_optimize = args[i];
            } else {
                println("Unknown argument: ", arg);
                throw usage;
            }

            if (i == n_args) {
                throw usage;
            }

            i += 1;
        }
    } catch(e) {
        throw usage;
    }
    
    improvement_limit = nil;
    if (week == nil) {
        week = 1;
    }
}

function main(args) {

    parse_args(args);

    input();

    lastBestValue = inf - 1 ;
    lastBestRunningTime = 0;
    lastImprovementRunningTime = 0;
    lastSolutionWritten = false;
    lastObjectiveValue = inf - 1;


    if (file_to_optimize != nil) {
        println("Optimizing solution ", file_to_optimize);
        optimize_bilevel(file_to_optimize);
        // lns(file_to_optimize);
        return;
    }

    if (eval_file != nil) {
        println("Evaluating file ", eval_file, " for problem ", problem_file);
        evaluate_flexibility(eval_file, false);
        return;
    } 
    
    println("Solving file ", problem_file, " for week ", week);
    if (initfile != nil) {
        println("   with initsol : ", initfile);
    }
    if (outfile != nil) {
        println("   with outfile : ", outfile);
    }
    if (timelimit != nil) {
        println("   with time limit : ", timelimit, " seconds");
    }

    solve();
    // solve_multi(10, 10);
}