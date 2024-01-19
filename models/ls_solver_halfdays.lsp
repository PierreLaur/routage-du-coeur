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

    max_palette_capacity = 800 ;
    demi_palette_capacity = 0.45 * max_palette_capacity ;
    
    norvegienne_capacity = 40 ;
    n_norvegiennes = 10 ;

    // The coefficients of the linear proxy for tour duration
    // Determined in fit_duration_proxy.ipynb
    duration_coefficients = {0.0435, 560.0} ;
    wait_at_centres = 15 ; // in minutes
    wait_at_pdrs = 40 ;
    max_duration = 200 ; // in minutes
    
    centres = csv.parse("data/centres.csv");
    // centres = csv.parse("data/centres_good_assignment.csv");
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
        demands["a"][i] = ceil(centres.rows[centre+1][4] * 1.3); 
        demands["f"][i] = ceil(centres.rows[centre+1][5] * 1.3); // Add 15% for robustness
        demands["s"][i] = ceil(centres.rows[centre+1][6] * 1.3);
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
        // frais[i] = "Oui";
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

    // Deliver quantity for each product type (ambiant/frais/surgelé)
    serve_a[h in 0...2][d in 0...n_days][v in 0...m][c in 0...n] <- (v == 1 || v >= 5) ? 0 : int(0, demands["a"][c]) ;
    serve_f[h in 0...2][d in 0...n_days][v in 0...m][c in 0...n] <- (v == 1 || v >= 5) ? 0 : 
                                                        (frais[v] == "Oui" ? int(0, demands["f"][c]) : 0) ;
    serve_s[h in 0...2][d in 0...n_days][v in 0...m][c in 0...n] <- (v == 1 || v >= 5) ? 0 : int(0, demands["s"][c]) ;

    constraint cover[h in 0...2][d in 0...n_days][v in 0...m](livraisons[h][d][v]);
    constraint cover[h in 0...2][d in 0...n_days][v in 0...m](ramasses[h][d][v]);

    for [h in 0...2] {
        for [d in 0...n_days] {
            for [v in 0...m] {

                n_livraisons[h][d][v] <- count(livraisons[h][d][v]) ;
                n_ramasses[h][d][v] <- count(ramasses[h][d][v]) ;

                local liv <- n_livraisons[h][d][v] > 0 ;
                local ram <- n_ramasses[h][d][v] > 0 ;

                visits_l[h][d][v][c in 0...n] <- contains(livraisons[h][d][v], c) ;
                visits_r[h][d][v][p in 0...n_pdr] <- contains(ramasses[h][d][v], p);
                
                for [c in 0...n] {
                    load[h][d][v][c] <- visits_l[h][d][v][c] * (serve_a[h][d][v][c] + serve_f[h][d][v][c] + serve_s[h][d][v][c]) ;
                    palettes[h][d][v][c] <- ceil(visits_l[h][d][v][c] * serve_a[h][d][v][c] / max_palette_capacity) ;
                    demi_palettes[h][d][v][c] <- ceil(visits_l[h][d][v][c] * serve_f[h][d][v][c] / demi_palette_capacity);

                    norvegiennes[h][d][v][c] <- int(0, n_norvegiennes);
                    demi_palettes_s[h][d][v][c] <- frais[v] == "Oui" ? 
                                        ceil(
                                            (visits_l[h][d][v][c] * serve_s[h][d][v][c] - 
                                            norvegiennes[h][d][v][c] * norvegienne_capacity) / demi_palette_capacity
                                            ) : 
                                        0;

                    constraint visits_l[h][d][v][c] * serve_s[h][d][v][c] <= 
                                    norvegiennes[h][d][v][c] * norvegienne_capacity +
                                    demi_palettes_s[h][d][v][c] * demi_palette_capacity ;
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
                            duration_coefficients[0] * 
                                    (route_dist[h][d][v] - (ram ? 
                                        dist_to_depot[ramasses[h][d][v][n_ramasses[h][d][v]-1] + n] : 
                                        (liv ? dist_to_depot[livraisons[h][d][v][n_livraisons[h][d][v]-1]] : 0))) + 
                            duration_coefficients[1] + 
                            wait_at_centres * 60 * n_livraisons[h][d][v] +
                            wait_at_pdrs * 60 * n_ramasses[h][d][v]
                            ;

                constraint tour_duration[h][d][v] <= max_duration * 60 ;
                constraint n_livraisons[h][d][v] + n_ramasses[h][d][v] <= 4; // This is basically redundant
            }

            constraint sum[c in 0...n][v in 0...m](norvegiennes[h][d][v][c]) <= n_norvegiennes ;

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
        for [v in 0...m] {

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
                if (place["type"] == "livraison") {
                    index = - 1 ;
                    for [c in 0...n] {
                        if (index_to_centre[c] == place["index"] - 1) {
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
                        constraint visits_l[d][v][index] == 1 ;
                        i += 1 ;
                        // constraint norvegiennes[d][v][index] >= place["norvegiennes"];
                        if (place["delivery"][0] == 0) {
                            constraint serve_a[d][v][index] == 0 ;
                        }
                        if (place["delivery"][1] == 0) {
                            constraint serve_f[d][v][index] == 0 ;
                        }
                        if (place["delivery"][2] == 0) {
                            constraint serve_s[d][v][index] == 0 ;
                        }
                    } else {
                        livraisons[d][v].value.add(index);
                        norvegiennes[d][v][index].value = place["norvegiennes"];

                        if (place["delivery"][0] > demands["a"][index]) {
                            serve_a[d][v][index].value = 0;
                        } else {
                            serve_a[d][v][index].value = place["delivery"][0];
                        }

                        if (frais[v] == "Oui") {
                            if (place["delivery"][1] > demands["f"][index]) {
                                serve_f[d][v][index].value = 0;
                            } else {
                                serve_f[d][v][index].value = place["delivery"][1];
                            }
                        }

                        if (place["delivery"][2] > demands["s"][index]) {
                            serve_s[d][v][index].value = 0;
                        } else {
                            serve_s[d][v][index].value = place["delivery"][2];
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

function set_current_tours() {
    /* Constrains the solution to follow the current tours */
    current_tours = json.parse("data/tours_tournees_actuelles_w" + week + ".json");
    for [d in 0...n_days] {
        for [v in 0...m] {
            key = "(" + d + ", " + v + ")" ;
            if (current_tours[key] == nil) continue;
            init_liv = {} ;
            init_ram = {} ;
            i = 0 ;
            j = 0 ;
            for [k in 1...current_tours[key].count()-1] {
                for [c in 0...n] {
                    if (index_to_centre[c] == current_tours[key][k] - 1) {
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


function check_solution_flexibility() {
    /*
    Checks the flexibility of the tours provided as input

    Simulates many demands from the given distribution, then for each one
    loosely constrains the solution to contain the same tours as the initial solution file
    and tries to find feasible solutions with a limited number of simple changes
    TODO
    */
}

function check_solution_robustness(n_runs, sigma) {
    /*
    Checks the robustness of the tours provided as input

    Simulates many demands from the given distribution, then for each one
    constrains the solution to be almost identical (we can only add products to existing tours, no new stops or tours)
    Prints how many feasible solutions it finds
    */

    if (initfile == nil) return;


    for [p in {"a", "f", "s"}] {
        for [c in 0...n] {
            original_demands[p][c] = demands[p][c] * (1/1.15) ;
        }
    }
    rng = random.create() ;

    feasible = 0 ;
    st = datetime.now() ;
    for [run in 0...n_runs] {
        println("Simulating run ", run, "/",n_runs,"... time=",datetime.now()-st,"s feasible=",feasible,"/",run) ;
        for [p in {"a", "f", "s"}] {
            for [c in 0...n] {
                demands[p][c] = max(1, ceil(original_demands[p][c] * rng.nextNormal(1, sigma))) ;
            }
        }

        with(ls = localsolver.create()) {
            // ls.addCallback("ITERATION_TICKED", callback);

            model();
            set_initial_solution(true, nil); // Force the tours to be identical to the given solution
            ls.model.close();
            set_initial_solution(false, nil); // Give the given solution as a hint
            ls.param.timeLimit = 3 ;
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


function write_solution() {
    if (outfile == nil) return;
    local outf = io.openWrite(outfile);

    sol = {};
    sol["total_distance"] = sum[d in 0...n_days][v in 0...m](route_dist[d][v].value);
    sol["fuel_consumption"] = fuel_consumption.value;
    sol["tours"] = {};

    for [d in 0...n_days] {
        for [v in 0...m] {
            if (livraisons[d][v].value.count() == 0 && ramasses[d][v].value.count() == 0) continue;

            key = d + ", " + v ;
            sol["tours"][key] = {};

            for [node in livraisons[d][v].value] {
                place = {};
                node_index = index_to_centre[node];

                place["index"] = node_index + 1 ;
                place["type"] = "livraison";

                name = matrix.rows[node_index + 1][0];
                place["name"] = name;

                delivery = {serve_a[d][v][node].value, serve_f[d][v][node].value, serve_s[d][v][node].value};
                place["delivery"] = delivery;

                place["palettes"] = {
                        palettes[d][v][node].value,
                        demi_palettes[d][v][node].value,
                        demi_palettes_s[d][v][node].value
                };
                place["norvegiennes"] = serve_s[d][v][node].value > demi_palettes_s[d][v][node].value * demi_palette_capacity ? 
                                            norvegiennes[d][v][node].value : 0;

                sol["tours"][key].add(place);
            }

            for [node in ramasses[d][v].value] {
                place = {};

                place["index"] = node + n_centres;
                place["type"] = "ramasse";

                name = matrix.rows[node + n_centres][0];
                place["name"] = name;

                sol["tours"][key].add(place);
            }
        }
    }
    json.dump(sol, outf);
}

function output(){

    // write_solution();

    // println("\nTotal distance : ", sum[d in 0...n_days][v in 0...m](route_dist[d][v].value)/1000);
    println("Total fuel consumption : ", fuel_consumption.value);

}

function callback(ls, cbType) {
    local stats = ls.statistics;
    local obj <- ls.model.objectives[0];
    if (ls.solution.status == "FEASIBLE" && obj.value < lastBestValue) {
        lastBestRunningTime = stats.runningTime;
        lastBestValue = obj.value;
        lastSolutionWritten = false ;
    }
    if (stats.runningTime - lastBestRunningTime > 5 && !lastSolutionWritten) {
        println(">>>>>>> No improvement during 5 seconds: writing current solution to file");
        write_solution();
        lastSolutionWritten = true ;
    }
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
    // check_solution_robustness(15, 0.3) ;

    with(ls = localsolver.create()) {
        // ls.addCallback("ITERATION_TICKED", callback);

        model();

        // fix(0.4);
        // set_current_tours();

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