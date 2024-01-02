use csv;
use json;
use io;
use localsolver;
use random;

function input() {

    usage = "Usage : localsolver ls_solver.lsp <init_solution or 'nil'> <desired output file path> <week number (1 or 2)>";
    if (outfile == nil) throw usage;
    if (week == nil) throw usage;

    max_palette_capacity = 800 ;
    
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

    // Demands (Add 15% for robustness)
    for [i in 0...n] {
        centre = index_to_centre[i] ;
        demands["a"][i] = ceil(centres.rows[centre+1][4] * 1.15); 
        demands["f"][i] = ceil(centres.rows[centre+1][5] * 1.15); // Add 15% for robustness
        demands["s"][i] = ceil(centres.rows[centre+1][6] * 1.15); // Add 15% for robustness
    }

    // Pickup weights (per pickup)
    for [p in 0...n_pdr] {
        weight[p] = ceil(pdr.rows[p][6] * 1.15);
    }

    for [i in 0...m] {
        row = vehicles.rows[i];
        // Capacities in kg
        capacities[i] = row[1];
        // Sizes in number of pallets
        sizes[i] = row[2];
        frais[i] = row[4];
        // frais[i] = "Oui";
    }

    // When do we have to pickup each site ?
    jours_map = {"Lundi": 0, "Mardi": 1, "Mercredi": 2, "Jeudi": 3, "Vendredi": 4};
    use_pl = {};
    for [i in 0...n_pdr] {
        jours = pdr.rows[i][3 + week].split(", ") ;
        j_de_ramasse[i] = {};
        for [k in 0...jours.count()] {
            
            // PL can be written to specifically ask a pickup at p on day d with the PL
            // in which case use_pl[d] = p
            if (jours[k].indexOf("(PL") != -1) {
                jours[k] = jours[k].replace("(PL)", "") ;
                use_pl[jours_map[jours[k]]] = i;
            }

            j_de_ramasse[i].add(jours_map[jours[k]]);
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

            if (tours_init[key] == nil) 
                continue;

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
                        } else {
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
                        constraint palettes[d][v][index] == place["palettes"];
                        constraint serve_a[d][v][index] >= place["delivery"][0] ;
                        constraint serve_f[d][v][index] >= place["delivery"][1] ;
                        constraint serve_s[d][v][index] >= place["delivery"][2] ;
                    } else {
                        livraisons[d][v].value.add(index);
                        palettes[d][v][index].value = place["palettes"];

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

function model() {

    livraisons[d in 0...n_days][v in 0...m] <- list(n);
    ramasses[d in 0...n_days][v in 0...m] <- list(n_pdr);

    // Deliver quantity for each product type (ambiant/frais/surgelé)
    serve_a[d in 0...n_days][v in 0...m][c in 0...n] <- int(0, demands["a"][c]) ;
    serve_f[d in 0...n_days][v in 0...m][c in 0...n] <- frais[v] == "Oui" ? int(0, demands["f"][c]) : 0 ;
    serve_s[d in 0...n_days][v in 0...m][c in 0...n] <- int(0, demands["s"][c]) ;

    constraint cover[d in 0...n_days][v in 0...m](livraisons[d][v]);
    constraint cover[d in 0...n_days][v in 0...m](ramasses[d][v]);

    for [d in 0...n_days] {
        for [v in 0...m] {

            n_livraisons[d][v] <- count(livraisons[d][v]) ;
            n_ramasses[d][v] <- count(ramasses[d][v]) ;

            local liv <- n_livraisons[d][v] > 0 ;
            local ram <- n_ramasses[d][v] > 0 ;

            visits_l[d][v][c in 0...n] <- contains(livraisons[d][v], c) ;
            visits_r[d][v][p in 0...n_pdr] <- contains(ramasses[d][v], p);
            
            for [c in 0...n] {
                load[d][v][c] <- visits_l[d][v][c] * (serve_a[d][v][c] + serve_f[d][v][c] + serve_s[d][v][c]) ;
                palettes[d][v][c] <- int(0, sizes[v]);

                // Use enough pallets to carry A+F
                // S can be carried aside
                constraint palettes[d][v][c] * max_palette_capacity >= visits_l[d][v][c] * (serve_a[d][v][c] + serve_f[d][v][c]);
            } 

            // Delivery capacity constraints
            constraint sum[c in 0...n](palettes[d][v][c]) <= sizes[v] ;
            constraint sum[c in 0...n](load[d][v][c]) <= capacities[v] ;

            // Pickup capacity constraints
            constraint sum[p in 0...n_pdr](weight[p] * visits_r[d][v][p]) <= capacities[v];
            // We assume a single palette is used for each pickup
            constraint sum[p in 0...n_pdr](visits_r[d][v][p]) <= sizes[v];

            // Deliver first, then pickup
            route_dist[d][v] <- 
                + (liv ? dist_from_depot[livraisons[d][v][0]] 
                        + sum(1...n_livraisons[d][v], i => dist_matrix[livraisons[d][v][i - 1]][livraisons[d][v][i]]) 
                        
                        + (ram ? 
                            dist_matrix[livraisons[d][v][n_livraisons[d][v]-1]][ramasses[d][v][0] + n]
                            : dist_to_depot[livraisons[d][v][n_livraisons[d][v]-1]])

                        : (ram ? dist_from_depot[ramasses[d][v][0] + n] : 0)) 
                + (ram ? dist_to_depot[ramasses[d][v][n_ramasses[d][v]-1] + n] : 0) 
                        + sum(1...n_ramasses[d][v], i => dist_matrix[ramasses[d][v][i - 1] + n][ramasses[d][v][i] + n])
                ;
        }
    }

    // Pickup each site enough times each week
    for [p in 0...n_pdr] {
        for [d in 0...n_days] {

            // How many times do we have to pickup p on day d ?
            n_ramasses = 0 ;
            for [j in j_de_ramasse[p]] {
                if (j==d) {
                    n_ramasses += 1 ;
                }
            }

            // Use the PL for one pickup, if specified
            if (use_pl[d] == p) {
                n_ramasses -= 1 ;
                constraint visits_r[d][0][p] ;
            }

            constraint sum[v in 1...m : frais[v] == "Oui"](visits_r[d][v][p]) == n_ramasses ; 
        }
    }

    // Meet demands for each product type
    for [c in 0...n] {
        constraint sum[d in 0...n_days][v in 0...m](visits_l[d][v][c] * serve_a[d][v][c]) >= demands["a"][c];
        constraint sum[d in 0...n_days][v in 0...m](visits_l[d][v][c] * serve_f[d][v][c]) >= demands["f"][c];
        constraint sum[d in 0...n_days][v in 0...m](visits_l[d][v][c] * serve_s[d][v][c]) >= demands["s"][c];
    }

    total_distance <- sum[d in 0...n_days][v in 0...m](route_dist[d][v]);

    minimize total_distance;
}

function write_solution() {
    if (outfile == nil) return;
    local outf = io.openWrite(outfile);

    sol = {};
    sol["total_distance"] = total_distance.value;
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

                pals = palettes[d][v][node].value;
                place["palettes"] = pals;

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

    write_solution();

    println("\nTotal distance : ", total_distance.value/1000);

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

    with(ls = localsolver.create()) {
        ls.addCallback("ITERATION_TICKED", callback);

        model();

        // set_initial_solution(true);
        // fix(0.5);

        ls.model.close();
        if (timelimit != nil) {
            ls.param.timeLimit = timelimit ;
            // ls.param.timeLimit = 2;
        }

        if (initfile != nil) {
            set_initial_solution(false, nil);
        }

        ls.solve();
        output();
    }
}