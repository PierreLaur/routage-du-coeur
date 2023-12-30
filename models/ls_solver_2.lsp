use csv;
use json;
use io;
use localsolver;
use random;

// TODO : doc

function input() {

    usage = "Usage : localsolver ls_solver.lsp <init_solution or 'nil'> <desired output file path> <week number (1 or 2)>";
    if (outfile == nil) throw usage;
    if (week == nil) throw usage;

    max_palette_capacity = 800 ;
    
    centres = csv.parse("data/centres.csv");
    pdr = csv.parse("data/points_de_ramasse.csv");
    matrix = csv.parse("data/euclidean_matrix.csv");
    vehicles = csv.parse("data/vehicules.csv");

    n = centres.nbRows - 1;
    m = vehicles.nbRows;
    n_pdr = pdr.nbRows;
    n_nodes = n + n_pdr;
    n_days = 5 ;

    for [node in 0...n_nodes] {
        dist_from_depot[node] = matrix.rows[0][node+2] ;
        dist_to_depot[node] = matrix.rows[node+1][1] ;
        for [node2 in 0...n_nodes] {
            dist_matrix[node][node2] = matrix.rows[node+1][node2+2];
        }
    }

    for [i in 0...n] {
        demands["a"][i] = ceil(centres.rows[i+1][4] * 1.15); // Add 15% for robustness
        demands["f"][i] = ceil(centres.rows[i+1][5] * 1.15); // Add 15% for robustness
        demands["s"][i] = ceil(centres.rows[i+1][6] * 1.15); // Add 15% for robustness
    }

    for [i in 0...n] {
        semaine_livraison = centres.rows[i+1][7] ;
        if(semaine_livraison != 0 && semaine_livraison != week) {
            demands["a"][i] = 0;
            demands["f"][i] = 0;
            demands["s"][i] = 0;
        };
    }

    for [i in 0...m] {
        row = vehicles.rows[i];
        capacities[i] = row[1];
        sizes[i] = row[2];
        frais[i] = row[4];
    }

    jours_map = {"Lundi": 0, "Mardi": 1, "Mercredi": 2, "Jeudi": 3, "Vendredi": 4};
    use_pl = {};
    for [i in 0...n_pdr] {
        jours = pdr.rows[i][3 + week].split(", ") ;
        j_de_ramasse[i] = {};
        for [k in 0...jours.count()] {
            
            if (jours[k].indexOf("(PL") != -1) {
                jours[k] = jours[k].replace("(PL)", "") ;
                use_pl[jours_map[jours[k]]] = i;
            }

            j_de_ramasse[i].add(jours_map[jours[k]]);
        }
    }
}


function set_initial_solution(force, specific_vd) {

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

            tour_init = {};
            if (!force) {
                tour[d][v].value.clear();
            }
            i = 0;
            for [place in tours_init[key]] {
                if (place["type"] == "livraison") {
                    index = place["index"] - 1;

                    tour_init.add(index);
                    if (force) {
                        constraint tour[d][v][i] == index;
                        constraint visits[d][v][index] == 1 ;
                        i += 1 ;
                        constraint palettes[d][v][index] == place["palettes"];
                        constraint serve_a[d][v][index] >= place["delivery"][0] ;
                        constraint serve_f[d][v][index] >= place["delivery"][1] ;
                        constraint serve_s[d][v][index] >= place["delivery"][2] ;
                    } else {
                        tour[d][v].value.add(index);
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
                    tour_init.add(place["index"] - 1);
                    if (force){
                        constraint tour[d][v][i] == place["index"]-1;
                        i += 1;
                    } else {
                        tour[d][v].value.add(place["index"] - 1);
                    }
                }
            }
            if (force) {
                constraint count(tour[d][v]) == tour_init.count();
            }
        }
    }
    
}

function model() {

    for [d in 0...n_days] {
        for [v in 0...m] {

            tour[d][v] <- list(n_nodes);

            for [node in 0...n_nodes] {
                visits[d][v][node] <- contains(tour[d][v], node) ;
                id[d][v][node] <- indexOf(tour[d][v], node);
            }

            for [pdr in n...n_nodes] {
                for [c in 0...n] {
                    constraint or(
                            !visits[d][v][c], 
                            !visits[d][v][pdr], 
                            id[d][v][pdr] >= id[d][v][c]
                        ) ;
                }
            }

            
            for [c in 0...n] {
                serve_a[d][v][c] <- int(0, demands["a"][c]) ;
                serve_f[d][v][c] <- frais[v] == "Oui" ? int(0, demands["f"][c]) : 0 ;
                serve_s[d][v][c] <- int(0, demands["s"][c]) ;

                load[d][v][c] <- serve_a[d][v][c] + serve_f[d][v][c] + serve_s[d][v][c] ;

                constraint visits[d][v][c] * demands["a"][c] >= serve_a[d][v][c] ;
                constraint visits[d][v][c] * demands["f"][c] >= serve_f[d][v][c] ;
                constraint visits[d][v][c] * demands["s"][c] >= serve_s[d][v][c] ;

                palettes[d][v][c] <- int(0, sizes[v]);
                constraint palettes[d][v][c] * max_palette_capacity >= serve_a[d][v][c] + serve_f[d][v][c];
            } 

            constraint sum(
                0...n,
                c => load[d][v][c]
            ) <= capacities[v] ;

            constraint sum[c in 0...n](palettes[d][v][c]) <= sizes[v] ;

            n_stops[d][v] <- count(tour[d][v]) ;

            route_dist[d][v] <- 
                (n_stops[d][v] > 0) ?
                sum(1...n_stops[d][v], i => 
                    dist_matrix[tour[d][v][i - 1]][tour[d][v][i]])
                    + dist_from_depot[tour[d][v][0]]
                    + dist_to_depot[tour[d][v][n_stops[d][v] - 1]] : 0
                
                ;
        }
    }

    // // Symmetry breaking
    // for [d in 0...n_days-1] {
    //     for [v in 0...m] {
    //         constraint n_stops[d][v] >= n_stops[d+1][v] ; 
    //     }
    // }
    // for [d in 0...n_days] {
    //     for [v in 2...9-1] {
    //         constraint n_stops[d][v] >= n_stops[d][v+1] ;
    //     }
    //     for [v in 10...m-1] {
    //         constraint n_stops[d][v] >= n_stops[d][v+1] ;
    //     }
    // }


    for [pdr in 0...n_pdr] {
        // Visit points de ramasse enough times each week
        for [d in 0...n_days] {

            n_ramasses = 0 ;
            for [j in j_de_ramasse[pdr]] {
                if (j==d) {
                    n_ramasses += 1 ;
                }
            }

            if (use_pl[d] == pdr) {
                n_ramasses -= 1 ;
                constraint visits[d][0][n + pdr] ;
            }

            constraint sum[v in 1...m : frais[v] == "Oui"](visits[d][v][n + pdr]) == n_ramasses ; 
        }
    }

    // constraint cover[d in 0...n_days][v in 0...m](tour[d][v]);


    // Meet demands
    for [c in 0...n] {
        constraint sum(
            0...n_days,
            d => sum(
                0...m,
                v => serve_a[d][v][c]
            )
        ) >= demands["a"][c];

        constraint sum(
            0...n_days,
            d => sum(
                0...m,
                v => serve_f[d][v][c]
            )
        ) >= demands["f"][c];

        constraint sum(
            0...n_days,
            d => sum(
                0...m,
                v => serve_s[d][v][c]
            )
        ) >= demands["s"][c];
    }

    total_distance <- sum
        [d in 0...n_days](
        sum[v in 0...m](route_dist[d][v])
    );

    minimize total_distance;
}

function param() {

    // lsVerbosity = 1;
    // lsTimeLimit = 120;
    // lsNbThreads = 32;
    // lsSeed = 1 ;
}


function write_solution() {
    if (outfile == nil) return;
    local outf = io.openWrite(outfile);

    sol = {};
    sol["total_distance"] = total_distance.value;
    sol["tours"] = {};

    for [d in 0...n_days] {
        for [v in 0...m] {
            if (tour[d][v].value.count() == 0) continue;

            key = d + ", " + v ;
            sol["tours"][key] = {};

            for [node in tour[d][v].value] {
                place = {};
                place["index"] = node+1;
                name = matrix.rows[node + 1][0];
                place["name"] = name;

                if (node < n) {
                    place["type"] = "livraison";

                    delivery = {serve_a[d][v][node].value, serve_f[d][v][node].value, serve_s[d][v][node].value};
                    place["delivery"] = delivery;

                    pals = palettes[d][v][node].value;
                    place["palettes"] = pals;
                } else {
                    place["type"] = "ramasse";

                }
                sol["tours"][key].add(place);
            }
        }
    }
    json.dump(sol, outf);
}

function output(){

    write_solution();

    for [d in 0...n_days] {
        println("DAY ", d);
        for [v in 0...m] {
            if (count(tour[d][v].value) == 0) continue;
            println("\tVehicle ", v, " (", vehicles.rows[v][0], ")");
            for [node in tour[d][v].value] {
                if (node < n) {
                    pal = palettes[d][v][node].value ;
                    println("\t\t", pal, " Palettes to ", matrix.rows[node + 1][0]);
                } else {
                    println("\t\t Ramasse at ", matrix.rows[node + 1][0]);
                }
            }
            println();
        }
    }
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

        if (initfile != nil) {
            set_initial_solution(false, nil);
        }
        param();

        ls.solve();
    }
    output();
}