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

function fix(percentage) {
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
                    index = place["index"] - 1;

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
                    ramasse_init.add(place["index"] - n - 1);
                    if (force){
                        constraint ramasses[d][v][j] == place["index"]-1-n;
                        j += 1;
                    } else {
                        ramasses[d][v].value.add(place["index"] - n - 1);
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

    for [d in 0...n_days] {
        for [v in 0...m] {

            livraisons[d][v] <- list(n);
            ramasses[d][v] <- list(n_pdr);

            for [c in 0...n] {
                visits_l[d][v][c] <- contains(livraisons[d][v], c) ;
            }
            for [p in 0...n_pdr] {
                visits_r[d][v][p] <- contains(ramasses[d][v], p);
            }
            
            for [c in 0...n] {
                serve_a[d][v][c] <- int(0, demands["a"][c]) ;
                serve_f[d][v][c] <- frais[v] == "Oui" ? int(0, demands["f"][c]) : 0 ;
                // serve_f[d][v][c] <- int(0, demands["f"][c]);
                serve_s[d][v][c] <- int(0, demands["s"][c]) ;

                load[d][v][c] <- serve_a[d][v][c] + serve_f[d][v][c] + serve_s[d][v][c] ;

                constraint visits_l[d][v][c] * demands["a"][c] >= serve_a[d][v][c] ;
                constraint visits_l[d][v][c] * demands["f"][c] >= serve_f[d][v][c] ;
                constraint visits_l[d][v][c] * demands["s"][c] >= serve_s[d][v][c] ;
                // constraint visits[d][v][c] * (demands["a"][c] + demands["f"][c] + demands["s"][c]) 
                //         >= serve_a[d][v][c] + serve_f[d][v][c] + serve_s[d][v][c] ;
                // constraint visits[d][v][c] * (demands["a"][c] + demands["f"][c] + demands["s"][c]) >= load[d][v][c] ;

                // constraint visits_l[d][v][c] == iif(load[d][v][c] > 0, 1, 0) ; // dominance breaking (but makes it hard to find anything feasible)

                palettes[d][v][c] <- int(0, sizes[v]);
                constraint palettes[d][v][c] * max_palette_capacity >= serve_a[d][v][c] + serve_f[d][v][c];
            } 

            constraint sum(
                0...n,
                c => load[d][v][c]
            ) <= capacities[v] ;

            constraint sum[c in 0...n](palettes[d][v][c]) <= sizes[v] ;

            n_livraisons[d][v] <- count(livraisons[d][v]) ;
            n_ramasses[d][v] <- count(ramasses[d][v]) ;

            local liv <- n_livraisons[d][v] > 0 ;
            local ram <- n_ramasses[d][v] > 0 ;

            route_dist[d][v] <- 
                sum(1...n_livraisons[d][v], i => dist_matrix[livraisons[d][v][i - 1]][livraisons[d][v][i]])
                + sum(1...n_ramasses[d][v], i => dist_matrix[ramasses[d][v][i - 1] + n][ramasses[d][v][i] + n])
                + (liv ? dist_from_depot[livraisons[d][v][0]] : 0) 
                + (ram ? dist_to_depot[ramasses[d][v][n_ramasses[d][v]-1] + n] : 0) 
                + (liv && ram ? dist_matrix[livraisons[d][v][n_livraisons[d][v]-1]][ramasses[d][v][0] + n] : 0) 
                + (liv && !ram ? dist_to_depot[livraisons[d][v][n_livraisons[d][v]-1]] : 0) 
                + (!liv && ram ? dist_from_depot[ramasses[d][v][0] + n] : 0) 
                ;
        }
    }



    // Symmetry breaking
    // for [d in 0...n_days-1] {
    //     for [v in 0...m] {
    //         constraint n_livraisons[d][v] >= n_livraisons[v][d+1] ; 
    //     }
    // }
    // for [d in 0...n_days] {
    //     for [v in 2...9-1] {
    //         constraint n_livraisons[d][v] >= n_livraisons[v+1][d] ;
    //     }
    //     for [v in 10...m-1] {
    //         constraint n_livraisons[d][v] >= n_livraisons[v+1][d] ;
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
                constraint visits_r[d][0][pdr] ;
            }

            constraint sum[v in 1...m : frais[v] == "Oui"](visits_r[d][v][pdr]) == n_ramasses ; 
        }
    }


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

    total_distance <- sum[d in 0...n_days](
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
            if (livraisons[d][v].value.count() == 0 && ramasses[d][v].value.count() == 0) continue;

            key = d + ", " + v ;
            sol["tours"][key] = {};

            for [node in livraisons[d][v].value] {
                place = {};

                place["index"] = node+1;
                place["type"] = "livraison";

                name = matrix.rows[node + 1][0];
                place["name"] = name;

                delivery = {serve_a[d][v][node].value, serve_f[d][v][node].value, serve_s[d][v][node].value};
                place["delivery"] = delivery;

                pals = palettes[d][v][node].value;
                place["palettes"] = pals;

                sol["tours"][key].add(place);
            }

            for [node in ramasses[d][v].value] {
                place = {};

                place["index"] = node+1+n;
                place["type"] = "ramasse";

                name = matrix.rows[node + 1 + n][0];
                place["name"] = name;

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
            if (count(livraisons[d][v].value) == 0 && count(ramasses[d][v].value) == 0) continue;
            println("\tVehicle ", v, " (", vehicles.rows[v][0], ")");
            for [node in livraisons[d][v].value] {
                pal = palettes[d][v][node].value ;
                println("\t\t", pal, " Palettes to ", matrix.rows[node + 1][0]);
            }
            for [node in ramasses[d][v].value] {
                println("\t\t Ramasse at ", matrix.rows[node + 1][0]);
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