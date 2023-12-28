use csv;
use json;
use io;
use localsolver;

function input() {

    usage = "Usage : localsolver ls_solver.lsp outfile=<your output file name>";
    if (outfile == nil) throw usage;

    // TODO : ignore centres & max palette capa as argument
    // ignore_centres = {2, 12, 7, 14, 21, 10, 3, 11}; // ALL SEMI HEBDO
    // ignore_centres = {2, 12, 7, 14, 21}; // WEEK 1 + FREE
    // ignore_centres = {2, 12, 10, 3, 11}; // WEEK 2 + FREE
    ignore_centres = {7, 14, 21}; // WEEK 1
    // ignore_centres = {10, 3, 11}; // WEEK 2
    // ignore_centres = {11}; // FRONTON
    // ignore_centres = {}; // NONE
    max_palette_capacity = 800 ;
    
    centres = csv.parse("data/centres.csv");
    pdr = csv.parse("data/points_de_ramasse.csv");
    matrix = csv.parse("data/euclidean_matrix.csv");
    vehicles = csv.parse("data/vehicules.csv");

    current_tours = json.parse("data/tours_tournees_actuelles_w1.json");

    n = centres.nbRows - 1;
    m = vehicles.nbRows;
    n_pdr = pdr.nbRows;
    n_nodes = n + n_pdr;
    n_days = 5 ;

    for [node in 0...n * 3]{
        a = node < n ? node : 
            node < 2 * n ? node - n : node - 2 * n;
        dist_from_depot[node] = matrix.rows[0][a+2] ;
        dist_to_depot[node] = matrix.rows[a+1][1] ;
        for [node2 in 0...n_nodes] {
            b = node2 < n ? node2 : 
                node2 < 2 * n ? node2 - n : node2 - 2 * n;
            dist_matrix[node][node2] = matrix.rows[a+1][b+2];
        }
    }

    for [i in 0...n] {
        demands["a"][i] = centres.rows[i+1][4] * 1.15; // Add 15% for robustness
        demands["f"][i] = centres.rows[i+1][5] * 1.15; // Add 15% for robustness
        demands["s"][i] = centres.rows[i+1][6] * 1.15; // Add 15% for robustness
    }

    for [i in ignore_centres] {
        demands["a"][i-1] = 0;
        demands["f"][i-1] = 0;
        demands["s"][i-1] = 0;
    }

    for [i in 0...m] {
        row = vehicles.rows[i];
        capacities[i] = row[1];
        sizes[i] = row[2];
        frais[i] = row[4];
    }

    for [i in 0...n_pdr] {
        freqs_pdr[i] = pdr.rows[i][4];
    }

    // this is cheating but i want to try
    for [c in 0...n] {
        if (demands["f"][c] > capacities[3]) {
            demands["f"][c] = capacities[3] ;
        }
    }
}

function set_initial_solution(livraisons, ramasses) {

    // FIXME
    for [d in 0...n_days] {
        for [v in 0...m] {
            tour = current_tours["("+v+", "+d+")"] ;
            local liv_init = {} ;
            local ram_init = {} ;
            if (tour != nil) {
                for [i in 1...count(tour)-1] {
                    if (tour[i] - 1 >= n) {
                        ram_init.add(tour[i] - n - 1);
                        // constraint contains(ramasses[d][v], tour[i]-n-1);
                    } else {
                        liv_init.add(tour[i] - 1);
                        // constraint contains(livraisons[d][v], tour[i]-1);
                    }
                }
                // println(v, " ", d, " ", liv_init);
                // println(v, " ", d, " ", ram_init);
                constraint livraisons[d][v] == liv_init;
                constraint ramasses[d][v] == ram_init;
            } else {
                // constraint count(livraisons[d][v]) == 0;
                // constraint count(ramasses[d][v]) == 0;
            }
        }
    }
}

function model() {


    tour[vd in 0...m * n_days] <- list(n * 3);

    for [vd in 0...m * n_days] {

        v = vd % m ;

        // Visit everyone 3 times

            for [node in 0...n*3] {
                visits[vd][node] <- contains(tour[vd], node) ;
                // id[d][v][node] <- indexOf(tour[d][v], node);
            }
            for [node in 0...3] {
                constraint or(
                    and(visits[vd][node], visits[vd][node+n], visits[vd][node+2*n]),
                    and(!visits[vd][node], !visits[vd][node+n], !visits[vd][node+2*n])
                ) ;
            }

            // for [pdr in n...n_nodes] {
            //     for [c in 0...n] {
            //         constraint or(
            //                 !visits[d][v][c], 
            //                 !visits[d][v][pdr], 
            //                 id[d][v][pdr] >= id[d][v][c]
            //             ) ;
            //     }
            // }
            
        for [c in 0...n] {
            // serve_a[vd][c] <- float(0, demands["a"][c]) ;
            // serve_f[vd][c] <- frais[v] == "Oui" ? float(0, demands["f"][c]) : 0 ;
            // // serve_f[vd][c] <- float(0, demands["f"][c]);
            // serve_s[vd][c] <- float(0, demands["s"][c]) ;

            // load[vd][c] <- serve_a[vd][c] + serve_f[vd][c] + serve_s[vd][c] ;

            // constraint visits[vd][c] * demands["a"][c] >= serve_a[vd][c] ;
            // constraint visits[vd][c] * demands["f"][c] >= serve_f[vd][c] ;
            // constraint visits[vd][c] * demands["s"][c] >= serve_s[vd][c] ;
            // constraint visits[vd][c] * (demands["a"][c] + demands["f"][c] + demands["s"][c]) 
            //         >= serve_a[vd][c] + serve_f[vd][c] + serve_s[vd][c] ;
            // constraint visits[vd][c] * (demands["a"][c] + demands["f"][c] + demands["s"][c]) >= load[vd][c] ;

            // constraint visits_l[vd][c] == iif(load[vd][c] > 0, 1, 0) ; // dominance breaking (but makes it hard to find anything feasible)

            // palettes[vd][c] <- int(0, sizes[v]);
            // constraint palettes[vd][c] * max_palette_capacity >= serve_a[d][v][c] + serve_f[d][v][c];
        } 

            // constraint sum(
            //     0...n,
            //     c => load[d][v][c]
            // ) <= capacities[v] ;


        local c <- count(tour[vd]) ;

        route_dist[vd] <- 
            sum(1...c, i => dist_matrix[tour[vd][i - 1]][tour[vd][i]])
            + (c>0 ? dist_from_depot[tour[vd][0]] : 0) 
            + (c>0 ? dist_to_depot[tour[vd][c-1]] : 0) 
            ;
    }

    // set_initial_solution(livraisons, ramasses);

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

    // for [pdr in 0...n_pdr] {
    //     // Visit points de ramasse enough times each week
    //     constraint sum(
    //         0...n_days,
    //         d => sum(
    //             0...m,
    //             v => visits[d][v][n + pdr]
    //         )
    //     ) >= freqs_pdr[pdr] ;

    //     // Don't visit the same one twice a day
    //     for [d in 0...n_days] {
    //         constraint sum(0...m,
    //             v => visits[d][v][n + pdr]
    //         ) <= 1 ;
    //     }
    // }


    // Meet demands
    constraint cover[vd in 0...n_days * m](tour[vd]);

    total_distance <- sum(
        0...n_days * m,
        vd => route_dist[vd]
        ) ;

    minimize total_distance;
}

function param() {
    // set_initial_solution(tours);
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
    if (stats.runningTime - lastBestRunningTime > 999 && !lastSolutionWritten) {
        println(">>>>>>> No improvement during 5 seconds: writing current solution to file");
        write_solution();
        lastSolutionWritten = true ;
    }
}

function main(args) {
    outfile=args[0];
    lastBestValue = inf - 1 ;
    lastBestRunningTime = 0;
    lastSolutionWritten = false;

    input();

    with(ls = localsolver.create()) {
        ls.addCallback("ITERATION_TICKED", callback);

        model();

        ls.model.close();
        param();
        ls.solve();
    }
    output();
}