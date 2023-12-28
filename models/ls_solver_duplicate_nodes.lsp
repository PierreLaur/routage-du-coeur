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
    n_days = 5 ;

    n_pdr = 0;
    node = 0 ;
    exclusion_groups = {};
    for [i in 0...pdr.nbRows] {
        freqs_pdr[i] = pdr.rows[i][4];
        n_pdr += freqs_pdr[i];
        exclusion_groups[i] = {};
        while (node < n_pdr) {
            exclusion_groups[i].add(node);
            node += 1;
        } 
    }

    n_nodes = n + n_pdr;
    for [node in 0...n_nodes]{
        a = node < n ? node : decode_pdr(node-n)+n ;
        dist_from_depot[node] = matrix.rows[0][a+2] ;
        dist_to_depot[node] = matrix.rows[a+1][1] ;
        for [node2 in 0...n_nodes] {
            b = node2 < n ? node2 : decode_pdr(node2-n)+n ;
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


}

function set_initial_solution(livraisons, ramasses) {
    for [v in 0...m] {
        for [d in 0...n_days] {
            tour = current_tours["("+v+", "+d+")"] ;
            local liv_init = {} ;
            local ram_init = {} ;
            if (tour != nil) {
                for [i in 1...count(tour)-1] {
                    if (tour[i] - 1 >= n) {
                        ram_init.add(tour[i] - n - 1);
                        // constraint contains(ramasses[v][d], tour[i]-n-1);
                    } else {
                        liv_init.add(tour[i] - 1);
                        // constraint contains(livraisons[v][d], tour[i]-1);
                    }
                }
                // println(v, " ", d, " ", liv_init);
                // println(v, " ", d, " ", ram_init);
                constraint livraisons[v][d] == liv_init;
                constraint ramasses[v][d] == ram_init;
            } else {
                // constraint count(livraisons[v][d]) == 0;
                // constraint count(ramasses[v][d]) == 0;
            }
        }
    }
}

function decode_pdr(p) {
    i = 0;
    j = 0;
    while (i<=p) {
        j += 1;
        i += freqs_pdr[j - 1] ; 
    }
    return j - 1;
}

function model() {

    for [v in 0...m] {
        for [d in 0...n_days] {

            livraisons[v][d] <- list(n);
            ramasses[v][d] <- list(n_pdr);

            for [c in 0...n] {
                visits_l[v][d][c] <- contains(livraisons[v][d], c) ;
            }
            for [p in 0...n_pdr] {
                visits_r[v][d][p] <- contains(ramasses[v][d], p);
            }


            
            for [c in 0...n] {
                serve_a[v][d][c] <- float(0, demands["a"][c]) ;
                serve_f[v][d][c] <- frais[v] == "Oui" ? float(0, demands["f"][c]) : 0 ;
                // serve_f[v][d][c] <- float(0, demands["f"][c]);
                serve_s[v][d][c] <- float(0, demands["s"][c]) ;

                load[v][d][c] <- serve_a[v][d][c] + serve_f[v][d][c] + serve_s[v][d][c] ;

                constraint visits_l[v][d][c] * demands["a"][c] >= serve_a[v][d][c] ;
                constraint visits_l[v][d][c] * demands["f"][c] >= serve_f[v][d][c] ;
                constraint visits_l[v][d][c] * demands["s"][c] >= serve_s[v][d][c] ;
                // constraint visits[v][d][c] * (demands["a"][c] + demands["f"][c] + demands["s"][c]) 
                //         >= serve_a[v][d][c] + serve_f[v][d][c] + serve_s[v][d][c] ;
                // constraint visits[v][d][c] * (demands["a"][c] + demands["f"][c] + demands["s"][c]) >= load[v][d][c] ;

                // constraint visits_l[v][d][c] == iif(load[v][d][c] > 0, 1, 0) ; // dominance breaking (but makes it hard to find anything feasible)

                palettes[v][d][c] <- int(0, sizes[v]);
                constraint palettes[v][d][c] * max_palette_capacity >= serve_a[v][d][c] + serve_f[v][d][c];
            } 

            constraint sum(
                0...n,
                c => load[v][d][c]
            ) <= capacities[v] ;


            n_livraisons[v][d] <- count(livraisons[v][d]) ;
            n_ramasses[v][d] <- count(ramasses[v][d]) ;

            local liv <- n_livraisons[v][d] > 0 ;
            local ram <- n_ramasses[v][d] > 0 ;

            route_dist[v][d] <- 
                sum(1...n_livraisons[v][d], i => dist_matrix[livraisons[v][d][i - 1]][livraisons[v][d][i]])
                + sum(1...n_ramasses[v][d], i => 
                        dist_matrix[ramasses[v][d][i - 1] + n]
                                    [ramasses[v][d][i] + n])
                + (liv ? dist_from_depot[livraisons[v][d][0]] : 0) 
                + (ram ? dist_to_depot[ramasses[v][d][n_ramasses[v][d]-1] + n] : 0) 
                + (liv && ram ? dist_matrix[livraisons[v][d][n_livraisons[v][d]-1]]
                                            [ramasses[v][d][0] + n] : 0) 
                + (liv && !ram ? dist_to_depot[livraisons[v][d][n_livraisons[v][d]-1]] : 0) 
                + (!liv && ram ? dist_from_depot[ramasses[v][d][0] + n] : 0) 
                ;
        }
    }

    constraint partition[v in 0...m][d in 0...n_days](ramasses[v][d]);
    
    // Don't visit the same pdr twice in the same day
    for [d in 0...n_days] {
        for [g in exclusion_groups] {
            constraint sum(0...m, 
                v => sum(
                    0...g.count(),
                    p => visits_r[v][d][g[p]]
                )
            ) <= 1 ;
        }
    }

    // set_initial_solution(livraisons, ramasses);

    // Symmetry breaking
    // for [d in 0...n_days-1] {
    //     for [v in 0...m] {
    //         constraint n_livraisons[v][d] >= n_livraisons[v][d+1] ; 
    //     }
    // }
    // for [d in 0...n_days] {
    //     for [v in 2...9-1] {
    //         constraint n_livraisons[v][d] >= n_livraisons[v+1][d] ;
    //     }
    //     for [v in 10...m-1] {
    //         constraint n_livraisons[v][d] >= n_livraisons[v+1][d] ;
    //     }
    // }

    // // Visit points de ramasse enough times each week
    // for [pdr in 0...n_pdr] {
    //     constraint sum(
    //         0...n_days,
    //         d => sum(
    //             0...m,
    //             v => visits_r[v][d][pdr]
    //         )
    //     ) >= freqs_pdr[pdr] ;
    // }

    // Meet demands
    for [c in 0...n] {
        constraint sum(
            0...n_days,
            d => sum(
                0...m,
                v => serve_a[v][d][c]
            )
        ) >= demands["a"][c] ;

        constraint sum(
            0...n_days,
            d => sum(
                0...m,
                v => serve_f[v][d][c]
            )
        ) >= demands["f"][c] ;

        constraint sum(
            0...n_days,
            d => sum(
                0...m,
                v => serve_s[v][d][c]
            )
        ) >= demands["s"][c] ;
    }

    total_distance <- sum(
        0...n_days,
        d => sum(
            0...m,
            v => route_dist[v][d]
        )
    );

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

    outf.println("total_distance=", total_distance.value);

    outf.println("tours={");
    for [d in 0...n_days] {
        for [v in 0...m] {
            if (livraisons[v][d].value.count() == 0 && ramasses[v][d].value.count() == 0) continue;
            outf.print("\t'",v, ", ", d, "' : [\n\t");

            for [node in livraisons[v][d].value] {

                outf.print("\t{");

                outf.print("\n\t\t\t");
                outf.print("'index' : ", node+1, ",");

                outf.print("\n\t\t\t");
                outf.print("'type' : 'livraison',");

                outf.print("\n\t\t\t");
                name = matrix.rows[node + 1][0];
                outf.print("'name' : '", name, "',");

                outf.print("\n\t\t\t");
                delivery = "(" + serve_a[v][d][node].value + ", " + serve_f[v][d][node].value + ", " + serve_s[v][d][node].value + ")";
                outf.print("'delivery' : ", delivery, ",");

                outf.print("\n\t\t\t");
                pals = palettes[v][d][node].value;
                outf.print("'palettes' : ", pals, ",");

                outf.print("\n\t\t},\n\t");
            }

            for [node in ramasses[v][d].value] {
                actual_node = decode_pdr(node) + 1 + n;

                outf.print("\t{");

                outf.print("\n\t\t\t");
                outf.print("'index' : ", +actual_node, ",");

                outf.print("\n\t\t\t");
                outf.print("'type' : 'ramasse',");

                outf.print("\n\t\t\t");
                name = matrix.rows[actual_node][0];
                outf.print("'name' : '", name, "',");

                outf.print("\n\t\t},\n\t");
            }

            outf.println("],");
        }
    }
    outf.println("}");
}

function output(){

    write_solution();

    for [d in 0...n_days] {
        println("DAY ", d);
        for [v in 0...m] {
            if (count(livraisons[v][d].value) == 0 && count(ramasses[v][d].value) == 0) continue;
            println("\tVehicle ", v, " (", vehicles.rows[v][0], ")");
            for [node in livraisons[v][d].value] {
                pal = palettes[v][d][node].value ;
                println("\t\t", pal, " Palettes to ", matrix.rows[node + 1][0]);
            }
            for [node in ramasses[v][d].value] {
                println("\t\t Ramasse at ", matrix.rows[decode_pdr(node) + 1 + n][0]);
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