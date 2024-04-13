# routage-du-coeur
Optimisation des tournées de véhicules des Restos du Coeur de Haute-Garonne

## Usage

Solve a problem file with Hexaly :  
```localsolver models/ls_solver.lsp <problem file> <init solution file or "nil"> <output file path or "nil"> <time limit>```  

Solve a problem file with OR-Tools (CP-SAT solver) :  
```python solve.py <problem file> --initsol <init solution file> --outfile <output file path>```  
options :
- ```--initsol <init solution file>``` (set an initial solution)  
- ```--outfile <output file path>``` (write the solution to a specified json file)
- ```--time_limit, -t <int>``` (time limit in seconds)

Check a solution file :  
```python -m utils.tests <problem file> <solution file>```

Create all possible JSON problem files from combinations of demands, week assignments and allowed days in problems/ subdirectories, using parameters in data/params.json, with the given number of scenarios :  
```python -m utils.make_problems  <n_scenarios>```

Print a solution as txt file/html map/dashboard :  
```python -m utils.plots <solution file> <-t> <-m> <-d>```

Run tests :  
```pytest```

[Example Solution](solutions/example.txt)  
[Example map](solutions/example.html)  
[Example dashboard](solutions/example.pdf)  