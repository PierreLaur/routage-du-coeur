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

```problem.py``` defines useful data structures, including dataclasses for problem instances (variations and scenarios) and for solutions, both of which can be written to JSON  

Create all possible JSON problem files from combinations of demands, week assignments and allowed days in problems/ subdirectories, using parameters in data/params.json, with the given number of scenarios :  
```python -m utils.make_problems  <n_scenarios>```

```tests.py``` defines a solution checking procedure. set the test instances in the test_instances variable and run ```pytest utils/tests.py``` to solve with time_limit=10 and check the resulting solutions.  
Run ```python -m utils.tests <problem> <solution>``` to test a specific problem and solution file.

Print the solution as txt file/html map/dashboard :  
```python -m utils.plots <solution file> <-t> <-m> <-d>```

[Example Solution](solutions/example.txt)  
[Example map](solutions/example.html)  
[Example dashboard](solutions/example.pdf)  