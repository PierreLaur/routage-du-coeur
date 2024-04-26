# routage-du-coeur

Optimisation des tournées de véhicules des Restos du Coeur de Haute-Garonne

## Usage

Solve a problem file :  
```python solve.py <problem file>```  
options :
- ```--solver ``` choose a solver : "hexaly" or "ortools" - default ortools  
- ```--week ``` choose the week (1 or 2) - default 1  
- ```--initsol ``` set an initial solution with a json solution file  
- ```--outfile ``` write the solution to a specified json file  
- ```--time_limit, -t ``` (time limit in seconds)

Check a solution file :  
```python -m utils.check_solution <problem file> <solution file>```

Validate a solution file (generate 100 instances, check that the solution is valid with few changes in every instance):  
```python validate.py <problem file> <solution file>```

Create all possible JSON problem files from combinations of demands, week assignments and allowed days in problems/ subdirectories, using parameters in data/params.json :  
```python -m utils.make_problems```

Print a solution as txt file/html map/dashboard :  
```python -m utils.plots <solution file> <-t> <-m> <-d>```

Run tests :  
```pytest```

[Example Solution](solutions/example.txt)  
[Example map](solutions/example.html)  
[Example dashboard](solutions/example.pdf)  