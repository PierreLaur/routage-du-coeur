# routage-du-coeur
Optimisation des tournées de véhicules des Restos du Coeur de Haute-Garonne

## Usage

Solve (requires Hexaly solver) :  
```localsolver models/ls_solver.lsp <init solution file or "nil"> <output file path> <week number (1 or 2)>```

With OR-Tools (CP-SAT solver) :  
```python solve.py <week number (1 or 2)> --outfile <output file path>```  
options :
- ```--infile <init solution file>``` (set an initial solution)  
- ```--improve <init solution file>``` (set a file as both input and output)

Print the solution as txt file/html map/dashboard :  
```python utils/plots.py <solution file> <-t> <-m> <-d>```

[Example Solution](solutions/test.txt)  
[Example map](solutions/test.html)  
[Example dashboard](solutions/test.pdf)  