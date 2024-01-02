# routage-du-coeur
Optimisation des tournées de véhicules des Restos du Coeur de Haute-Garonne

## Usage

With Hexaly :
```localsolver models/ls_solver.lsp <init solution file or "nil"> <output file path> <week number (1 or 2)>```

Alternatively, with OR-Tools :
```python solve.py --infile <init solution file> --outfile <output file path> <week number (1 or 2)>```