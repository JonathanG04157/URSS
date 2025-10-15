# URSS — Enumeration and Generation of Simple Unlabelled Graphs

**Author:** Jonathan Glanfield

A concise, didactic implementation for enumerating and generating simple *unlabelled* graphs. The project demonstrates counting by symmetry (Pólya / Burnside), McKay-style canonical labelling for isomorphism testing, and a recursive generator that extends graphs by orbit representatives under automorphism groups.

## Contents
- `PET.py` — tools for counting unlabelled graphs using cycle-index / Pólya methods.  
- `McKay.py` — canonical labelling (equitable refinement, individualisation, search tree) and basic automorphism detection.  
- `Enumeration.py` — recursive generator producing non-isomorphic simple graphs on *n* vertices.  
- `Enumeration_and_Generation_of_Simple_Unlabelled_Graphs.pdf` — project report with theory and complexity analysis.

## Requirements
- Python 3.8+  
- `networkx`  
- `sympy`

Install dependencies:
```bash
pip install networkx sympy
```

## Quick usage
These modules expose functions rather than a CLI. Example in a Python REPL or small driver script:
```python
from Enumeration import enumerate_graph
G_list = enumerate_graph(4)   # list of non-isomorphic graphs on 4 vertices
print(len(G_list))            # should match theoretical counts (e.g. 11 for n=4)
```

## Notes
- Implementations are primarily educational and **not optimised** for large *n*. The PDF discusses theoretical bounds and practical behaviour.  
- For large-scale enumeration consider specialised tools such as Brendan McKay’s `nauty` / `Traces`.

## References
See the included PDF for references (Pólya, Burnside, McKay, Harary & Palmer, Erdős–Rényi) and full algorithmic detail.

