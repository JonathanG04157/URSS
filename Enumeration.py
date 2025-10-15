import networkx as nx
from itertools import combinations
from McKay_Latex import canonical_label_and_aut_group, adj_matrix
from PET_Latex import count_graphs

# Function to enumerate all non-isomorphic graphs on n vertices
def enumerate_graph(n):
    # Base case: n = 0 or 1 has only the empty graph
    if n <= 1:
        G = nx.empty_graph(n)
        return [G]

    # Build up from graphs on n-1 vertices
    old_graphs = enumerate_graph(n - 1)
    all_graphs = []

    for graph in old_graphs:
        V = list(graph.nodes())
        _, aut = canonical_label_and_aut_group(graph)

        seen = []
        # Try all ways of connecting the new vertex to subsets of V
        for r in range(len(V) + 1):
            for subset in combinations(V, r):
                # Skip subsets equivalent under automorphisms
                if any({gen[v] for v in subset} in seen for gen in aut):
                    continue

                seen.append(set(subset))

                # Add new vertex and edges
                graph_to_add = graph.copy()
                new_v = n - 1
                graph_to_add.add_node(new_v)
                for u in subset:
                    graph_to_add.add_edge(u, new_v)

                all_graphs.append(graph_to_add)

    # Filter out isomorphic duplicates via canonical labelling
    matrices = []
    all_graphs_iso = []

    for g in all_graphs:
        label, _ = canonical_label_and_aut_group(g)
        matrix = adj_matrix(g, label)

        if matrix not in matrices:
            matrices.append(matrix)
            all_graphs_iso.append(g)

    return all_graphs_iso