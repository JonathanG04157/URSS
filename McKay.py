import networkx as nx

# Check if the partition is discrete (each cell has only one vertex)
def is_discrete(partition):
    return all(len(cell) == 1 for cell in partition)


# Equitable refinement of the partition
def refine(G, partition):
    P = [sorted(cell) for cell in partition]
    changed = True

    # Keep refining until partition stabilises (equitable)
    while changed:
        changed = False
        new_P = []

        for cell in P:
            if len(cell) <= 1:
                new_P.append(cell)
                continue

            buckets = {}

            # Assign vertices to buckets based on number of neighbours in other cells
            for v in cell:
                sig = tuple(len(set(G.neighbors(v)).intersection(other)) for other in P if other is not cell)
                buckets.setdefault(sig, []).append(v)

            # Split the cell if vertices have different number of neighbours in other cells
            if len(buckets) > 1:
                changed = True
                for group in sorted(buckets):
                    new_P.append(sorted(buckets[group]))
            else:
                new_P.append(cell)

        P = new_P

    return P

# Select a non-trivial cell from the partition
def select_nontrivial_cell(partition):
    for cell in partition:
        if len(cell) > 1:
            return cell
    return None

# Individualise a partition by a vertex v
def individualize(partition, v):
    new = []

    for cell in partition:
        if v in cell:
            rest = [u for u in cell if u != v]
            new.append([v])
            if rest:
                new.append(sorted(rest))
        else:
            new.append(cell)

    return new

# Convert a partition to an ordering
def partition_to_ordering(partition):
    return [v for cell in partition for v in cell]

# Apply a labelling to a graph, stored as an adjacency matrix
def adj_matrix(G, ordering):
    idx = {v: i for i, v in enumerate(ordering)}
    n = len(ordering)
    mat = [[0] * n for _ in range(n)]

    for i, v in enumerate(ordering):
        for u in set(G.neighbors(v)):
            mat[i][idx[u]] = 1

    return mat

# Main search
def search(G, part, aut_gens, best_order=None, best_mat=None):
    # Refine the current partition
    P = refine(G, part)

    # If fully discrete, we have a complete ordering
    if is_discrete(P):
        return partition_to_ordering(P)

    # Choose a cell to individualise further
    cell = select_nontrivial_cell(P)

    for v in sorted(cell):
        # Create a refined partition by individualising v
        Pi = individualize(P, v)

        # Recursively search with the new partition
        cand = search(G, Pi, aut_gens, best_order, best_mat)
        mat = adj_matrix(G, cand)


        if best_order is None:
            best_order, best_mat = cand, mat
        else:
            if mat == best_mat and cand != best_order:
                # New automorphism found
                gen = {best_order[i]: cand[i] for i in range(len(cand))}
                if gen not in aut_gens:
                    aut_gens.append(gen)
            elif mat > best_mat:
                # Lexicographically better labeling
                best_order, best_mat = cand, mat

    return best_order

# Canonical labelling + automorphism group
def canonical_label_and_aut_group(G):
    if G.nodes() == None or len(G.nodes()) == 0:
        return [], []

    init = list(G.nodes())
    gens = [{v: v for v in init}]
    lab = search(G, [init], gens)

    return lab, gens

# Isomorphism check
def isomorphic(G1, G2):
    lab1, _ = canonical_label_and_aut_group(G1)
    lab2, _ = canonical_label_and_aut_group(G2)
    A1, A2 = adj_matrix(G1, lab1), adj_matrix(G2, lab2)

    return A1 == A2