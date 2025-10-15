from sympy.combinatorics.named_groups import SymmetricGroup
from sympy.combinatorics.permutations import Permutation
from itertools import combinations
from math import factorial, comb

# Function to count the number of non-isomorphic graphs on n vertices
def count_graphs(n):
    if n == 0:
        return 1

    # Prepare all possible vertex‐pairs and their indices
    vertex_pairs = list(combinations(range(n), 2))
    index_of_pairs = {pair: idx for idx, pair in enumerate(vertex_pairs)}

    # The symmetric group acting on n vertices
    group = SymmetricGroup(n)
    total = 0

    # Sum over the conjugacy classes (Burnside’s Lemma)
    for conjugacy_class in group.conjugacy_classes():
        f = next(iter(conjugacy_class))

        # Build the induced permutation on edges
        induced_index = []
        for (u, v) in vertex_pairs:
            new_u, new_v = f(u), f(v)
            image = tuple(sorted((new_u, new_v)))
            induced_index.append(index_of_pairs[image])

        induced_group = Permutation(induced_index)

        # Count how many edges are fixed
        cycle_lengths = [len(c) for c in induced_group.cyclic_form]
        moved_edges = sum(cycle_lengths)
        total_edges = comb(n, 2)
        fixed_edges = total_edges - moved_edges

        # Number of edge‐orbits under this group element
        edge_orbits = len(cycle_lengths) + fixed_edges

        # Each orbit can be present/absent → 2^orbits
        term = 2 ** edge_orbits

        total += len(conjugacy_class) * term

    # Divide by |G| = n! to count distinct graphs
    return total // factorial(n)