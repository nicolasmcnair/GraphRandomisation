import numpy as np
from random import randrange, choice
from collections import deque
from itertools import product, permutations, chain
try:
    from itertools import izip_longest
    zip_function = izip_longest
except ImportError:
    from itertools import zip_longest
    zip_function = zip_longest

def validate_adjacency_matrix(matrix):
    """
    Ensures that the supplied adjacency matrix is in the right format (square, numpy array with no negative weights),
    and will provide a Euler circuit (all vertices have edges, and each vertex has the same number of in/out edges).
    Returns True if the adjacency matrix is valid, otherwise raises an error.
    """
    # Check matrix is numpy array
    if type(matrix).__module__ != np.__name__:
        raise TypeError('Expected numpy matrix')
    # Check adjacency matrix is square
    elif matrix.ndim != 2 or (matrix.shape[0] != matrix.shape[1]):
        raise ValueError('Matrix is not square')
    # Check whether adjacency matrix is Eulerian (i.e., sum of each node row equals the sum for the column of that node)
    elif any([matrix[:,x].sum() != matrix[:,x].sum() for x in xrange(matrix.shape[0])]):
        raise ValueError('Matrix is not Eulerian')
    # Check that each node has at least one edge
    elif ~matrix.any(axis=0).any() or ~matrix.any(axis=1).any():
        raise ValueError('One or more nodes\\vertices possess no edges')
    # Check that there are no negative values in the adjacency matrix
    elif (matrix < 0).any():
        raise ValueError('No negative weights permitted')
    else:
        return True

def generate_adjacency_matrix_and_nodes(sequence_info):
    """
    Generate adjacency matrix from supplied sequence information. The argument sequence_info must be a dictionary with
    at least a 'conditions' key which has either a positive integer value representing the number of conditions,
    or a container consisting of the different conditions. Other possible keys include: 'order', with a positive integer
    value indicating the degree of ordering to control for, 'repetitions', with a positive integer value indicating the
    number of times to repeat each condition order, and 'omit_self_adjacencies', which is a boolean indicating whether
    to exclude conditions from being adjacent to each other in the sequence (N.B., this should only ever be used with
    1st-degree order.
    """
    sequence_info = dict({'order':1,'repetitions':1,'omit_self_adjacencies':False},**sequence_info)
    # Generate nodes based on conditions and order
    nodes = range(sequence_info['conditions']) if isinstance(sequence_info['conditions'],int) else range(len(sequence_info['conditions']))
    n_nodes = len(nodes)
    # Create adjacency matrix (accounting for order depth)
    if sequence_info['order'] == 1:
        adjacency_matrix = np.full((n_nodes,n_nodes), sequence_info['repetitions'], dtype='int32')
        if sequence_info['omit_self_adjacencies']:
            adjacency_matrix[np.diag_indices(n_nodes)] = 0
    elif sequence_info['order'] > 1:
        nodes = list(product(nodes, repeat=sequence_info['order']))
        if sequence_info['omit_self_adjacencies']:
            nodes = [node for node in nodes if not any(x == y for x,y in zip(node[1:],node[:-1]))]
        n_nodes = len(nodes)
        adjacency_matrix = np.array([[sequence_info['repetitions'] if nodes[i][1:] == nodes[j][:-1] else 0 for j in range(n_nodes)] for i in range(n_nodes)], dtype='int32')
    else:
        raise ValueError('sequence_info[\'order\'] must have a value of 1 or greater.')
    return adjacency_matrix, nodes

def count_euler_circuits(adjacency_matrix):
    """Return the number of euler circuits that can be produced from the supplied adjacency matrix."""
    from scipy.misc import factorial
    # Check that matrix is correct
    validate_adjacency_matrix(adjacency_matrix)
    degree_matrix = np.diag(adjacency_matrix.sum(axis=1,dtype='int64'))
    laplacian_matrix = (degree_matrix - adjacency_matrix)
    # Return euler count
    return int(np.linalg.det(laplacian_matrix[:-1,:-1])) * np.prod(factorial(degree_matrix.sum(axis=0) - 1))

def count_graph_sequences(adjacency_matrix):
    """Return the number of possible conditions sequences that can be produced from the supplied adjacency matrix."""
    from scipy.misc import factorial
    return count_euler_circuits(adjacency_matrix) * adjacency_matrix.sum() * (1.0 / np.prod(factorial(adjacency_matrix)))

def get_cardinality_vector(sequence):
    """Return the cardinality vector of the supplied condition sequence."""
    conditions = sorted(set(sequence))
    # Get location of each condition in the sequence
    index_matrix = [[idx for idx,item in enumerate(sequence[1:]) if item == condition] for condition in conditions]
    # Get difference of the gap between each instance of a condition and the optimal gap distance
    gap_matrix = [j for i in [[abs((index_matrix[condition][idx] - index_matrix[condition][idx-1]) - len(conditions)) for idx in xrange(1,len(index_matrix[condition]))] for condition in conditions] for j in i]
    return [gap_matrix.count(position) for position in xrange(max(gap_matrix) + 1)]

def compare_cardinality_vectors(vector_a,vector_b):
    """Compare two cardinality vectors and return False if the first vector is less uniform than the second vector, otherwise return True."""
    return next((x for x in [a-b for a,b in zip_function(vector_a,vector_b,fillvalue=0)][::-1] if x),1) > 0

def get_trials_in_sequence(sequence_info):
    """Return the total number of trials that would be produced using the supplied sequence information."""
    sequence_info = dict({'order':1,'repetitions':1,'omit_self_adjacencies':False},**sequence_info)
    return sequence_info['repetitions'] * get_trials_per_repeat(sequence_info) + sequence_info['order']

def get_trials_per_repeat(sequence_info):
    """Return the number of trials that would be produced per repetition using the supplied sequence information."""
    sequence_info = dict({'order':1,'repetitions':1,'omit_self_adjacencies':False},**sequence_info)
    n_nodes = sequence_info['conditions'] if isinstance(sequence_info['conditions'],int) else len(sequence_info['conditions'])
    return n_nodes * pow(n_nodes - sequence_info['omit_self_adjacencies'],sequence_info['order'])

def generate_sequence_from_adjacency_matrix(adjacency_matrix,nodes,start=None):
    """Generate condition sequence from supplied adjacency matrix and nodes, using optional starting node if specified."""
    # Check that matrix is correct
    validate_adjacency_matrix(adjacency_matrix)

    n_nodes = len(nodes)
    first_order = True if isinstance(nodes[0], int) else False

    # Create initial walk through the arborescence
    arborescence = np.zeros(adjacency_matrix.shape, dtype='int32')
    current_node = randrange(n_nodes) if start is None else start
    nodes_visited = {current_node}

    # Walk backwards through the adjacency matrix, selecting nodes at random
    while len(nodes_visited) < n_nodes:
        source = choice(np.nonzero(adjacency_matrix[:,current_node])[0])

        # Add source to arborescence if not already there
        if source not in nodes_visited:
            arborescence[source,current_node] = 1

        # Add source to nodes visited, then set as current_node
        nodes_visited.add(source)
        current_node = source

    remaining_edges = adjacency_matrix - arborescence
    out_edges = [None] * n_nodes

    for node in range(n_nodes):
        out_edges[node] = deque(np.append(np.random.permutation(list(chain.from_iterable([[idx] * repeat for idx,repeat in enumerate(remaining_edges[node,:])]))), np.nonzero(arborescence[node,:])[0]))
    
    # Set first node in sequence to root of arboresence
    seq = [np.where(np.sum(arborescence,axis=1) == 0)[0][0]]

    # Walk the circuit
    while any(out_edges):
        seq += [out_edges[seq[-1]].popleft()]
    
    # If order > 1, then piece together the condition sequence from the circuit
    if not first_order:
        seq = list(nodes[seq[0]]) + [x[-1] for x in [nodes[node] for node in seq[1:]]]
    return seq

def generate_sequence(sequence_info,n_sequences=1):
    """Generate condition sequence from supplied sequence info, using optional starting node if specified."""
    sequence_info = dict({'order':1,'repetitions':1,'omit_self_adjacencies':False,'start':None},**sequence_info)
    if n_sequences > 1:
        adjacency_matrix,nodes = generate_adjacency_matrix_and_nodes(sequence_info)
        return [generate_sequence_from_adjacency_matrix(adjacency_matrix,nodes,start=sequence_info['start']) for _ in xrange(n_sequences)]
    else:
        return generate_sequence_from_adjacency_matrix(*generate_adjacency_matrix_and_nodes(sequence_info),start=sequence_info['start'])

def generate_sequence_using_matrix(adjacency_matrix,nodes,n_sequences=1,start=None):
    """Return a list of condition sequences from supplied adjacency matrix and nodes, using optional starting node if specified."""
    if n_sequences > 1:
        return [generate_sequence_from_adjacency_matrix(adjacency_matrix,nodes,start) for _ in xrange(n_sequences)]
    else:
        return generate_sequence_from_adjacency_matrix(adjacency_matrix,nodes,start)

def generate_uniform_sequence_using_matrix(adjacency_matrix,nodes,n_sequences=100,start=None):
    """Generate a number of condition sequences from supplied adjacency matrix and nodes, and return the sequence that has the most uniformity."""
    current_sequence = generate_sequence_from_adjacency_matrix(adjacency_matrix,nodes,start)
    current_cardinality_vector = get_cardinality_vector(current_sequence)
    for _ in xrange(1,n_sequences):
        new_sequence = generate_sequence_from_adjacency_matrix(adjacency_matrix,nodes,start)
        new_cardinality_vector = get_cardinality_vector(new_sequence)
        if compare_cardinality_vectors(current_cardinality_vector,new_cardinality_vector):
            current_sequence,current_cardinality_vector = new_sequence,new_cardinality_vector
    return current_sequence

def generate_uniform_sequence(sequence_info,n_sequences=100):
    """Generate a number of condition sequences from supplied sequence info, and return the sequence that has the most uniformity."""
    sequence_info = dict({'order':1,'repetitions':1,'omit_self_adjacencies':False,'start':None},**sequence_info)
    return generate_uniform_sequence_using_matrix(*generate_adjacency_matrix_and_nodes(sequence_info),n_sequences=n_sequences,start=sequence_info['start'])
