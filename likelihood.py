import math


########################################################
############## Programming tasks #######################
########################################################


"""
In this assignment you will be asked to implement the dynamic programming approach for the complete tree likelihood
computation
"""

"""
First lets agree on the tree representation you will be working with.
The given tree will be represented by the dictionary. The keys of the dictionary are the tree nodes.
The leaf nodes have corresponding nucleotides as values.
All the other nodes have a tuple of the corresponding child nodes. Each child node is represented by a tuple itself.
A child node is a tuple with two values. First value is the child node and the second value is the distance.
tree_example = {0: ((5, 1), (6, 2)),
                    5: ((1, 1), (2, 1)),
                    6: ((3, 1), (4, 2)),
                    1: "A",
                    2: "A",
                    3: "G",
                    4: "C"}
Tree example represents the tree from the first exercise.
Nodes 1, 2, 3 and 4 are the leaf nodes with the corresponding values A, A, G and C.
Node 5 has two child nodes 1 and 2 and the node 6 has two child nodes 3 and 4.
The root node 0 has child nodes 5 and 6.
In order to complete the programming assignment implement the following helper functions:
"""


def tree_to_order(tree):
    """ Exercise 2a
    This function takes the dictionary which represent the tree structure
    and returns a list with the order of nodes which can be computed in the dynamic programming fashion
    for the tree_example it can for example be:
    [1, 2, 3, 4, 5, 6, 0] or
    [2, 4, 1, 5, 3, 6, 0] and so on ...
    Most importantly for each node in the order list all child nodes have to be either leaf nodes or pre-computed
    i.e. found in the order list before.
    """
    order = []
    while tree != {}:
        not_computed = {}
        for key,val in tree.items():
            if type(val) is str:
                order += [key]
            elif val[0][0] in order and val[1][0] in order:
                order += [key]
            else:
                not_computed[key] = val
        tree = not_computed
    return order
    
tree_example = {0: ((5, 1), (6, 2)),
                    5: ((1, 1), (2, 1)),
                    6: ((3, 1), (4, 2)),
                    1: "A",
                    2: "A",
                    3: "G",
                    4: "C"}


def compute_p_self(alpha, t):
    """ Exercise 2b
    See formula (2) in the first exercise
    """
    return 1/4 * (1 + 3*math.exp(-4*alpha*t))
    

def compute_p_other(alpha, t):
    """ Exercise 2c
    See formula (3) in the first exercise
    """
    return 1/4 * (1 - math.exp(-4*alpha*t))


def init_l_matrix(tree):
    """ Exercise 2d
    We will use the L matrix to compute all the values in the dynamic programming fashion
    The L matrix is two-dimensional array
    The first dimension is related to nodes and has the size equal to the number of nodes in the tree.
    Keep the same index as the nodes in the tree. Do not use the previously calculated order here.
    The second dimension is related to 4 nucleotides A C G and T. It is important to keep this order of these for the
    tests.
    init_l_matrix has to return the corresponding array filled in with zeros except the corresponding leaf nodes.
    In the leaf nodes there has to be 1 for the corresponding nucleotide. See formula (4) in the first exercise.
    """
    L = []
    for i in range(len(tree)): # nodes
        if type(tree[i]) is tuple:
            L += [[0] * 4]
            continue
        L += [[]]
        for base in "ACGT":
            if tree[i] == base:
                L[i] += [1]
            else:
                L[i] +=[0]
    return L


def find_parent(curr_node,tree):
    f = False
    for key,val in tree.items():
        to_compute = []
        if type(val) is tuple:
            for tupl in val:
                to_compute += [tupl[0]]
                if tupl[0] == curr_node:
                    f = True
            if f :
                return key,val,to_compute + [key]


def compute_l_matrix(tree, alpha):
    """ Exercise 2e
    Compute the whole matrix using dynamic programming. See formula (5) in the first exercise
    Do not forget to use:
        tree_to_order
        compute_p_self
        compute_p_other
        init_l_matrix
    """
    order = tree_to_order(tree)
    L = init_l_matrix(tree)
    while order != []:
        if len(order) != 1 : #root needs to be computed
            parent,childs,computed = find_parent(order[0],tree)
        else:
            parent,childs,computed = order[0],tree[order[0]],order
        while computed != []: # remove computed values
            if computed[0] in order:
                order.remove(computed[0])
            computed = computed[1:]
        for parent_base in range(4): # L_parent ACGT
            calc = [0 for _ in range(len(childs))]
            for child in range(len(childs)): # loop through childs
                for base in range(4): # loop through ACGT
                    if parent_base == base:
                        calc[child] += compute_p_self(alpha,childs[child][1]) * L[childs[child][0]][base]
                    else:
                        calc[child] += compute_p_other(alpha,childs[child][1]) * L[childs[child][0]][base]
            product = 1
            for val in calc: # product of all childs
                product *= val
            L[parent][parent_base] = product
    return L


def likelihood_computation(tree, alpha):
    """ Exercise 2f
    Compute the likelihood of the tree
    See formula (6) in the first exercise
    """
    L = compute_l_matrix(tree,alpha)
    order = tree_to_order(tree)
    product = 0
    for l_0 in L[order[-1]]:
        product += 1/4 * l_0
    return product

t1 = {
            0: "A",
            1: ((0, 1.5), (8, 1.5)),
            2: ((1, 4.25), (5, 2.25)),
            3: "C",
            4: ((3, 3), (6, 3)),
            5: ((4, 0.5), (7, 3.5)),
            6: "G",
            7: "A",
            8: "C"
        }

print(compute_l_matrix(t1,0.03))
print(likelihood_computation(t1,0.03))
