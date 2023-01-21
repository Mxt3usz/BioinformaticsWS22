from math import inf

########################################################
############## Programming tasks #######################
########################################################

from helpers import compute_distance
from helpers import Node
from helpers import Tree

def tutorial():
    """
    In all the tasks we will use the same structure of the given data which can be seen as three parts:
    1. matrix of distances between nodes
    2. names of the nodes
    3. wpgma/upgma mode
    These three components are merged into a "distance info" tuple for our convenience.
    """
    matrix_dist_1 = [[0, 3, 12, 12, 9], [3, 0, 13, 13, 10], [12, 13, 0, 6, 7], [12, 13, 6, 0, 7], [9, 10, 7, 7, 0]]
    nodes_1 = ["a", "b", "c", "d", "e"]
    weight_1 = "wpgma"
    distance_info1 = matrix_dist_1, nodes_1, weight_1


    """
    You are provided with one helper function and two helper classes to perform the tree assembly:
    compute_distance - computes the distance between two objects, it works for both Nodes and Trees and for both wpgma and upgma
    modes.
    Node - a helper class to convert initial data into a node
    Tree - a tree representation class which can be printed in Newick's representation
    
    Lets start creating nodes and trees.
    """

    a = Node("a", distance_info1)
    b = Node("b", distance_info1)
    first_tree = a + b
    print("Nodes a and b", a, b)
    print("Distance between nodes a and b is", compute_distance(a, b))
    print(type(a))
    print("Result of merging two nodes a and is a tree", first_tree)
    print(type(a+b))

    c = Node("c", distance_info1)
    d = Node("d", distance_info1)
    second_tree = c + d
    print("Distance between nodes c and d is", compute_distance(c, d))
    print("Result of merging two nodes c and d is tree", second_tree)

    merged_trees = first_tree + second_tree
    print("Distance between trees  (a,b) and (c,d) is", compute_distance(first_tree, second_tree))
    print("Result of merging together two trees", merged_trees)


def convert_to_nodes(distance_info):
    """ Exercise 4 a
    Implement the function 'convert_to_nodes' which takes the distance information
    and converts it to a list of nodes. Keep in mind that the node names are located in
    distance_info[1]
    """
    return [Node(node,distance_info) for node in distance_info[1]]


def merge_best_pair(list_elements):
    """ Exercise 4 b
    Implement the function 'merge_best_pair' which takes a list of elements
    that can be both nodes and trees. It finds the best pair to merge based on
    the distance. Then it merges the two closest objects and returns a list with the merged
    object as well as the remaining ones.
    Use the 'compute_distance' function to determine the distance between two elements.
    """
    min_val = inf
    i,j = None,None
    for val1 in range(len(list_elements)):
        for val2 in range(val1+1,len(list_elements)):
            if compute_distance(list_elements[val1],list_elements[val2]) < min_val:
                min_val = compute_distance(list_elements[val1],list_elements[val2])
                i , j = val1,val2
    list_elements += [list_elements[i] + list_elements[j]]
    list_elements = [list_elements[ele] for ele in range(len(list_elements)) if ele != i and ele != j]
    return list_elements
    

def build_the_tree(distance_info):
    """ Exercise 4 c
    Implement the function 'build_the_tree' which takes the distance_info and
    outputs the final tree. Use your implementations of 'convert_to_nodes' and
    'merge_best_pair'
    """
    node_list = convert_to_nodes(distance_info)
    while len(node_list) != 1:
        node_list = merge_best_pair(node_list)
    return node_list[0]


if __name__ == '__main__':
    tutorial()