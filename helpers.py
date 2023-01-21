from math import inf


def compute_distance(obj1, obj2):
    if obj1 == obj2:
        return 0
    if isinstance(obj1, Node) and isinstance(obj2, Node):
        return obj1.distances[obj1.index][obj2.index]
    elif isinstance(obj1, Tree) and isinstance(obj2, Node):
        if obj1.weight_type == "wpgma":
            left_distance = obj1.left_branch.compute_distance(obj2)
            right_distance = obj1.right_branch.compute_distance(obj2)
            total_distance = (left_distance + right_distance) / 2
            return total_distance
        else:
            left_distance = obj1.left_branch.compute_distance(obj2) * obj1.left_branch.size
            right_distance = obj1.right_branch.compute_distance(obj2) * obj1.right_branch.size
            total_distance = (left_distance + right_distance) / (obj1.left_branch.size + obj1.right_branch.size)
            return total_distance
    elif isinstance(obj1, Node) and isinstance(obj2, Tree):
        if obj1.weight_type == "wpgma":
            left_distance = obj2.left_branch.compute_distance(obj1)
            right_distance = obj2.right_branch.compute_distance(obj1)
            total_distance = (left_distance + right_distance) / 2
            return total_distance
        else:
            left_distance = obj2.left_branch.compute_distance(obj1) * obj2.left_branch.size
            right_distance = obj2.right_branch.compute_distance(obj1) * obj2.right_branch.size
            total_distance = (left_distance + right_distance) / (obj2.left_branch.size + obj2.right_branch.size)
            return total_distance
    elif isinstance(obj1, Tree) and isinstance(obj2, Tree):
        if obj1.weight_type == "wpgma":
            left_distance = obj1.compute_distance(obj2.left_branch)
            right_distance = obj1.compute_distance(obj2.right_branch)
            total_distance = (left_distance + right_distance) / 2
            return total_distance
        else:
            min_tree, max_tree = sorted([obj1, obj2], key=lambda x: x.depth)
            left_distance = min_tree.compute_distance(max_tree.left_branch) * max_tree.left_branch.size
            right_distance = min_tree.compute_distance(max_tree.right_branch) * max_tree.right_branch.size
            total_distance = (left_distance + right_distance) / (max_tree.left_branch.size + max_tree.right_branch.size)
            return total_distance
    else:
        raise TypeError("Only supported between trees and nodes")


class Node:
    def __init__(self, node_name, distance_info):
        self.node_name = node_name
        self.distance_info = distance_info
        self.distances = distance_info[0]
        self.names = distance_info[1]
        self.index = self.names.index(node_name)
        self.size = 1
        self.weight_type = distance_info[2]

    def __repr__(self):
        return self.node_name

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        if not isinstance(other, Node):
            return False
        assert self.distance_info == other.distance_info, \
            "Can only compare Nodes from the same distance info"
        return self.node_name == other.node_name

    def compute_distance(self, other):
        if isinstance(other, Node):
            return self.distances[self.index][other.index]
        elif isinstance(other, Tree):
            return other.compute_distance(self)
        else:
            raise TypeError("Only supported between trees and nodes")

    def __add__(self, other):
        if isinstance(other, Node):
            half_distance = self.compute_distance(other)/2
            return Tree(self, half_distance, other, half_distance, self.distance_info)
        elif isinstance(other, Tree):
            return other + self
        else:
            raise TypeError("Only supported between trees and nodes")


class Tree:
    def __init__(self, left_branch, left_branch_distance, right_branch, right_branch_distance, distance_info):
        if isinstance(left_branch, Node) and isinstance(right_branch, Node):
            self.left_branch_distance = left_branch_distance
            self.right_branch_distance = right_branch_distance
            self.left_branch = left_branch
            self.right_branch = right_branch

        elif isinstance(left_branch, Node) and isinstance(right_branch, Tree):
            self.left_branch_distance = right_branch_distance
            self.right_branch_distance = left_branch_distance
            self.left_branch = right_branch
            self.right_branch = left_branch

        elif isinstance(left_branch, Tree) and isinstance(right_branch, Node):
            self.left_branch_distance = left_branch_distance
            self.right_branch_distance = right_branch_distance
            self.left_branch = left_branch
            self.right_branch = right_branch

        else:
            (left_branch, left_branch_distance), (right_branch, right_branch_distance) = \
                sorted([(left_branch, left_branch_distance), (right_branch, right_branch_distance)],
                       key=lambda x: x[1], reverse=True)
            self.left_branch_distance = left_branch_distance
            self.right_branch_distance = right_branch_distance
            self.left_branch = left_branch
            self.right_branch = right_branch

        self.distance_info = distance_info
        self.size = left_branch.size + right_branch.size
        self.weight_type = distance_info[2]

    def __repr__(self):
        return f"({self.left_branch}: {self.left_branch_distance}, {self.right_branch}: {self.right_branch_distance})"

    def __str__(self):
        return self.__repr__()

    def __add__(self, other):
        if not isinstance(other, Tree):
            if not isinstance(other, Node):
                raise TypeError("Only supported between trees and nodes")

        if isinstance(other, Node):
            if self.weight_type == "wpgma":
                distance = self.compute_distance(other)
                half_distance = distance/2
                self_distance = half_distance - self.depth
                return Tree(self, self_distance, other, half_distance, self.distance_info)
            else:
                left_distance = self.left_branch.compute_distance(other) * self.left_branch.size
                right_distance = self.right_branch.compute_distance(other) * self.right_branch.size
                total_distance = (left_distance + right_distance) / (self.left_branch.size + self.right_branch.size)
                half_distance = total_distance/2
                self_distance = half_distance - self.depth
                return Tree(self, self_distance, other, half_distance, self.distance_info)

        if isinstance(other, Tree):
            if self.weight_type == "wpgma":
                distance = self.compute_distance(other)
                half_distance = distance/2
                self_distance = half_distance - self.depth
                other_distance = half_distance - other.depth
                return Tree(self, self_distance, other, other_distance, self.distance_info)
            else:
                min_tree, max_tree = sorted([self, other], key=lambda x: x.depth)
                left_distance = min_tree.compute_distance(max_tree.left_branch) * max_tree.left_branch.size
                right_distance = min_tree.compute_distance(max_tree.right_branch) * max_tree.right_branch.size
                total_distance = (left_distance + right_distance) / (max_tree.left_branch.size + max_tree.right_branch.size)
                half_distance = total_distance/2
                min_tree_distance = half_distance - min_tree.depth
                max_tree_distance = half_distance - max_tree.depth
                return Tree(min_tree, min_tree_distance, max_tree, max_tree_distance, self.distance_info)

    def __eq__(self, other):
        if not isinstance(other, Tree):
            return False
        return all([(self.left_branch == other.left_branch and self.right_branch == other.right_branch) or (self.left_branch == other.right_branch and self.right_branch == other.left_branch),
                    self.left_branch_distance == other.left_branch_distance,
                    self.right_branch_distance == other.right_branch_distance])

    @property
    def depth(self):
        if isinstance(self.left_branch, Node):
            return self.left_branch_distance
        else:
            return self.left_branch_distance + self.left_branch.depth

    def compute_distance(self, other):
        if not isinstance(other, Tree):
            if not isinstance(other, Node):
                raise TypeError("Only supported between trees and nodes")

        if isinstance(other, Node):
            if self.weight_type == "wpgma":
                left_distance = self.left_branch.compute_distance(other)
                right_distance = self.right_branch.compute_distance(other)
                total_distance = (left_distance + right_distance) / 2
                return total_distance
            else:
                left_distance = self.left_branch.compute_distance(other) * self.left_branch.size
                right_distance = self.right_branch.compute_distance(other) * self.right_branch.size
                total_distance = (left_distance + right_distance) / (self.left_branch.size + self.right_branch.size)
                return total_distance

        if isinstance(other, Tree):
            if self.weight_type == "wpgma":
                left_distance = self.compute_distance(other.left_branch)
                right_distance = self.compute_distance(other.right_branch)
                total_distance = (left_distance + right_distance) / 2
                return total_distance
            else:
                min_tree, max_tree = sorted([self, other], key=lambda x: x.depth)
                left_distance = min_tree.compute_distance(max_tree.left_branch) * max_tree.left_branch.size
                right_distance = min_tree.compute_distance(max_tree.right_branch) * max_tree.right_branch.size
                total_distance = (left_distance + right_distance) / (max_tree.left_branch.size + max_tree.right_branch.size)
                return total_distance


def convert_to_nodes_correct(distance_info):
    return [Node(x, distance_info) for x in distance_info[1]]


def merge_best_pair_correct(list_elements):
    list_distances = [compute_distance(x, y) if x != y else inf for x in list_elements for y in list_elements]
    min_nonzero = min(list_distances)
    index_min = list_distances.index(min_nonzero)
    row = index_min % len(list_elements)
    column = int(index_min // len(list_elements))
    merged = list_elements[column] + list_elements[row]
    return [merged] + [elem for index, elem in enumerate(list_elements) if index not in [row, column]]


def build_the_tree_correct(distance_info):
    list_elements = convert_to_nodes_correct(distance_info)
    while len(list_elements) > 1:
        list_elements = merge_best_pair_correct(list_elements)
    return list_elements[0]


def main():

    matrix_dist_2 = [
            [0, 17, 21, 31, 23],
            [17, 0, 30, 34, 21],
            [21, 30, 0, 28, 39],
            [31, 34, 28, 0, 43],
            [23, 21, 39, 43, 0]
        ]

    nodes_ = ["a", "b", "c", "d", "e"]
    weight_1 = "upgma"

    distance_info4 = matrix_dist_2, nodes_, weight_1

    a_node = Node("a", distance_info4)
    b_node = Node("b", distance_info4)

    c_node = Node("c", distance_info4)
    d_node = Node("d", distance_info4)
    e_node = Node("e", distance_info4)

    first_tree = a_node + b_node
    second_tree = c_node + d_node
    third_tree = second_tree + e_node
    third_tree_alternative = e_node + second_tree
    last_tree = first_tree + third_tree
    #print(first_tree)
    #print(second_tree)
    #print(third_tree)
    #print(third_tree_alternative)
    #print(last_tree)

    print(convert_to_nodes_correct(distance_info4))
    print(merge_best_pair_correct(convert_to_nodes_correct(distance_info4)))
    print(build_the_tree_correct(distance_info4))
    tree_one = a_node + b_node + e_node
    #print(tree_one)
    #print(compute_distance(tree_one, d_node))


if __name__ == "__main__":
    main()