import sys
import numpy as np

def find_all_trees(input_dir):
    # find all tree files in input_dir
    tree_files = os.listdir(input_dir)
    tree_files = [os.path.join(input_dir, t) for t in tree_files if t != '.DS_Store']
    n_trees = len(tree_files)
    return tree_files, n_trees


def is_tree(A, labels, filename):
    # a debugging function, which checks that the input
    # is a tree. Outputs bad labels
    n = A.shape[0]
    if sorted(np.sum(A, axis=0)) == [0] + [1]*(n-1):
        return True
    else:
        # print orphans
        for orphan in np.argwhere(np.sum(A, axis=0) == 0).reshape((-1,)):
            print labels[orphan], 'is an orphan'

        # print multiple parents
        for multi_parent in np.argwhere(np.sum(A, axis=0) > 1).reshape((-1,)):
            parents = [labels[g] for g in np.argwhere(A[    :, multi_parent]).reshape((-1,))]
            print labels[multi_parent], 'has multiple parents:',', '.join(parents)

        raise ValueError(filename + " is not a tree")

def distance(T1, T2):

    # main high-level function for computing the distance.
    # check shape
    assert T1.shape == T2.shape
    n = T1.shape[0]

    # make transitive closures
    T1trans = make_transitive_closure(T1, n)
    T2trans = make_transitive_closure(T2, n)

    # get intersection
    I = np.logical_and(T1trans, T2trans).astype(int)

    # -> tree
    Itree = dag_to_tree(I)
    Itrans = make_transitive_closure(Itree, n)

    s1 = np.sum(T1trans) - np.sum(np.logical_and(T1trans, Itrans))
    s2 = np.sum(T2trans) - np.sum(np.logical_and(T2trans, Itrans))
    D = s1 + s2

    # return the distance
    return D


def dag_to_tree(dag):

    def get_depth(i, depths, ancestors):
        if len(ancestors[i]) == 0:
            return 1
        else:
            return max([depths[a] for a in ancestors[i]]) + 1

    def get_root(D, ignore, n):
        Dtemp = np.empty_like(D)
        Dtemp[:] = D

        Dtemp[ignore, :] = 0
        Dtemp[:, ignore] = 0
        roots = np.argwhere(np.sum(Dtemp, axis=0) == 0).reshape((-1,))
        roots = [r for r in roots if r not in ignore]
        return roots[0]

    n = dag.shape[0]
    parents = get_parents(dag, n)
    ancestors = get_ancestors(parents, n)

    # get initial root
    x = get_root(dag, [], n)
    ignore, depths = [x], {x : 1}
    T = np.zeros((n, n), dtype=int)
    while len(ignore) < n:

        y = get_root(dag, ignore, n)

        # 2. get it's depth and store
        d = get_depth(y, depths, ancestors)
        depths[y] = d

        # 3. find any ancestor with depth = d - 1: any will do
        x = [a for a in ancestors[y] if depths[a] == d - 1]
        if len(x) > 1:
            print '    - warning: found two possible spanning trees: choosing one arbitrarily'
        x = x[-1]

        # 4. draw an edge twixt root and y
        T[x, y] = 1

        # 4. Remove y
        ignore.append(y)

    return T


def make_transitive_closure(M, n):
    parents = get_parents(M, n)

    Mdash = np.empty_like(M)
    Mdash[:] = M
    for seed in range(n):
        curr = seed
        while parents[curr] != []:
            for p in parents[curr]:
                Mdash[p, seed] = 1
                curr = p
    return Mdash


def get_ancestors(parents, n):

    node_ancestors = dict()
    for i in range(n):
        curr = i
        node_ancestors[curr] = []
        while parents[curr] != []:
            for p in parents[curr]:
                if p not in node_ancestors[i]:
                    node_ancestors[i].append(p)
                curr = p

    return node_ancestors


def get_parents(M, n):
    node_parents = dict()
    for i in range(n):
        if 1 in M[:, i]:
            node_parents[i] = list(np.argwhere(M[:, i] == 1).reshape((-1,)))
        else:
            node_parents[i] = []

    return node_parents


def make_bush(n, root_index):
    # makes a bush of size n
    bush = np.zeros((n, n), dtype=int)
    bush[root_index, :] = 1
    bush[root_index, root_index] = 0
    return bush


def read_tree(tree_file):
    # read tree file and compute label set
    labels = set()
    T = []
    for line in open(tree_file, 'r'):
        line = line.strip()
        parent, child = line.split(',')
        parent, child = parent.strip(), child.strip()

        T.append([parent, child])

        labels.add(parent)
        labels.add(child)

    return T, labels
    
def label_jaccard(labels):
    # compute the label jaccard
    n_trees = len(labels)
    Jaccard = np.zeros((n_trees, n_trees))
    for i in range(n_trees):
        for j in range(i+1, n_trees):
            intersection = labels[i].intersection(labels[j])
            union = labels[i].union(labels[j])
            Jaccard[i,j] = (len(intersection) / float(len(union)))

    return Jaccard

def add_missing_nodes(A, root_index):
    # find empty rows and add in nodes if needed
    to_insert = np.argwhere(np.sum(A, axis=0) == 0).reshape((-1,))
    if len(to_insert) > 1:
        print '    - warning: adding ' + str(len(to_insert)) + ' nodes to ' + tree_file
    for root in to_insert:
        if root != root_index:
            A[root_index, root] = 1
    return A


def tree_files_to_adjacency_matrices(tree_files):

    # first get the label set for each tree
    tree_label_sets, tree_labels = [], []
    for tree_file in tree_files:
        t, l = read_tree(tree_file)
        tree_labels.append(t)
        tree_label_sets.append(l)

    # so that we can compute the jaccard overlap
    Jaccard = label_jaccard(tree_label_sets)

    # now get all genres
    complete_labels = set()
    for l in tree_label_sets:
        complete_labels = complete_labels.union(l)

    # convert to a list so we can index
    complete_labels = sorted(list(complete_labels))
    n_labels = len(complete_labels)

    # Now go though the files again, building adjacency 
    # matrices indexed by complete_labels. Initialise
    # the root index to be None
    As = []
    root_index = None
    for T in tree_labels:
        # Adjacency matrix for this tree
        A = np.zeros((n_labels, n_labels), dtype=int)

        # convert parent, child labels to indices
        for parent_child in T:
            parent, child = parent_child      
            A[complete_labels.index(parent), complete_labels.index(child)] = 1

        # get/set root index
        root_index = tree_to_root(A, root_index)

        # Post-processing: add in any missing nodes
        A = add_missing_nodes(A, root_index)

        # check A is a tree (and report if not), store
        is_tree(A, complete_labels, tree_file)
        As.append(A)

    return As, root_index, Jaccard, complete_labels

def tree_to_root(A, root=None):
    # gets the root of a tree from the adjacency matrix, and
    # does some checks
    roots = np.argwhere(np.sum(A, axis=0) == 0).reshape((-1,))
    
    # do a few sanity checks. Need a single root
    if len(roots) > 1:
        raise ValueError("Tree has multiple roots")

    # if root is unset, set    
    if root == None:
        return roots[0]

    # if root is set but doesn't match
    if roots[0] != root:
        raise ValueError("Trees have different roots")
    
    # else it's ok
    return roots[0]


def die_with_usage():
    print """
    SuMoTED
    =======

    Subtree Moving Tree Edit Distance code. Computes the pairwise distance
    between trees in a directory using the SuMoTED algorithm. Usage:

      $ python SuMoTED.py directory

    here directory is a path to a directory of trees. Trees should be specified
    in a file listing (parent, child) relationships. Trivial example found in /data/toy:

                   A                       A 
                 /   \\                   /   \\
        T1  =   B     C       T2  =     D     C
                |                             |
                D                             B

   then /data/toy/ should contain two files named T1 and T2, with content:
       T1            T2
       A, B          A, D
       A, C          A, C
       B, D          C, B

    (as it does). We can convert T1 to T2 with a local upward move of D to 
    be a child of A (cost 1), then a local downward move of B to be a child of 
    C (cost 1). Therefore, the un-normalised distance between T1 and T2 is 2, which
    is actually the maximal distance between trees with these nodes, so they
    have a normalised similarity = 0.0. This can be verified by running:

      $ python SuMoTED.py ./data/toy

    """
    sys.exit()


# Main -----------------------------------------------------------------------
if __name__ == "__main__":

    import sys, os
    if len(sys.argv) - 1 != 1:
        die_with_usage()

    # find all the trees in the input directory
    print ''
    print '  finding all trees in the directory...'
    tree_dir = sys.argv[1]
    tree_files, n_trees = find_all_trees(tree_dir)

    # read in trees and store as adjacency matrices. This also
    # pulls out the root index, stores the label names and 
    # computes the Jaccard similarity between the trees
    print '  building adjacency matrices...'
    trees, root_index, Jaccard, labels = tree_files_to_adjacency_matrices(tree_files)

    # make bush (for normalisation)
    print '  making bush (for normalisation)...'
    tree_size = trees[0].shape[0]
    bush = make_bush(tree_size, root_index)
    print ''

    # now compute similarities
    Dmat = np.zeros((n_trees, n_trees), dtype=int)
    Dmat_norm = np.ones((n_trees, n_trees))

    # initialise intersection
    for it1 in range(n_trees):
        # compute bush distance here. We don't need
        # the transitive closures
        t1_bush_dist = distance(trees[it1], bush)

        for it2 in range(it1 + 1, n_trees):
            print "  computing distance between", tree_files[it1], 'and', tree_files[it2]

            # as above
            t2_bush_dist = distance(trees[it2], bush)

            # main distance
            Dmat[it1, it2] = distance(trees[it1], trees[it2])

            # check if norm can be computed
            norm = t1_bush_dist + t2_bush_dist
            if norm > 0:
                Dmat_norm[it1, it2] = 1 - Dmat[it1, it2] / float(norm)
            else:
                Dmat_norm[it1, it2] = 1.0

    print ''
    print "  Label similarities:"
    print Jaccard

    print ''
    print "  Number of tree edits:"
    print Dmat.T

    print ''
    print "  Normalised tree similarities:"
    print Dmat_norm
    print ''
