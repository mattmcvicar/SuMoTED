import sys
import numpy as np


def open_tree(path, genre_codes):

    # Make into matrix
    n = len(genre_codes)
    M = np.zeros((n, n), dtype=int)

    for line in open(path, 'r'):
        if line[0] == '#':
            continue
        if len(line.strip()) == 0:
            continue

        parent, child = line.split(',')
        parent, child = parent.strip(), child.strip()

        M[genre_codes[parent], genre_codes[child]] = 1

    return M


def is_tree(A, genres, filename):
    n = A.shape[0]
    if sorted(np.sum(A, axis=0)) == [0] + [1]*(n-1):
        return True
    else:
        # print orphans
        for orphan in np.argwhere(np.sum(A, axis=0) == 0).reshape((-1,)):
            print genres[orphan], 'is an orphan'

        # print multiple parents
        for multi_parent in np.argwhere(np.sum(A, axis=0) > 1).reshape((-1,)):
            parents = [genres[g] for g in np.argwhere(A[    :, multi_parent]).reshape((-1,))]
            print genres[multi_parent], 'has multiple parents:',', '.join(parents)

        raise ValueError(filename + " is not a tree")



def distance(T1, T2):

    # check shape
    assert T1.shape == T2.shape
    n = T1.shape[0]

    # for the paper I wanted to know the maximum depth
    # parents = get_parents(T1, n)
    # ancestors = get_ancestors(parents, n)
    # maxdepth = max([len(b) for a, b in ancestors.items()])
    # print maxdepth

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

    # normalisation
    return D, T1trans, T2trans


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
            print '    - warning: found two possible trees: choosing one arbitrarily'
        x = x[-1]
        #x = [a for a in ancestors[y] if depths[a] == d - 1][0]

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

    bush = np.zeros((n, n), dtype=int)
    bush[root_index, :] = 1
    bush[root_index, root_index] = 0
    return bush


def find_all_trees(input_dir):

    # find all tree files
    tree_files = os.listdir(input_dir)
    tree_files = [os.path.join(input_dir, t) for t in tree_files if t != '.DS_Store']
    n_trees = len(tree_files)

    return tree_files, n_trees


def tree_files_to_adjacency_matrices(tree_files):

    # first do one pass to get the complete set of labels
    file_genres = []
    for tree_file in tree_files:
        file_genre = set()
        for line in open(tree_file, 'r'):
            line = line.strip()
            parent, child = line.split(',')
            parent, child = parent.strip(), child.strip()

            file_genre.add(parent)
            file_genre.add(child)

        file_genres.append(file_genre)

    # compute the label jaccard
    n_trees = len(tree_files)
    Jaccard = np.zeros((n_trees, n_trees))
    for i in range(n_trees):
        for j in range(i+1, n_trees):
            intersection = file_genres[i].intersection(file_genres[j])
            union = file_genres[i].union(file_genres[j])
            Jaccard[i,j] = (len(intersection) / float(len(union)))

    # now get all genres
    genres = set()
    for f in file_genres:
        genres = genres.union(f)

    genres = sorted(list(genres))
    n_genres = len(genres)

    root_index = genres.index("ALL_MUSIC")
    # Now go though again, building adjacency matrices
    As = []
    for tree_file in tree_files:
        A = np.zeros((n_genres, n_genres), dtype=int)
        for line in open(tree_file, 'r'):
            parent, child = line.split(',')
            parent, child = parent.strip(), child.strip()

            if parent == child:
                continue

            A[genres.index(parent), genres.index(child)] = 1

        # Post-processing: any row which is empty needs to be filled in
        # as a parent of the root
        to_insert = np.argwhere(np.sum(A, axis=0) == 0).reshape((-1,))
        if len(to_insert) > 1:
            print '    - warning: adding ' + str(len(to_insert)) + ' nodes to ' + tree_file
        for root in to_insert:
            if root != root_index:
                A[root_index, root] = 1

        # check A is a tree (and report if not)
        is_tree(A, genres, tree_file)
        As.append(A)

    return As, root_index, Jaccard, genres


# Main -----------------------------------------------------------------------
if __name__ == "__main__":

    import sys, os
    if len(sys.argv) - 1 < 1:
        raise ValueError("You must specify the path to a directory of trees")

    # find all the trees in the input directory
    print ''
    print '  finding all trees in the directory...'
    tree_dir = sys.argv[1]
    tree_files, n_trees = find_all_trees(tree_dir)

    # turn into adjacency matrices (this will also add nodes
    # if necessary)
    print '  building adjacency matrices...'
    trees, root_index, Jaccard, genres = tree_files_to_adjacency_matrices(
        tree_files)

    # make bush (for normalisation)
    print '  making bush (for normalisation)...'
    tree_size = trees[0].shape[0]
    bush = make_bush(tree_size, root_index)
    print ''

    # now compute similarities
    Dmat = np.zeros((n_trees, n_trees), dtype=int)
    Dmat_norm = np.ones((n_trees, n_trees))

    # initialise intersection
    I = np.ones((tree_size, tree_size), dtype=int)
    for it1 in range(n_trees):
        # compute bush distance here
        t1_bush_dist, _, _ = distance(trees[it1], bush)

        for it2 in range(it1 + 1, n_trees):

            # print progress
            print "  computing distance between", tree_files[it1], 'and', tree_files[it2]

            # as above
            t2_bush_dist, _, _ = distance(trees[it2], bush)

            # main distance
            Dmat[it1, it2], transT1, transT2 = distance(
                trees[it1], trees[it2])

            I = np.logical_and(I, transT1, dtype=int)
            I = np.logical_and(I, transT2, dtype=int)

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
