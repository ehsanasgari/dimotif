import re
import scipy.cluster.hierarchy as hac
import matplotlib.pyplot as plt


class HierarchicalClutering(object):
    '''
    classdocs
    '''


    def __init__(self, distance_matrix, labels_out):
        '''
        Constructor
        '''
        z = hac.complete(distance_matrix)
        hac.dendrogram(z,labels=labels_out)
        self.tree = hac.to_tree(z,False)
        self.nwk=self.get_newick(self.tree, "", self.tree.dist, labels_out)
        plt.show()


    def get_newick(self, node, newick, parentdist, leaf_names):
        '''
        :param node:
        :param newick:
        :param parentdist:
        :param leaf_names:
        :return: the Newick format based on the provided distance matrix
        '''
        if node.is_leaf():
            return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
        else:
            if len(newick) > 0:
                newick = "):%.2f%s" % (parentdist - node.dist, newick)
            else:
                newick = ");"
            newick = self.get_newick(node.get_left(), newick, node.dist, leaf_names)
            newick = self.get_newick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
            newick = "(%s" % (newick)
            return newick
