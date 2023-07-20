#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys

import pandas as pd
import numpy as np
import scipy.linalg as slin
import scipy.optimize as sopt
from scipy.special import expit as sigmoid
import csv
import sys
import igraph as ig
import os
import timeit
import random


# In[2]:


import glob


# In[3]:


def notears_linear(X,Z,lambda1,lambda2,loss_type, max_iter=100, h_tol=1e-8, rho_max=1e+16, w_threshold=0.3):

    ## divided W into covariates and network
    full_data=np.concatenate([X,Z],axis=1)
    def _divideParameters(W):
        W_net=W[0:d]
        W_cov=W[d:]
        return W_net,W_cov

    ## Nothing needs to be changed here
    def _loss(W):
        
        M = np.matmul(full_data,W)
        if loss_type == 'l2':
            R = X - M
            loss = 0.5 / full_data.shape[0] * (R ** 2).sum()
            G_loss = - 1.0 / full_data.shape[0] * np.matmul(full_data.T,R)
        else:
            raise ValueError('unknown loss type')
        return loss, G_loss

    ## Nothing needs to be changed here
    def _h(W):
        
        W_net,W_cov=_divideParameters(W)
        E = slin.expm(W_net * W_net)  # (Zheng et al. 2018)
        h = np.trace(E) - d
        #     # A different formulation, slightly faster at the cost of numerical stability
        #     M = np.eye(d) + W * W / d  # (Yu et al. 2019)
        #     E = np.linalg.matrix_power(M, d - 1)
        #     h = (E.T * M).sum() - d
        G_h = E.T * W_net * 2
        G_h_expand=np.concatenate([G_h,np.zeros([m,d])],axis=0)
        return h, G_h_expand

    ## Nothing needs to be changed here except for the dimension
    def _adj(w):
        
        return (w[:((d+m) * d)] - w[((d+m) * d):]).reshape([(d+m), d])

    def _func(w):
        
        W = _adj(w)
        W_net,W_cov=_divideParameters(W)
        loss, G_loss = _loss(W)
        h, G_h = _h(W)
        obj = loss + 0.5 * rho * h * h + alpha * h + lambda1 * W_net.sum()+lambda2* W_cov.sum()
        G_smooth = G_loss + (rho * h + alpha) * G_h
        G_nonsmooth=np.concatenate([lambda1*np.ones([d,d]),lambda2*np.ones([m,d])])
        g_obj = np.concatenate((G_smooth + G_nonsmooth, - G_smooth + G_nonsmooth), axis=None)
        return obj, g_obj

    n, d = X.shape
    ## Added: m=number of covariates
    m=Z.shape[1]
    w_est, rho, alpha, h = np.zeros(2 * (d+m) * d), 1.0, 0.0, np.inf  # double w_est into (w_pos, w_neg)
    bnds = [(0, 0) if i == j else (0, None) for _ in range(2) for i in range(d+m) for j in range(d)]
    if loss_type == 'l2':
        full_data = full_data - np.mean(full_data, axis=0, keepdims=True)
    mycount=0
    for _ in range(max_iter):
        mycount=mycount+1
        print(mycount)
        w_new, h_new = None, None
        while rho < rho_max:
            sol = sopt.minimize(_func, w_est, method='L-BFGS-B', jac=True, bounds=bnds)
            w_new = sol.x
            h_new, _ = _h(_adj(w_new))
            if h_new > 0.25 * h:
                rho *= 10
            else:
                break
        w_est, h = w_new, h_new
        alpha += rho * h
        if h <= h_tol or rho >= rho_max:
            break
    W_est = _adj(w_est)
    #W_est[np.abs(W_est) < w_threshold] = 0
    return W_est




# In[4]:


# A class to represent a graph object
class Graph:
    # Constructor
    def __init__(self, edges, n):
 
        # A list of lists to represent an adjacency list
        self.adjList = [[] for _ in range(n)]
 
        # add edges to the directed graph
        for (src, dest) in edges:
            self.adjList[src].append(dest)
 
 
# Perform DFS on the graph and set the departure time of all vertices of the graph
def DFS(graph, v, discovered, departure, time):
 
    # mark the current node as discovered
    discovered[v] = True
 
    # do for every edge (v, u)
    for u in graph.adjList[v]:
        # if `u` is not yet discovered
        if not discovered[u]:
            time = DFS(graph, u, discovered, departure, time)
 
    # ready to backtrack
    # set departure time of vertex `v`
    departure[v] = time
    time = time + 1
 
    return time
 
 
# Returns true if the given directed graph is DAG
def isDAG(graph, n):
 
    # keep track of whether a vertex is discovered or not
    discovered = [False] * n
 
    # keep track of the departure time of a vertex in DFS
    departure = [None] * n
 
    time = 0
 
    # Perform DFS traversal from all undiscovered vertices
    # to visit all connected components of a graph
    for i in range(n):
        if not discovered[i]:
            time = DFS(graph, i, discovered, departure, time)
 
    # check if the given directed graph is DAG or not
    for u in range(n):
 
        # check if (u, v) forms a back-edge.
        for v in graph.adjList[u]:
 
            # If the departure time of vertex `v` is greater than equal
            # to the departure time of `u`, they form a back edge.
 
            # Note that `departure[u]` will be equal to `departure[v]`
            # only if `u = v`, i.e., vertex contain an edge to itself
            if departure[u] <= departure[v]:
                return False
 
    # no back edges
    return True
 
 
if __name__ == '__main__':
 
    # List of graph edges as per the above diagram
    edges = [(0, 1), (0, 3), (1, 2), (1, 3), (3, 2), (3, 4), (3, 0), (5, 6), (6, 3)]
 
    # total number of nodes in the graph (labelled from 0 to 6)
    n = 7
 
    # build a graph from the given edges
    graph = Graph(edges, n)
 
    # check if the given directed graph is DAG or not
    if isDAG(graph, n):
        print('The graph is a DAG')
    else:
        print('The graph is not a DAG')


# In[5]:


def dagThreshold(G,threshold):
    n=G.shape[0]
    edgelist=[]
    for i in range(0,n):
        for j in range(0,n):
            if abs(G[i,j])>threshold:
                edgelist.append([i,j])
    graph = Graph(edgelist, n)
    return isDAG(graph,n)

def findmin(G):
    mymin=0.01
    while True:
        if dagThreshold(G,mymin)==True:
            G[np.abs(G)<=mymin]=0
            return G,mymin
        else:
            mymin=mymin+0.01


    
  


# In[6]:


for file in os.listdir(path):
    if file.endswith(".csv"):
        filename=os.path.splitext(file)[0]
        filepath=os.path.join(path,file)
        dataset=pd.read_csv(filepath,",")
        ngene=dataset.shape[1]-2
        X_data=dataset.iloc[:,range(1,ngene+1)]
        Y_data=dataset.iloc[:,range(ngene+1,ngene+2)]
        genelist=dataset.columns[1:ngene+1]
        X=X_data.to_numpy()
        Y=Y_data.to_numpy()
        for hyper in [0.1,0.05,0.01,0.005,0.001]:
            random.seed(10)
            W_est=notears_linear(X,Y,lambda1=hyper,lambda2=hyper,loss_type="l2",h_tol=1e-8,w_threshold=0.3)
            W_coef=W_est[range(0,ngene),:]
            W_cova=W_est[range(ngene,ngene+1),:]
            W_coef_copy=W_coef.copy()
            for i in range(0,W_coef_copy.shape[0]):
                W_coef_copy[i,i]=0
            finalfit=findmin(W_coef_copy)[0]
            NOA=sum(np.abs(finalfit.flatten()>0))
            output1=pd.DataFrame(finalfit, columns = dataset.columns[range(1,len(dataset.columns)-1)].tolist(),index=dataset.columns[range(1,len(dataset.columns)-1)].tolist())
            print(output1.shape)
            #top10gene=np.flip(np.argsort(W_cova_reshape))[0:10]+1
            output1.to_csv("./Result3/"+filename+"NOA"+str(NOA)+"hyper"+str(hyper)+".csv")
            W_cov_df=pd.DataFrame(W_cova,columns=dataset.columns[range(1,len(dataset.columns)-1)].tolist())
            W_cov_df.to_csv("./Result3/"+filename+str(hyper)+"Cov.csv")


# In[ ]:


finalfit=findmin(W_coef_copy)[0]
output1=pd.DataFrame(finalfit, columns = dataset.columns[range(1,len(dataset.columns)-1)].tolist(),index=dataset.columns[range(1,len(dataset.columns)-1)].tolist())
#W_cov_df.to_csv("chlodowncov.csv")
W_cov_df=pd.DataFrame(W_cova,columns=dataset.columns[range(1,len(dataset.columns)-1)].tolist())
W_cov_df.to_csv("chlodowncova.csv")


# In[ ]:




        
        


# In[5]:


#index=list(range(0,84))
#A=random.sample(range(84,702),84)
#index_combined=index+A


# In[6]:





# In[ ]:





# In[1]:





# In[9]:


#sum(np.ndarray.flatten(W_est)!=0)
W_coef[np.abs(W_coef) < 0.25]=0
W_cova[np.abs(W_cova) < 0.15]=0


# In[10]:


W_coef_df=pd.DataFrame(W_coef,columns=genelist,index=genelist)
W_cov_df=pd.DataFrame(W_cova,columns=genelist)


# In[11]:


W_coef_df.to_csv("chlodowncoef.csv")
W_cov_df.to_csv("chlodowncov.csv")


# In[11]:


class DAGLongestPath:
    """Calculate the longest path in a directed acyclic graph (DAG) in terms of node weights
    
    Use this class to get (one of) the paths with the largest sum of node weights
    in a directed acyclic graph (DAG). After constructing the empty object,
    use `add_node(label, weight)` and `add_edge(label1, label2)` to build the graph, 
    and then call `longest_path` to retrieve the path and the sum of the weights.
    This latter operation is destructive and will delete the graph.
    """
    
    def __init__(self):
        """Construct a new empty graph."""
        self.nodes = {}  # Dictionary {<label>:<weight>, ...}
        self.edges = {}  # Dictionary of sets dict{ <source_label>: set{<target_label>, ...}, ...}
        self.rev_edges = {}  # Dictionary of sets
        self.unseen_sources = set()  # Labels of all nodes not processed yet that have no incoming edges
        self.longest_in_weight = {}  # Dictionary {<label>:<weight>, ...}
        self.longest_in_route = {}   # Dictionary {<label>:[<label>, ...], ...}
        self.longest_route = None;   # The longest route (in weights) we have seen
        self.longest_route_weight = None;  # The larges weight we have seen
    
    def add_node(self, label, weight):
        """Add a node to a graph.
        
        # Arguments
            label: a scalar label for the node
            weight: a nonnegative number
        """
        if weight < 0: raise ValueError("weight cannot be negative")
        self.nodes[label] = weight
        self.edges[label] = set()
        self.rev_edges[label] = set()
        self.unseen_sources.add(label)
        
    def add_edge(self, source, target):
        """Add an edge to a graph.
        
        # Arguments
            source: the label of the source node; it should already exist in the graph
            target: the label of the target node; it should already exist in the graph
        """
        if source not in self.nodes: raise ValueError("source {} not a node".format(source))
        if target not in self.nodes: raise ValueError("target {} not a node".format(target))
        self.edges[source].add(target)
        self.rev_edges[target].add(source)
        self.unseen_sources.discard(target)
        
    def __del_edges_from(self, source):
        """Private method to delete all outgoing edges from a node."""
        targets = self.edges[source]
        self.edges[source] = set()
        for target in targets:
            self.rev_edges[target].discard(source)
            if len(self.rev_edges[target]) == 0: # no incoming edges
                self.unseen_sources.add(target)
                
    def __print(self):
        """Private method to print information about the graph."""
        print("Nodes, Edges")
        for id, w in self.nodes.items():
            print("  {}{} = {} -> {}".format(
                's' if id in self.unseen_sources else ' ', 
                id, 
                w,
                ",".join([str(x) for x in self.edges[id]])
            ))
        print("Rev-Edges")
        for id, source in self.rev_edges.items():
            print("  {} <- {}".format(id, ",".join([str(x) for x in source])))
        print("Longest in")
        for id, w in self.nodes.items():
            print("  {} : {} = {}".format(
                id,
                str(self.longest_in_weight.get(id, 0)),
                ",".join([str(x) for x in self.longest_in_route.get(id, [])])
            ))        
        print("")
        
        
    def longest_path(self):
        """Return the longest path in the graph in terms of the node weights.
        
        Warning: This operation is destructive and will delete the graph.
        
        # Returns
            An array of the route (array of labels), and the sum of the weights along the route.
        """
        while len(self.unseen_sources) > 0:
            sourcenode = self.unseen_sources.pop()
            
            new_weight = self.longest_in_weight.get(sourcenode, 0) + self.nodes[sourcenode]
            new_route = self.longest_in_route.get(sourcenode, []) + [sourcenode]

            if len(self.edges[sourcenode]) == 0: # no outgoing edges; isolated node
                if self.longest_route is None or self.longest_route_weight < new_weight:
                    self.longest_route = new_route
                    self.longest_route_weight = new_weight
                continue
            
            # There are outgoing edges            
            for target in self.edges[sourcenode]:
                
                if self.longest_in_weight.get(target, 0) < new_weight:
                    self.longest_in_weight[target] = new_weight
                    self.longest_in_route[target] = new_route
                
            self.__del_edges_from(sourcenode)
                
        return (self.longest_route, self.longest_route_weight)


if __name__ == '__main__':

    dag = DAGLongestPath()
    for i in range(6): dag.add_node(i, 1)
    dag.add_edge(3, 5)
    dag.add_edge(2, 5)
    dag.add_edge(2, 4)
    dag.add_edge(1, 2)
    print(dag.longest_path())


# In[44]:


path="/Users/taoxu/Desktop/COVID19PAPER/Result/ResultProcessed/Coefficient"
for file in os.listdir(path):
    if file.endswith(".csv"):
        A=pd.read_csv("/Users/taoxu/Desktop/COVID19PAPER/Result/ResultProcessed/Coefficient/"+file,index_col=0)
        with open(file+'.txt', 'w') as f:
            dag=DAGLongestPath()
            for i in range(A.shape[0]):
                dag.add_node(i,1)
            for j in range(A.shape[0]):
                for k in range(A.shape[0]):
                    if abs(A.iloc[j,k])>0.001:
            
                        dag.add_edge(j,k)
            result=dag.longest_path()
            colnames=list(A.columns)
            geneindex=result[0]
            for l in geneindex:
                f.write(colnames[l]+" ")
            
        
        
        
    


# In[ ]:




