# NAME: hmm_example.py

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
#matplotlib inline

"""
A Markov chain (model) describes a stochastic process where the assumed probability 
of future state(s) depends only on the current process state and not on any the states 
that preceded it (shocker).

Let's get into a simple example. Assume you want to model the future probability that 
you land in a CpG island given its current state. To do this we need to 
specify the state space, the initial probabilities, and the transition probabilities.

Imagine you have a DNA sequence. We define the state space as the four diffent bases A,T,C and G. 
We will set the initial probabilities to 0.25%, 0.25%, 0.25% and 0.25% respectively.
"""

# create state space and initial state probabilities

states = ['a', 't', 'c' ,'g']
pi = [0.25, 0.25, 0.25, 0.25]
state_space = pd.Series(pi, index=states, name='states')

# create transition matrix
# equals transition probability matrix of changing states given a state
# matrix is size (M x M) where M is number of states

q_df = pd.DataFrame(columns=states, index=states)
q_df.loc[states[0]] = [0.180, 0.274, 0.426, 0.120]
q_df.loc[states[1]] = [0.171, 0.368, 0.274, 0.188]
q_df.loc[states[2]] = [0.161, 0.339, 0.375, 0.125]
q_df.loc[states[3]] = [0.079, 0.355, 0.384, 0.182]

q = q_df.values


"""
Now that we have the initial and transition probabilities setup we can create a 
Markov diagram using the Networkx package.

To do this requires a little bit of flexible thinking. Networkx creates Graphs 
that consist of nodes and edges. In our example the possible states are 
the nodes and the edges are the lines that connect the nodes. The transition 
probabilities are the weights. They represent the probability of transitioning 
to a state given the current state.

Something to note is networkx deals primarily with dictionary objects. With that 
said, we need to create a dictionary object that holds our edges and their weights.
"""

from pprint import pprint 

# create a function that maps transition probability dataframe 
# to markov edges and weights

def _get_markov_edges(Q):
    edges = {}
    for col in Q.columns:
        for idx in Q.index:
            edges[(idx,col)] = Q.loc[idx,col]
    return edges

edges_wts = _get_markov_edges(q_df)
"""
Now we can create the graph. To visualize a Markov model we need to 
use nx.MultiDiGraph(). A multidigraph is simply a directed graph which can have 
multiple arcs such that a single node can be both the origin and destination. 

In the following code, we create the graph object, add our nodes, edges, and 
labels, then draw a bad networkx plot while outputting our graph to a dot file. 
"""

# create graph object
G = nx.MultiDiGraph()

# nodes correspond to states
states = ['a', 't', 'c', 'g']
G.add_nodes_from(states)

# edges represent transition probabilities
for k, v in edges_wts.items():
    tmp_origin, tmp_destination = k[0], k[1]
    G.add_edge(tmp_origin, tmp_destination, weight=v, label=v) 

pos = nx.drawing.nx_pydot.graphviz_layout(G, prog='dot')
nx.draw_networkx(G, pos)
# In Windows: dot -Tps filename.dot -o outfile.ps


# create edge labels for jupyter plot but is not necessary
edge_labels = {(n1,n2):d['label'] for n1,n2,d in G.edges(data=True)}
nx.draw_networkx_edge_labels(G , pos, edge_labels=edge_labels)
nx.drawing.nx_pydot.write_dot(G, 'cpg_markov.dot')


"""
Lets us assume that you are traversing through the DNA sequence.
Consider a situation you encounter more CpG pairs along the DNA and you wanted to model 
the probability of you landing in a CpG island

In this situation the true state of the sequence is unknown, thus hidden from you. 
One way to model this is to assume hidden state. Let's walk through an example.

First we create our state space - CpG or Not-Cpg. We assume they are equiprobable.  
"""

# create state space and initial state probabilities

hidden_states = ['CpG', 'Not-Cpg']
pi = [0.5, 0.5]
state_space = pd.Series(pi, index=hidden_states, name='states')

# Next we create our transition matrix for the hidden states. 
# create hidden transition matrix
# a or alpha = transition probability matrix of changing states given a state
# matrix is size (M x M) where M is number of states

a_df = pd.DataFrame(columns=hidden_states, index=hidden_states)
a_df.loc[hidden_states[0]] = [0.7, 0.3]
a_df.loc[hidden_states[1]] = [0.4, 0.6]

a = a_df.values

"""
Now we create the emission or observation probability matrix. This matrix is size M x O where M is the number 
of hidden states and O is the number of possible observable states. 

The emission matrix tells us the probability that we are in one of the hidden 
states, given the current, observable state. 

Let's keep the same observable states from the previous example. We can be in
either A, T,C or G. For now we make our best guess to fill in 
the probabilities. 
"""

# create matrix of observation (emission) probabilities
# b or beta = observation probabilities given state
# matrix is size (M x O) where M is number of states 
# and O is number of different possible observations

observable_states = ['a', 't', 'c', 'g']

b_df = pd.DataFrame(columns=observable_states, index=hidden_states)
b_df.loc[hidden_states[0]] = [0.155, 0.341, 0.350, 0.154]
b_df.loc[hidden_states[1]] = [0.262, 0.246, 0.239, 0.253]

b = b_df.values


# Now we create the graph edges and the graph object. 
# create graph edges and weights

hide_edges_wts = _get_markov_edges(a_df)
#pprint(hide_edges_wts)

emit_edges_wts = _get_markov_edges(b_df)
#  pprint(emit_edges_wts)
print()

# create graph object
G = nx.MultiDiGraph()

# nodes correspond to states
G.add_nodes_from(hidden_states)


# edges represent hidden probabilities
for k, v in hide_edges_wts.items():
    tmp_origin, tmp_destination = k[0], k[1]
    G.add_edge(tmp_origin, tmp_destination, weight=v, label=v)

# edges represent emission probabilities
for k, v in emit_edges_wts.items():
    tmp_origin, tmp_destination = k[0], k[1]
    G.add_edge(tmp_origin, tmp_destination, weight=v, label=v)
      

pos = nx.drawing.nx_pydot.graphviz_layout(G, prog='neato')
nx.draw_networkx(G, pos)

# create edge labels for jupyter plot but is not necessary
emit_edge_labels = {(n1,n2):d['label'] for n1,n2,d in G.edges(data=True)}
nx.draw_networkx_edge_labels(G , pos, edge_labels=emit_edge_labels)
nx.drawing.nx_pydot.write_dot(G, 'pet_dog_hidden_markov.dot')
# In Windows: dot -Tps filename.dot -o outfile.ps

print("========================================================================")
print("                   CpG Island using HMMs")
print("========================================================================") 


def viterbi(pi, a, b, obs_seq):
    
    nStates = np.shape(b)[0]
    T = np.shape(obs_seq)[0]
    
    # init blank path
    path = np.zeros(T)
    # delta --> highest probability of any path that reaches state i
    delta = np.zeros((nStates, T))
    # phi --> argmax by time step for each state
    phi = np.zeros((nStates, T))
    
    # init delta and phi 
    delta[:, 0] = pi * b[:, obs_seq[0]]
    phi[:, 0] = 0

    # the forward algorithm extension
    for t in range(1, T):
        for s in range(nStates):
            delta[s, t] = np.max(delta[:, t-1] * a[:, s]) * b[s, obs_seq[t]] 
            phi[s, t] = np.argmax(delta[:, t-1] * a[:, s])
            
    
    # find optimal path
    print('-'*50)
    path[T-1] = np.argmax(delta[:, T-1])
    #p('init path\n    t={} path[{}-1]={}\n'.format(T-1, T, path[T-1])) #LPW
    for t in range(T-2, -1, -1): 
        path[t] = phi[int(path[t+1]), [t+1]]
        
    return path, delta, phi

"""
 
"""

# observation sequence of DNA
# observations are encoded numerically



obs_map = { 0:'a', 1:'t',2:'c',3:'g' }

filepath = "dna_seq.txt"

fp=open(filepath,'r+')
i=1
for line in fp.readlines():

    line = line.splitlines()

    char_array = []
    for entry in line:
        for c in entry:
            char_array.append(c)


    obs = np.array(char_array)

    inv_obs_map = dict((v,k) for k, v in obs_map.items())
    obs_seq = [inv_obs_map[v] for v in list(obs)]


    path, delta, phi = viterbi(pi, a, b, obs_seq)

    # Let's take a look at the result. 
    state_map = {0:'I', 1:'N'}
    state_path = [state_map[v] for v in path]
    
    print(' '.join(str(o) for o in obs))
    print(' '.join(str(p) for p in state_path))

print()    
 


"""
References

    https://en.wikipedia.org/wiki/Andrey_Markov
    https://www.britannica.com/biography/Andrey-Andreyevich-Markov
    https://www.reddit.com/r/explainlikeimfive/comments/vbxfk/eli5_brownian_motion_and_what_it_has_to_do_with/
    http://www.math.uah.edu/stat/markov/Introduction.html
    http://setosa.io/ev/markov-chains/
    http://www.cs.jhu.edu/~langmea/resources/lecture_notes/hidden_markov_models.pdf
    https://github.com/alexsosn/MarslandMLAlgo/blob/master/Ch16/HMM.py
    http://hmmlearn.readthedocs.io
    http://www.blackarbs.com/blog/introduction-hidden-markov-models-python-networkx-sklearn/2/9/2017
"""
