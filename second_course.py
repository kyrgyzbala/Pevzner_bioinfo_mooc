import sys
import numpy as np


def string_composition(k, text):
	output=[]
	for i in range(0,len(text)-k+1):
		output.append(text[i:i+k])
	return output


def string_spell(kmers):
	return kmers[0]+"".join(kmer[-1] for kmer in kmers[1:])


def overlap_graph(kmers):
	output=[]
	for kmer1 in kmers:
		for kmer2 in kmers:
			if kmer1!=kmer2 and kmer1[1:]==kmer2[:-1]:
				output.append((kmer1, kmer2))
	return output


def overlap_graph_map(kmers):
	graph = overlap_graph(kmers)
	graph_map = {}
	for nodes in graph:
		if nodes[0] in graph_map:
			graph_map[nodes[0]].append(nodes[1])
		else:
			graph_map[nodes[0]] = [nodes[1]]
	return graph_map


def get_span(graph):
	path=[]
	path_len = 0

	def span(start, graph, path):
		if len(path)==16:			
			return path
		else:
			for node in graph[start]:
				if node not in path:
					return span(node, graph, path+[node])

	for start in graph.keys():		
		init_path = [start]
		peth_len = 1
		path = span(start, graph, init_path)
		if path:
			return path


def debruijn(k, text):
	kmers = string_composition(k-1, text)
	graph={}
	for i in range(len(kmers)-1):
		cur_kmer = kmers[i]
		if cur_kmer in graph:
			graph[cur_kmer] += [kmers[i+1]]
		else:
			graph[cur_kmer] = [kmers[i+1]]

	return graph


def debruijn_kmers(kmers):
	graph = {}
	for kmer in kmers:
		k = kmer[:-1]
		v = kmer[1:]
		if k in graph:
			graph[k] += [v]
		else:
			graph[k] = [v]
	return graph


def debruijn_paired(kmers):
	graph = {}
	fmt="%s|%s"
	for kmer in kmers:
		[pair1, pair2] = kmer.split('|')
		cur_node  = fmt%(pair1[:-1], pair2[:-1])
		next_node = fmt%(pair1[1:],  pair2[1:])
		if cur_node in graph:
			graph[cur_node] += [new_node]
		else:
			graph[cur_node] = [next_node]
	return graph

def find_euler_cycle(graph):
	path = []
	# nodes = graph.keys()
	nodes = sorted(graph.keys())
	stack = [nodes[0]]
	visited_edges = []

	while stack:
		cur_node = stack[-1]
		new_edges = [v for v in graph[cur_node] if '%s-%s'%(cur_node, v) not in visited_edges]
		if new_edges:
			new_node = new_edges[0]
			visited_edges.append('%s-%s'%(cur_node, new_node))
			stack.append(new_node)
		else:			
			path.append(stack.pop())

	path.reverse()
	return path


def find_euler_path(graph):
	nodes = []
	for k,v in graph.items():
		nodes += [k] + v
	nodes = sorted(set(nodes))
	M = np.zeros((len(nodes), len(nodes)))
	
	for k in graph.keys():
		ind_y = nodes.index(k)
		for n in graph[k]:
			ind_x = nodes.index(n)
			M[ind_y, ind_x]=1
	
	unbalanced_indices = [i for i in range(len(nodes)) if np.sum(M[i,:])!=np.sum(M[:,i])]

	assert(len(unbalanced_indices) in [2,0])

	if len(unbalanced_indices)==2:
		ind1 = unbalanced_indices[0]
		ind2 = unbalanced_indices[1]

		first, second = nodes[ind1], nodes[ind2]
		if sum(M[:,ind1])>sum(M[ind1,:]):
			pivot = second
			if first in graph:
				graph[first].append(second)
			else:
				graph[first]=[second]
		elif sum(M[:,ind2])>sum(M[ind2,:]):
			pivot = first
			if second in graph:
				graph[second].append(first)
			else:
				graph[second]=[first]

		euler_cycle = find_euler_cycle(graph)
		euler_cycle = euler_cycle[:-1]
		pivot = euler_cycle.index(pivot)
		return euler_cycle[pivot:] + euler_cycle[:pivot]
	
	euler_cycle = find_euler_cycle(graph)
	euler_cycle = euler_cycle[:-1]
	return euler_cycle


def read_graph(lines):
	graph={}
	for l in lines:
		first = l.split('->')[0].strip()
		second = l.split('->')[1].strip()
		graph[first] = second.split(',')
	
	return graph


def universal_string(k):
	alp = ['0','1']
	kmers = []

	for i in range(k):
		if i==0:
			kmers = alp
		else:
			prev_list = [t for t in kmers if len(t)==i]
			new_list = [t+u for t in prev_list for u in alp]
			kmers = new_list
	
	return string_spell(find_euler_cycle(debruijn_kmers(kmers))[:-k+1])


def paired_composition(k, d, text):

	composition = []
	for i in range(len(text)-2*k-d+1):
		composition.append([text[i:i+k], text[i+k+d:i+2*k+d]])
	return composition


def string_spelled_by_gapped_patterns(patterns, k, d):
	first_patterns = [p.split('|')[0] for p in patterns]
	second_patterns = [p.split('|')[1] for p in patterns]
	prefix_string = string_spell(first_patterns)
	suffix_string = string_spell(second_patterns)
	for i in range(k+d+1, len(prefix_string)):
		if prefix_string[i]!=suffix_string[i-k-d]:
			return None
	
	return prefix_string+suffix_string[-k-d:]

# def span(start, graph, path):
# 		if len(path)==16:			
# 			return path
# 		else:
# 			for node in graph[start]:
# 				if node not in path:
# 					return span(node, graph, path+[node])


def maximal_non_branching_paths(graph):
	paths = []

	nodes = []
	for k,v in graph.items():
		nodes += [k] + v
	nodes = sorted(set(nodes))
	
	M = np.zeros((len(nodes), len(nodes)))
	
	for k in graph.keys():
		ind_y = nodes.index(k)
		for n in graph[k]:
			ind_x = nodes.index(n)
			M[ind_y, ind_x]=1
	
	for node in nodes:
		node_ind = nodes.index(node)
		node_in = np.sum(M[:, node_ind])
		node_out = np.sum(M[node_ind, :])
		if not node_in==node_out==1 and node_out>0:
			for adj_node in graph[node]:
				non_branching_path = [node, adj_node]
				adj_node_ind = nodes.index(adj_node)
				while np.sum(M[:, adj_node_ind])==np.sum(M[adj_node_ind, :])==1:
					adj_node = graph[adj_node][0]
					adj_node_ind = nodes.index(adj_node)
					non_branching_path.append(adj_node)
				paths.append(non_branching_path)

	visited = set(sum(paths, []))
	left_nodes = set(graph.keys())-visited
	while left_nodes:
		node = left_nodes.pop()
		path = [node]
		if node in graph:
			adj_node = graph[node][0]
			path.append(adj_node)
			while adj_node in graph and adj_node!=node:
				adj_node = graph[adj_node][0]
				path.append(adj_node)
		paths.append(path)
		visited.update(path)
		left_nodes = set(graph.keys())-visited
	return paths



if __name__=='__main__':
	lines=[l.strip() for l in open('data.txt').readlines()]
	kmers=[l.strip() for l in lines[1:]]

	# print string_spell(find_euler_path(debruijn_kmers(kmers)))

	# k = int(lines[0].split()[0])
	# d = int(lines[0].split()[1])
	# for i in range(10):
		# print find_euler_path(debruijn_paired(kmers))
	# print string_spelled_by_gapped_patterns(find_euler_path(debruijn_paired(kmers)), k, d)
	for path in maximal_non_branching_paths(read_graph(lines)):
		print " -> ".join(path)