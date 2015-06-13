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
	

if __name__=='__main__':
	lines=[l.strip() for l in open('data.txt').readlines()]
	kmers=[l.strip() for l in lines[1:]]
	k=int(lines[0])
	# text=lines[1].strip()
 eb42f19cb391b3ea5945d6c821a0604dddccfefd
	# graph = debruijn_kmers(lines)
	# keys = sorted(graph.keys())
	# for k in keys:
	# 	print '%s -> %s'%(k, ",".join(sorted(graph[k])))

	# graph = read_graph(lines)
	# euler_path = find_euler_path(graph)
	# print "->".join(euler_path)

	
	# print "->".join(euler_path)
	# print euler_path
	# print string_spell(find_euler_path(debruijn_kmers(kmers)))
	print universal_string(9)
