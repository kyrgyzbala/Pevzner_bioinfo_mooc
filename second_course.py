import sys

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


def universal_string():
	alp = ['0','1']
	all_kmers = [a+b+c+d for a in alp for b in alp for c in alp for d in alp]
	graph = overlap_graph_map(all_kmers)
	path = get_span(graph)
	print string_spell(path)	

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

def is_universal(text):
	alp = ['0','1']
	all_kmers = [a+b+c for a in alp for b in alp for c in alp]
	for kmer in all_kmers:
		if kmer not in text:
			return False
	return True


if __name__=='__main__':
	# lines=[l.strip() for l in open('data.txt').readlines()]
	# # k=int(lines[0])
	# # text=lines[1].strip()
	# graph = debruijn_kmers(lines)
	# keys = sorted(graph.keys())
	# for k in keys:
	# 	print '%s -> %s'%(k, ",".join(sorted(graph[k])))
	nums = ['1110001011',
			'0111010010',
			'0011101000',
			'1001101100',
			'0011100100',
			'0111010001']
	for n in nums:
		print n, is_universal(n)
	
