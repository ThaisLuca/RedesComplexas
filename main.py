
from __future__ import division
from graph_tool.all import *
from collections import Counter
import os, sys
from pylab import *
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import math
import random as rd


LES_MISERABLES="/lesmis/lesmis.gml"

POLITICAL_BLOGS="/polblogs/polblogs.gml"
POWER="/power/power.gml"
ASTRO_PH="/astro-ph/astro-ph.gml"
INTERNET="/internet/internet.gml"

def plot_distribution(graph, all_measure, xlabel, filename, metric='degree'):
    freq = freq_relative(graph, all_measure, metric)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('CDF')
    plt.xlabel(xlabel)
    plt.plot(range(len(freq)), freq, 'o', clip_on=False)
    plt.savefig('graficos/'+filename+'_cdf.jpg')
    plt.clf()


def plot_ccdf(graph, all_measure, xlabel, filename, metric='degree'):
    _ccdf = ccdf(graph, all_measure, metric)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('CCDF')
    plt.xlabel(xlabel
    plt.plot(range(len(_ccdf)), _ccdf, 'o', clip_on=False)
    plt.savefig('graficos/'+filename+'_ccdf.jpg')
    plt.clf()

    
def freq_relative(graph, all_measure, metric='Degree'):
    if metric == 'Degree':
        degree_distribution = np.bincount(list(all_measure))
        return degree_distribution/len(graph.get_vertices())
    elif metric == 'Distance':
        distance_distribution = np.bincount(list(all_measure))
        return distance_distribution/comb(get_num_vertex(graph), 2)
    else:
        all_measure = np.array(all_measure)
        all_sum = float(all_measure.sum())
        return all_measure.cumsum(0)/all_sum  

def ccdf(graph, all_measure, metric='degree'):
    return 1 - freq_relative(graph, all_measure, metric)
	
def new_graph(file):
	path = os.getcwd() + "/datasets" + file
	return load_graph(path)

def draw_graph(g, name):
	name = name +".pdf"
	graph_draw(g, vertex_text=g.vertex_index, vertex_font_size=10, output_size=(800, 800), output=name)	
	
def get_vertices_degree(g):
	degrees = []
	for v in g.get_vertices():
		d = 0
		a = g.vertex(v)
		d = a.out_degree()
		#d = a.in_degree()
		degrees.append(d)
	return degrees

def degree_measures(g):
	print("-- Grau: ")
	degrees = get_vertices_degree(g)
	print("   Maximo: %d" % max(degrees))
	print("   Minimo: %d" % min(degrees))
	print("   Grau medio (media = %.3f, desvio padrao = %.3f)" % vertex_average(g, "total"))
	print("   Mediana: %d" % median(degrees))
	print("\n")
	return degrees

def get_distances(g):
	
	source = np.random.randint(low=0,high=len(g.get_vertices()), size=100)
 	dist = []
 	for source in source:
  		target_vertex = np.random.randint(low=0,high=len(g.get_vertices()), size=100)
  		for target in target_vertex:
  			a = shortest_distance(g, source, target=target)
  			if not(a == 0 or isinf(a)):
   				dist.append(a)
			
	print("-- Distancia: ")
	print("   Maximo: %f" % max(dist))
	print("   Minimo: %f" % min(dist))
	print("   Distancia media (media = %.3f, desvio padrao = %.3f)" % (np.mean(dist), np.std(dist)))
	print("   Mediana: %f" % median(dist))
	print("\n")
	return dist
		
def get_components(g):
	v = []
	comp, hist = label_components(g)
	comp = Counter(comp.a)
	for key, value in comp.items():
		 v.append(value)
	
	print("-- Componentes: ")
	print("   Maximo: %f" % max(v))
	print("   Minimo: %f" % min(v))
	print("   Distancia media (media = %.3f, desvio padrao = %.3f)" % (np.mean(v), np.std(v)))
	print("   Mediana: %f" % median(v))
	print("\n")

def central_measure(degrees, n_vertices):
	c = []
	for d in degrees:
		c.append(d/(10-1))
	return c

def print_stats(x):
	print("   Maximo: %f" % max(x))
	print("   Minimo: %f" % min(x))
	print("   Media (media = %.3f, desvio padrao = %.3f)" % (np.mean(x), np.std(x)))
	print("   Mediana: %f" % median(x))
	print("\n")

def main():
	g = new_graph(INTERNET)
	print(g.is_directed())
	print("Graph: INTERNET\n")
	print("   Numero de vertices: %d" % g.num_vertices())
	print("   Numero de arestas: %d" % g.num_edges())
	print("   Densidade: %3f" % (2*g.num_edges()/(g.num_vertices()*g.num_vertices()-g.num_vertices())))
	dist, ends = pseudo_diameter(g)
	print("   Diametro: %3f " % dist)
	
	all_dist = get_distances(g)
	plot_distribution(g, all_dist, "Distancia", "political_blogs_distancias", 'Distance')
	plot_ccdf(g, all_dist, "Distancia", "political_blogs_distancias", 'Distance')	

	print("\n")
	
	plot_distribution(g, degree_measures(g), "Grau", "astro_graus")
	plot_ccdf(g, degree_measures(g), "Grau", "political_blogs_entrada")
	degree_measures(g)
	
	print("-- Betweeness: ")
	bv, be = betweenness(g)
	all_ = []
	for v in bv:
		if isnan(v): all_.append(0)
		else: all_.append(v)
	print_stats(all_)
	plot_distribution(g, all_, "Betweeness", "astro_betweeness", metric='betweeness')
	plot_ccdf(g, all_, "Betweeness", "internet_betweeness", metric='betweeness')
	return
	
	print("-- Closeness: ")
	c = closeness(g)
	all_c = []
	for v in c:
		if isnan(v): all_c.append(0)
		else: all_c.append(v)
	print_stats(all_c)
	plot_distribution(g, all_c, "Closeness", "astro_closeness", metric='closeness')
	plot_ccdf(g, all_c, "Closeness", "astro_closeness", metric='closeness')
	
	print("-- Katz:")
	k = katz(g)
	all_k = []
	for v in k:
		if isnan(v): all_k.append(0)
		else: all_k.append(v)
	print_stats(all_k)
	plot_distribution(g, all_k, "Katz", "astro_katz", metric='katz')
	plot_ccdf(g, all_k, "Katz", "astro_katz", metric='katz')
	
	
	print("-- Autovetor: ")
	ee, x = eigenvector(g)
	all_x = []
	for v in x:
		if isnan(v): all_x.append(0)
		else: all_x.append(v)
	print_stats(all_x)
	plot_distribution(g, all_x, "Centralidade de Autovetor", "astro_autovetor", metric='autovetor')
	plot_ccdf(g, all_x, "Autovetor", "astro_autovetor", metric='autovetor')
	
	print("-- Clusterizacao Local: ")
	all_cl = []
	cl = local_clustering(g)
	for v in cl:
		if isnan(v): all_cl.append(0)
		else: all_cl.append(v)
	print_stats(all_cl)
	plot_distribution(g, all_cl, "Clusterizacao Local", "astro_clus_local", metric='cluster')
	plot_ccdf(g, all_cl, "Clusterizacao Local", "astro_clus_local", metric='cluster')
	
	print("-- Clusterizacao Global: ")
	c = global_clustering(g)
	print("   Coeficiente global de clusterizacao: %f" % c[0])
	print("   Desvio Padrao: %f" % c[1])
	print("\n")

def plot_graph_metrics(g, x, name, metric):
	name = name + "_" + metric + ".pdf"
	g = GraphView(g, vfilt=label_largest_component(g))
	graph_draw(g, vertex_fill_color=x, vertex_size=prop_to_size(x, mi=5, ma=15), vcmap=matplotlib.cm.gist_heat, vorder=x, output=name)

