
from __future__ import division
from graph_tool.all import *
from collections import Counter
import os, sys
from pylab import *
import numpy as np
import matplotlib.pyplot as plt

LES_MISERABLES="/lesmis/lesmis.gml"

POLITICAL_BLOGS="/polblogs/polblogs.gml"
POWER="/power/power.gml"
COND_MAT_2005="/cond-mat-2005/cond-mat-2005.gml"
ASTRO_PH="/astro-ph/astro-ph.gml"

def plot(data,filename,degreetype):
    """ Plot Distribution """
    plt.plot(range(len(data)),data,'bo')
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel('Freq')
    plt.xlabel('Degree')
    plt.savefig(filename + '_' + degreetype + '_distribution.pdf', bbox_inches="tight")
    plt.clf()

    """ Plot CDF """
    s = float(data.sum())
    cdf = data.cumsum(0)/s
#    plt.plot(range(len(cdf)),cdf,'bo')
#    plt.xscale('log')
#    plt.ylim([0,1])
#    plt.ylabel('CDF')
#    plt.xlabel('Degree')
#    plt.savefig(filename + '_' + degreetype + '_cdf.pdf')
#    plt.clf()

    """ Plot CCDF """
    ccdf = 1-cdf
    plt.plot(range(len(ccdf)),ccdf,'bo')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([0,1])
    plt.ylabel('CCDF')
    plt.xlabel('Degree')
    plt.savefig(filename + '_' + degreetype + '_ccdf.pdf', bbox_inches="tight")
    plt.clf()

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
		if g.is_directed():
			d += a.in_degree()
		degrees.append(d)
	return degrees

def degree_measures(g):
	print("-- Grau: ")
	degrees = get_vertices_degree(g)
	print("   Maximo: %d" % max(degrees))
	print("   Minimo: %d" % min(degrees))
	print("   Grau medio (media = %.3f, desvio padrao = %.3f)" % vertex_average(g, "total"))
	print("   Mediana: %d" % median(degrees))
	print "\n"
	#plot_vertex_hist(g, "lemis")

	#indegree_distribution = np.bincount(degrees)

	#plot(indegree_distribution, "teste", 'totaldegree')
	
	print("-- Centralidade de Grau: ")
	cent = central_measure(degrees, g.num_vertices())
	print("   Maximo: %f" % max(cent))
	print("   Minimo: %f" % min(cent))
	print("   Centralidade media (media = %.3f, desvio padrao = %.3f)" % (np.mean(cent), np.std(cent)))
	print("   Mediana: %f" % median(cent))
	print "\n"
		
#Consertar
def get_components(g):
	comp, hist = label_components(g, attractors=True)
	return Counter(comp)
	

def central_measure(degrees, n_vertices):
	c = []
	for d in degrees:
		c.append(d/(10-1))
	return c

def print_stats(x):
	print("   Maximo: %f" % x.a.max())
	print("   Minimo: %f" % x.a.min())
	print("   Media (media = %.3f, desvio padrao = %.3f)" % (x.a.mean(), x.a.std()))
	print("   Mediana: %f" % median(x.a))
	print "\n"

def main():
	g = new_graph(ASTRO_PH)
	print(g.is_directed())
	print("Graph: ASTRO_PH\n")
	print("   Numero de vertices: %d" % g.num_vertices())
	print("   Numero de arestas: %d" % g.num_edges())
	print("   Densidade: %3f" % (2*g.num_edges()/(g.num_vertices()*g.num_vertices()-g.num_vertices())))
	dist, ends = pseudo_diameter(g)
	print("   Diametro: %3f " % dist)

	print("\n")
	
	degree_measures(g)
	
	print("-- Betweeness: ")
	bv, be = betweenness(g)
	print_stats(bv)
	
	print("-- Closeness: ")
	print_stats(closeness(g))
	
	print("-- Katz:")
	print_stats(katz(g))
	
	print("-- Autovetor: ")
	ee, x = eigenvector(g)
	print_stats(x)
	
	print("-- Clusterizacao Local: ")
	print_stats(local_clustering(g))
	
	print("-- PageRank: ")
	print_stats(pagerank(g))
	
	print("-- Clusterizacao Global: ")
	c = global_clustering(g)
	print("   Coeficiente global de clusterizacao: %f" % c[0])
	print("   Desvio Padrao: %f" % c[1])
	print "\n"
	
main()

def plot_graph_metrics(g, x, name, metric):
	name = name + "_" + metric + ".pdf"
	g = GraphView(g, vfilt=label_largest_component(g))
	graph_draw(g, vertex_fill_color=x, vertex_size=prop_to_size(x, mi=5, ma=15), vcmap=matplotlib.cm.gist_heat, vorder=x, output=name)

