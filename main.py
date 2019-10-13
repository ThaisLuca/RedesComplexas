
from __future__ import division
from graph_tool.all import *
import os, sys
from pylab import *
import numpy as np

LES_MISERABLES="/lesmis/lesmis.gml"
POLITICAL_BLOGS="/polblogs/polblogs.gml"

def new_graph(file):
	path = os.getcwd() + "/datasets" + file
	return load_graph(path)

def draw_graph(g, name):
	name = name +".png"
	graph_draw(g, vertex_text=g.vertex_index, vertex_font_size=14, output_size=(800, 800), output=name)	
	
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
	
	print("-- Centralidade de Grau: ")
	cent = central_measure(degrees, g.num_vertices())
	print("   Maximo: %f" % max(cent))
	print("   Minimo: %f" % min(cent))
	print("   Centralidade media (media = %.3f, desvio padrao = %.3f)" % (np.mean(cent), np.std(cent)))
	print("   Mediana: %f" % median(cent))
	print "\n"
	#plot_vertex_central(cent, "lemis")
	
def betweeness_measures(g):
	print("-- Betweeness: ")
	bv, be = betweenness(g)
	print("   Maximo: %f" % bv.a.max())
	print("   Minimo: %f" % bv.a.min())
	print("   Betweeness media (media = %.3f, desvio padrao = %.3f)" % (bv.a.mean(), bv.a.std()))
	print("   Mediana: %f" % median(bv.a))
	print "\n"
	#plot_betweness_hist(bv, "lemis")
	plot_graph_betweeness(g)
	
def closeness_measures(g):
	print("-- Closeness: ")
	c = closeness(g)
	print("   Maximo: %f" % c.a.max())
	print("   Minimo: %f" % c.a.min())
	print("   Betweeness media (media = %.3f, desvio padrao = %.3f)" % (c.a.mean(), c.a.std()))
	print("   Mediana: %f" % median(c.a))
	print "\n"
	
	
def plot_graph_betweeness(g):
	g = GraphView(g, vfilt=label_largest_component(g))
	vp, ep = betweenness(g)
	graph_draw(g, vertex_fill_color=vp, vertex_size=prop_to_size(vp, mi=5, ma=15), vcmap=matplotlib.cm.gist_heat, vorder=vp, output="lesmis_betweenness.pdf")
	
def central_measure(degrees, n_vertices):
	c = []
	for d in degrees:
		c.append(d/(n_vertices-1))
	return c

def main():
	g = new_graph(LES_MISERABLES)
	print "Graph: LES_MISERABLES\n"
	print("   Numero de vertices: %d" % g.num_vertices())
	print("   Numero de arestas: %d" % g.num_edges())
	print("\n")
	
	#degree_measures(g)
	#betweeness_measures(g)
	
	closeness_measures(g)
	
main()

def plot_vertex_hist(g, name):
	in_hist = vertex_hist(g, "total")
	y = in_hist[0]
	err = sqrt(in_hist[0])
	err[err >= y] = y[err >= y] - 1e-2

	figure(figsize=(6,4))
	errorbar(in_hist[1][:-1], in_hist[0], fmt="o", yerr=err, label="in")
	gca().set_yscale("log")
	gca().set_xscale("log")
	gca().set_ylim(1e-1, 1e5)
	gca().set_xlim(0.8, 1e3)
	subplots_adjust(left=0.2, bottom=0.2)
	xlabel("$k$")
	ylabel("$NP(k)$")
	tight_layout()
	name = name + "-deg-dist.pdf"
	savefig(name)

