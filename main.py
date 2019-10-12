
from graph_tool.all import *
import os, sys

LES_MISERABLES="/lesmis/lesmis.gml"
POLITICAL_BLOGS="/polblogs/polblogs.gml"

def new_graph(file):
	path = os.getcwd() + "/datasets" + file
	return load_graph(path)

def draw_graph(g, name):
	name = name +".png"
	graph_draw(g, vertex_text=g.vertex_index, vertex_font_size=14, output_size=(800, 800), output=name)
	
def get_max_min_vertex_degree(g):
	degrees = []
	print(g.get_out_degrees(g.get_vertices()))
	total_degrees = g.get_total_degrees(g.get_vertices())
	return max(total_degrees), min(total_degrees)
		

def main():
	g = new_graph(LES_MISERABLES)
	print "Graph: LES_MISERABLES"
	print "-- Number of vertices: ", g.num_vertices()
	max_degree, min_degree = get_max_min_vertex_degree(g)
	print "-- Maximum vertice degree: ", max_degree
	print "-- Minimum vertice degree: ", min_degree
	print "-- Vertices average (degree, standart desviation) ", vertex_average(g, "total") 
main()



