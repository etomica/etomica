package etomica.virial.cluster2.graph.algorithms;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.Nodes;

/**
 * 
 * This interface generalizes algorithms that decide whether a graph G has a
 * particular property P. The property checker algorithm only needs to look at
 * the graph nodes and vertices to decide whether G has property P.
 * 
 * An example of an implementing class is a class that checks if G is connected.
 * 
 */
public interface GraphProperty {

  public boolean check(Nodes nodes, Edges edges);
}
