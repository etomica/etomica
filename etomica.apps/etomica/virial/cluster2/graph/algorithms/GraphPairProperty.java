package etomica.virial.cluster2.graph.algorithms;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.Nodes;

/**
 * This interface generalizes algorithms that decide whether graphs G and G'
 * have a particular property P. The property checker algorithm must compare
 * the two graphs to decide whether G and G' have property P.
 * 
 * An example of an implementing class is a class that checks if G and G' are
 * isomorphs.
 * 
 */
public interface GraphPairProperty {

  /*
   * Assumes the two graphs are defined over a common set of nodes.
   * Or, at least, that there is a bijection from the nodes of each
   * graph to the nodes passed to this property checker.
   * 
   */
  public boolean check(Nodes nodes, Edges edges1, Edges edges2);
}
