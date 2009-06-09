package etomica.virial.cluster2.graph.algorithms.impl;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.Nodes;
import etomica.virial.cluster2.graph.algorithms.GraphPairProperty;

public class CheckIsomorphic implements GraphPairProperty {

  @Override
  public boolean check(Nodes nodes, Edges edges1, Edges edges2) {
    
    if (nodes.count() <= 1) {
      return true;
    }
    if (edges1.count() != edges2.count()) {
      return false;
    }
    // create a degrees list L(nodes, edges1)
    // for every node u, add d(u, edges1) to L
    // for every node u, delete d(u, edges2) from L
    // if L is not empty, the graphs are not isomorphic
    return false;
  }

}
