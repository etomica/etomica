package etomica.virial.cluster2.graph.algorithms.impl;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.Nodes;
import etomica.virial.cluster2.graph.algorithms.GraphAlgorithmsFactory;
import etomica.virial.cluster2.graph.algorithms.GraphProperty;

/**
 * A graph which is connected if there exists a path from any point to any other
 * point in the graph. A graph that is not connected is said to be disconnected.
 * This definition means that the null graph and singleton graph are considered
 * connected.
 * 
 * @author Demian Lessa
 * 
 */
public class CheckConnected implements GraphProperty {

  @Override
  public boolean check(Nodes nodes, Edges edges) {

    if ((nodes == null) || (edges == null)) {
      return false;
    }
    // by definition, a null graph and a singleton graph are connected
    if (nodes.count() <= 1) {
      return true;
    }
    // invariant: a connected graph has at least N-1 edges
    if (edges.count() < (nodes.count() - 1)) {
      return false;
    }
    // invoke a BF traversal starting at the first node (nodeID = 0) and
    // return true IFF all nodes in the graph are traversed
    return GraphAlgorithmsFactory.getBFTAlgo().traverseComponent(0, nodes,
        edges, null);
  }
}