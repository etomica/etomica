package etomica.graph.property;


import etomica.graph.model.Graph;
import etomica.graph.traversal.Biconnected;
import etomica.graph.traversal.Traversal;

public class IsBiconnected implements Property {

  private Traversal bcTraversal;

  public IsBiconnected() {

    this.bcTraversal = new Biconnected();
  }

  public boolean check(Graph graph) {

    // by definition, a null graph is not biconnected
    if (graph.getNodeCount() == 0) {
      return false;
    }
    // by definition, a singleton graph and a two node graph are biconnected
    if (graph.getNodeCount() == 1) {
      return true;
    }
    // invariant: a biconnected graph has at least N edges
    if (graph.getEdgeCount() < graph.getNodeCount()) {
      return false;
    }
    // traverse the graph starting at nodeID and return true IFF all nodes in the graph
    // are traversed in the same connected component
    return bcTraversal.traverseComponent((byte) 0, graph, null);
  }
}