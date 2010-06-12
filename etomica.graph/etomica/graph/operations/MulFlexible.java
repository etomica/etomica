package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.Node;

/**
 * Perform multiplication for flexible molecules.  The result of multiplication
 * for each pair of diagrams is a diagram composed by both original diagrams
 * without any connection between them (the new diagram is disconnected).
 *
 * @author Andrew Schultz
 */
public class MulFlexible implements Binary {

  public Set<Graph> apply(Set<Graph> argument, Set<Graph> argument2, Parameters params) {

    assert (params == null);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      for (Graph g2 : argument2) {
        Graph newGraph = apply(g, g2);
        result.add(newGraph);
      }
    }
    return result;
  }

  public Graph apply(Graph g1, Graph g2) {
    Graph result;
    if (g1.nodeCount() == 1) {
      result = g2.copy();
      result.coefficient().multiply(g1.coefficient());
    }
    else if (g2.nodeCount() == 1) {
      result = g1.copy();
      result.coefficient().multiply(g2.coefficient());
    }
    else {
      byte nodes1 = g1.nodeCount();
      result = GraphFactory.createGraph((byte)(nodes1 + g2.nodeCount()));
      // add edges from g1
      for (Node node1 : g1.nodes()) {
        result.getNode(node1.getId()).setType(node1.getType());
        result.getNode(node1.getId()).setColor(node1.getColor());
        for (Node node2 : g1.nodes()) {
          if (node2.getId() <= node1.getId() || !g1.hasEdge(node1.getId(), node2.getId())) continue;
          result.putEdge(node1.getId(), node2.getId());
          result.getEdge(node1.getId(), node2.getId()).setColor(g1.getEdge(node1.getId(), node2.getId()).getColor());
        }
      }
      // now add edges from g2
      for (Node node1 : g2.nodes()) {
        result.getNode((byte)(node1.getId()+nodes1)).setType(node1.getType());
        result.getNode((byte)(node1.getId()+nodes1)).setColor(node1.getColor());
        for (Node node2 : g2.nodes()) {
          if (node2.getId() <= node1.getId() || !g2.hasEdge(node1.getId(), node2.getId())) continue;
          result.putEdge((byte)(node1.getId()+nodes1), (byte)(node2.getId()+nodes1));
          result.getEdge((byte)(node1.getId()+nodes1), (byte)(node2.getId()+nodes1)).setColor(g2.getEdge(node1.getId(), node2.getId()).getColor());
        }
      }
      result.coefficient().multiply(g1.coefficient());
      result.coefficient().multiply(g2.coefficient());
    }
    result.setNumFactors(g1.factors().length);
    result.addFactors(g1.factors());
    result.addFactors(g2.factors());
    return result;
  }
}
