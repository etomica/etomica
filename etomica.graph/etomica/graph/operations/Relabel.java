package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Edge;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.Node;

public class Relabel implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {

    assert (params instanceof RelabelParameters);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      Graph newGraph = apply(g, (RelabelParameters) params);
      result.add(newGraph);
    }
    return result;
  }

  public Graph apply(Graph argument, RelabelParameters params) {

    Node[] nodes = new Node[argument.nodeCount()];
    for (byte nodeId = 0; nodeId < nodes.length; nodeId++) {
      nodes[nodeId] = argument.getNode(params.map(nodeId));
    }
    // copy the mapped nodes
    Graph result = GraphFactory.createGraph(nodes);
    // copy the coefficient
    result.coefficient().multiply(argument.coefficient());
    // copy the edges from and to mapped nodes
    for (Edge edge : argument.edges()) {
      Edge newEdge = result.putEdge(params.map(argument.getFromNode(edge.getId())), params.map(argument
          .getToNode(edge.getId())));
      newEdge.setColor(edge.getColor());
    }
    return result;
  }
}
