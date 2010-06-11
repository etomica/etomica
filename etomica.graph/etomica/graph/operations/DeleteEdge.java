package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Edge;
import etomica.graph.model.Graph;

/**
 * Deletes all edges of a given color.  This is useful for e=f+1.  After using
 * substitution with an explicit 1-bond, you can use this class to delete the
 * resulting 1-bonds.
 * 
 * @author Andrew Schultz
 */
public class DeleteEdge implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {

    assert (params instanceof DeleteEdgeParameters);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      Set<Graph> newSet = apply(g, (DeleteEdgeParameters) params);
      if (newSet != null) {
        result.addAll(newSet);
      }
    }
    Unary isoFree = new IsoFree();
    return isoFree.apply(result, params);
  }

  public Set<Graph> apply(Graph g, DeleteEdgeParameters params) {

    Set<Graph> result = new HashSet<Graph>();
    Graph newGraph = g.copy();
    for (Edge edge : g.edges()) {
      if (edge.getColor() == params.color()) {
        byte fromNode = g.getFromNode(edge.getId());
        byte toNode = g.getToNode(edge.getId());
        newGraph.deleteEdge(fromNode, toNode);
      }
    }
    result.add(newGraph);
    return result;
  }
}