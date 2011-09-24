package etomica.graph.operations;

import etomica.graph.model.Graph;

/**
 * Interface for a simple operation that takes a Graph and returns a Graph.
 * 
 * @author Andrew Schultz
 */
public interface GraphOp {
  public Graph apply(Graph g);
}