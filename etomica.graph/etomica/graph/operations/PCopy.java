package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Graph;

/*
 * PCopy (Positive Copy): GraphSet --> GraphSet
 *
 * Given a graph set S, returns a graph set S' such that every graph g in S has an exact
 * copy g' in S' and coefficient(g') = coefficient(g).
 */
public class PCopy implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {

    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      result.add(g.copy());
    }
    return result;
  }
}