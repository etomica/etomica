package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Graph;

public class NCopy implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {

    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      Graph ncopy = g.copy();
      ncopy.coefficient().switchSign();
      result.add(ncopy);
    }
    return result;
  }
}