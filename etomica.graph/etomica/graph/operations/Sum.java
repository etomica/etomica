package etomica.graph.operations;

import java.util.Set;

import etomica.graph.model.Graph;

public class Sum implements Binary {

  public Set<Graph> apply(Set<Graph> left, Set<Graph> right, Parameters params) {

    Binary union = new Union();
    Unary isoFree = new IsoFree();
    return isoFree.apply(union.apply(left, right, params), params);
  }
}