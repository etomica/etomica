package etomica.graph.operations;

import java.util.Set;

import etomica.graph.model.Graph;

public class Sub implements Binary {

  public Set<Graph> apply(Set<Graph> left, Set<Graph> right, Parameters params) {

    Binary union = new Union();
    Unary isoFree = new IsoFree();
    Unary negation = new NCopy();
    return isoFree.apply(union.apply(left, negation.apply(right, params), params), params);
  }
}