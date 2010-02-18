package etomica.graph.operations;

import java.util.Set;

import etomica.graph.model.Graph;

/*
 * Sum (GraphSet Sum): GraphSet x GraphSet --> GraphSet
 *
 * Given two graph sets S1 and S2, Sum(S1,S2) = IsoFree(Union(S1,S2)).
 */
public class Sum implements Binary {

  public Set<Graph> apply(Set<Graph> left, Set<Graph> right, Parameters params) {

    Binary union = new Union();
    Unary isoFree = new IsoFree();
    return isoFree.apply(union.apply(left, right, params), params);
  }
}