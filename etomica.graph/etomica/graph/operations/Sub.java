package etomica.graph.operations;

import java.util.Set;

import etomica.graph.model.Graph;

/*
 * Sub (GraphSet Sum): GraphSet x GraphSet --> GraphSet
 *
 * Given two graph sets S1 and S2, Sub(S1,S2) = IsoFree(Union(S1,NCopy(S2))).
 */
public class Sub implements Binary {

  public Set<Graph> apply(Set<Graph> left, Set<Graph> right, Parameters params) {

    Binary union = new Union();
    Unary isoFree = new IsoFree();
    Unary negation = new NCopy();
    return isoFree.apply(union.apply(left, negation.apply(right, params), params), params);
  }
}