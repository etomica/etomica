package etomica.graph.operations;

import java.util.Set;

import etomica.graph.model.Graph;

/*
 * Union (GraphSet Union): GraphSet x GraphSet --> GraphSet
 *
 * Given two graph sets S1 and S2, Union(S1, S2) is the set union of the positive copy of S1 and the positive copy of S2.
 */
public class Union implements Binary {

  public Set<Graph> apply(Set<Graph> left, Set<Graph> right, Parameters params) {

    PCopy pcopy = new PCopy();
    Set<Graph> result = pcopy.apply(left, params);
    result.addAll(pcopy.apply(right, params));
    return result;
  }
}