package etomica.graph.operations;

import java.util.Set;

import etomica.graph.model.Graph;

public class Union implements Binary {

  public Set<Graph> apply(Set<Graph> left, Set<Graph> right, Parameters params) {

    PCopy pcopy = new PCopy();
    Set<Graph> result = pcopy.apply(left, params);
    result.addAll(pcopy.apply(right, params));
    return result;
  }
}