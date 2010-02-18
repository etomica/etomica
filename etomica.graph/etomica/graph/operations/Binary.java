package etomica.graph.operations;

import java.util.Set;

import etomica.graph.model.Graph;

public interface Binary {

  public Set<Graph> apply(Set<Graph> left, Set<Graph> right, Arguments params);
}