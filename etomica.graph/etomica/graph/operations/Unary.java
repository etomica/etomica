package etomica.graph.operations;

import java.util.Set;

import etomica.graph.model.Graph;

public interface Unary {

  public Set<Graph> apply(Set<Graph> argument, Arguments params);
}
