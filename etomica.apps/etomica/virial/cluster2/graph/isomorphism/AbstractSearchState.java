package etomica.virial.cluster2.graph.isomorphism;

import etomica.virial.cluster2.graph.Graph;

public abstract class AbstractSearchState implements SearchState {

  private Graph firstGraph;
  private Graph secondGraph;

  public AbstractSearchState(Graph g1, Graph g2) {

    firstGraph = g1;
    secondGraph = g2;
  }

  @Override
  public Graph getG1() {

    return firstGraph;
  }

  @Override
  public Graph getG2() {

    return secondGraph;
  }
}
