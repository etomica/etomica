package etomica.graph.property;

import etomica.graph.model.Graph;

public class HasEdgeCount implements Property {

  private byte edgeCount;

  public HasEdgeCount(byte edgeCount) {

    this.edgeCount = edgeCount;
  }

  public boolean check(Graph graph) {

    return graph.edgeCount() == edgeCount;
  }
}
