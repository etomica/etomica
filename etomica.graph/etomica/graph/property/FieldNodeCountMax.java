package etomica.graph.property;

import etomica.graph.model.Graph;
import etomica.graph.model.Node;

public class FieldNodeCountMax implements Property {

  protected final int maxCount;

  public FieldNodeCountMax(int maxCount) {
    this.maxCount = maxCount;
  }

  public boolean check(Graph graph) {
    int fieldCount = 0;
    for (Node node : graph.nodes()) {
      if (node.getType() == 'F') {
        fieldCount++;
      }
    }
    return fieldCount <= maxCount;
  }
}
