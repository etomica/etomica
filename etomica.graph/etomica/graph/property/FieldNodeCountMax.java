/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.property;

import etomica.graph.model.Graph;
import etomica.graph.model.Metadata;
import etomica.graph.model.Node;

public class FieldNodeCountMax implements Property {

  protected final int maxCount;

  public FieldNodeCountMax(int maxCount) {
    this.maxCount = maxCount;
  }

  public boolean check(Graph graph) {
    int fieldCount = 0;
    for (Node node : graph.nodes()) {
      if (node.getType() == Metadata.TYPE_NODE_FIELD) {
        fieldCount++;
      }
    }
    return fieldCount <= maxCount;
  }
}
