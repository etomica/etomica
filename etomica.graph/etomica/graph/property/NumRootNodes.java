/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.property;

import etomica.graph.model.Graph;
import etomica.graph.model.Metadata;
import etomica.graph.model.Node;

/**
 * Simple class that returns the number of root nodes.
 */
public class NumRootNodes {

  public static int value(Graph g) {
    int n = 0;
    for (Node node : g.nodes()) {
      if (node.getType() == Metadata.TYPE_NODE_ROOT) n++;
    }
    return n;
  }
}
