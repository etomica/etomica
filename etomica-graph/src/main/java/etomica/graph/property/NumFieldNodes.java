/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.property;

import etomica.graph.model.Graph;
import etomica.graph.model.Metadata;

/**
 * Simple class that returns the number of field nodes.
 */
public class NumFieldNodes {

  public static int value(Graph g) {
    int n = 0, total = g.nodeCount();
    for (byte i=0; i<total; i++) {
      if (g.getNode(i).getType() == Metadata.TYPE_NODE_FIELD) n++;
    }
    return n;
  }
}
