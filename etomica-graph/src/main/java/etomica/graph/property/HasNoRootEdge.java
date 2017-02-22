/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.property;

import static etomica.graph.model.Metadata.TYPE_NODE_ROOT;

import etomica.graph.model.Edge;
import etomica.graph.model.Graph;

public class HasNoRootEdge implements Property {

  public boolean check(Graph graph) {

    for (Edge edge : graph.edges()) {
      if (graph.getNode(graph.getFromNode(edge.getId())).getType() == TYPE_NODE_ROOT
          && graph.getNode(graph.getToNode(edge.getId())).getType() == TYPE_NODE_ROOT) {
        return false;
      }
    }
    return true;
  }
}