/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.traversal;

import etomica.graph.model.Graph;
import static etomica.graph.model.Metadata.*;

public class BiconectedNodalPoint extends Biconnected {

  /**
   * A variation of BCTraversal that only traverses connected components starting with
   * field nodes.
   */
  @Override
  public byte traverseAll(Graph graph, TraversalVisitor visitor) {

    byte result = 0;
    if (setup(graph, visitor)) {
      byte nodeID = 0;
      while (nodeID < graph.nodeCount()) {
        // start the next traversal with the first unseen field node
        while (nodeID < graph.nodeCount()
            && (((getSeen() & BitmapUtils.bitOnMask(nodeID)) > 0) || (graph.nodes().get(nodeID).getType() == TYPE_NODE_FIELD))) {
          nodeID++;
        }
        // make sure the nodeID is valid
        if (nodeID < graph.nodeCount()) {
          // every connected component traversal starts at time 0;
          // this ensures that every connected component has its
          // own traversal root
          time = 0;
          // traverse all bicomponents starting with a field node
          localVisit(STATUS_START_COMPONENT);
          traverseBCC(nodeID, graph, true);
          localVisit(STATUS_VISITED_COMPONENT);
          result++;
        }
      }
    }
    if (seenAll()) {
      visit(STATUS_VISITED_ALL);
    }
    return result;
  }
}
