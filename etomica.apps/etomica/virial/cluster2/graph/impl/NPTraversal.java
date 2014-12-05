/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.GraphFactory;
import etomica.virial.cluster2.graph.GraphTraversal;
import etomica.virial.cluster2.graph.Nodes;
import etomica.virial.cluster2.graph.NodesVisitor;
import etomica.virial.cluster2.util.BitmapUtils;

public class NPTraversal extends BCTraversal {

  /**
   * A variation of BCTraversal that only traverses connected components
   * starting with field nodes.
   */
  @Override
  public void traverseAll(Nodes nodes, Edges edges, NodesVisitor visitor) {

    if (setup(nodes, edges, visitor)) {
      int nodeID = 0;
      while (nodeID < nodes.count()) {
        // start the next traversal with the first unseen field node
        while (nodeID < nodes.count()
            && (((getSeen() & BitmapUtils.bitOnMask(nodeID)) > 0) || (nodes.getAttributes(nodeID)
                .isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)))) {
          nodeID++;
        }
        // make sure the nodeID is valid
        if (nodeID < nodes.count()) {
          // every connected component traversal starts at time 0;
          // this ensures that every connected component has its
          // own traversal root
          time = 0;
          // traverse all bicomponents starting with a field node
          localVisit(GraphTraversal.START_COMPONENT);
          traverseBCC(nodeID, nodes, edges, true);
          localVisit(GraphTraversal.VISITED_COMPONENT);
        }
      }
    }
    if (seenAll()) {
      visit(VISITED_ALL);
    }
  }
}
