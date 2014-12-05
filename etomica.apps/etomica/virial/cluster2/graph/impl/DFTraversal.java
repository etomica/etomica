/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph.impl;

import java.util.LinkedList;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.GraphTraversal;
import etomica.virial.cluster2.graph.Nodes;

public class DFTraversal extends AbstractGraphTraversal {

  @Override
  protected void traverseComponent(int nodeID, Nodes nodes, Edges edges) {

    LinkedList<Integer> toExploreQ = new LinkedList<Integer>();
    // start a component traversal
    visit(GraphTraversal.START_COMPONENT);
    // visit the node and update the seen nodes
    visit(nodeID);
    // queue the visited node for exploration of its neighbors
    toExploreQ.add(nodeID);
    // done when:
    //   (1) all nodes seen OR
    //   (2) no new nodes to explore in this connected component
    while (!seenAll() && !toExploreQ.isEmpty()) {
      // retrieve the next node to explore from the queue (DO NOT REMOVE!!!)
      int explore = toExploreQ.getLast();
      // does the node we are exploring have an unseen neighbor?
      boolean unseenNeighbor = false;
      // visit the first unseen neighbor (if any) and enqueue it for traversal
      for (int i = 0; i < edges.getOutDegree(explore); i++) {
        int neighbor = edges.getOutNode(explore, i);
        if (unseenNeighbor(explore, neighbor)) {
          visit(neighbor);
          toExploreQ.add(neighbor);
          unseenNeighbor = true;
          // in DF, we follow the first unseen neighbor
          break;
        }
      }
      // dequeue a node only when all its neighbors are seen;
      // early dequeuing breaks backtracking of the DF traversal;
      if (!unseenNeighbor) {
        toExploreQ.removeLast();
      }
    }
    // this component has been visited
    visit(VISITED_COMPONENT);
  }
}