package etomica.virial.cluster2.graph.impl;

import java.util.LinkedList;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.Nodes;

public class BFTraversal extends AbstractGraphTraversal {

  protected void traverseComponent(int nodeID, Nodes nodes, Edges edges) {

    LinkedList<Integer> toExploreQ = new LinkedList<Integer>();
    // visit the node and update the seen nodes
    visit(nodeID);
    // queue the visited node for exploration of its neighbors
    toExploreQ.add(nodeID);
    // done when: (1) all nodes seen OR (2) no new nodes to explore
    while (!seenAll() && !toExploreQ.isEmpty()) {
      // remove the first node to explore from the queue
      int explore = toExploreQ.remove();
      // visit every unseen neighbor and enqueue them for traversal
      for (int neighbor = 0; neighbor < nodes.count(); neighbor++) {
        if (unseenNeighbor(explore, neighbor, edges)) {
          visit(neighbor);
          toExploreQ.add(neighbor);
        }
      }
    }
  }
}