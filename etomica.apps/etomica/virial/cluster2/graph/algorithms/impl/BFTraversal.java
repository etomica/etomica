package etomica.virial.cluster2.graph.algorithms.impl;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.Nodes;

public class BFTraversal extends AbstractGraphTraversal {

  protected void traverseComponent(int nodeID, Nodes nodes, Edges edges) {

    QueueOfInt toExploreQ = new QueueOfInt(nodes.count());
    // visit the node and update the seen nodes
    visit(nodeID);
    // queue the visited node for exploration of its neighbors
    toExploreQ.enqueue(nodeID);
    // done when: (1) all nodes seen OR (2) no new nodes to explore
    while (!seenAll() && !toExploreQ.isEmpty()) {
      // remove the first node to explore from the queue
      int explore = toExploreQ.dequeue();
      // visit every unseen neighbor and enqueue them for traversal
      for (int neighbor = 0; neighbor < nodes.count(); neighbor++) {
        if (unseenNeighbor(explore, neighbor, edges)) {
          visit(neighbor);
          toExploreQ.enqueue(neighbor);
        }
      }
    }
  }

//    // nothing to traverse
//    if (nodes.count() == 0) {
//      return;
//    }
//    // visit the root node
//    visitor.visit(0);
//    // just one node to traverse
//    if (nodes.count() == 1) {
//      return;
//    }
//    // one bit for each node we have to see
//    int goal = BitmapUtils.bitMask(nodes.count());
//    // one bit for each node seen
//    int seen = 0x00000001;
//    // nodes to explore
//    QueueOfInt toExploreQ = new QueueOfInt(nodes.count());
//    toExploreQ.enqueue(0);
//    // done when: (1) all nodes seen OR (2) no new nodes to explore
//    while ((seen != goal) && !toExploreQ.isEmpty()) {
//      // retrieve the first node to explore (remove from the queue)
//      int explore = toExploreQ.dequeue();
//      boolean unseenNeighbor = false;
//      for (int node = 0; node < nodes.count(); node++) {
//        // check if this is an unseen neighbor of the node being explored
//        unseenNeighbor = (node != explore)
//            && ((seen & BitmapUtils.bitOnMask(node)) == 0)
//            && edges.hasEdge(explore, node);
//        // if node is an unseen neighbor, enqueue it for traversal
//        if (unseenNeighbor) {
//          seen = seen | BitmapUtils.bitOnMask(node);
//          toExploreQ.enqueue(node);
//          // visit the node
//          visitor.visit(node);
//        }
//      }
//    }
//    return;
}