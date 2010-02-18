package etomica.graph.traversal;

import etomica.graph.model.Graph;
import static etomica.graph.model.Metadata.*;

public class BiconectedNodalPoint extends Biconnected {

  /**
   * A variation of BCTraversal that only traverses connected components starting with
   * field nodes.
   */
  @Override
  public void traverseAll(Graph graph, TraversalVisitor visitor) {

    if (setup(graph, visitor)) {
      byte nodeID = 0;
      while (nodeID < graph.getNodeCount()) {
        // start the next traversal with the first unseen field node
        while (nodeID < graph.getNodeCount()
            && (((getSeen() & BitmapUtils.bitOnMask(nodeID)) > 0) || (graph.nodes().get(nodeID).getType() == TYPE_NODE_FIELD))) {
          nodeID++;
        }
        // make sure the nodeID is valid
        if (nodeID < graph.getNodeCount()) {
          // every connected component traversal starts at time 0;
          // this ensures that every connected component has its
          // own traversal root
          time = 0;
          // traverse all bicomponents starting with a field node
          localVisit(STATUS_START_COMPONENT);
          traverseBCC(nodeID, graph, true);
          localVisit(STATUS_VISITED_COMPONENT);
        }
      }
    }
    if (seenAll()) {
      visit(STATUS_VISITED_ALL);
    }
  }
}
