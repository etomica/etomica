package etomica.virial.cluster2.graph.impl;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.GraphTraversal;
import etomica.virial.cluster2.graph.Nodes;
import etomica.virial.cluster2.graph.NodesVisitor;
import etomica.virial.cluster2.util.BitmapUtils;

public class BCTraversal extends AbstractGraphTraversal {

  private int time;
  private int[] visited;
  private int[] low;
  private Map<Integer, Set<Integer>> edgesMap = new HashMap<Integer, Set<Integer>>();
  private LinkedList<NodePair> edgesStack = new LinkedList<NodePair>();

  private class NodePair {

    public int fromNode;
    public int toNode;

    public NodePair(int x, int y) {

      fromNode = x;
      toNode = y;
    }

    @Override
    public boolean equals(Object other) {

      if (other instanceof NodePair) {
        NodePair p = (NodePair) other;
        return p.fromNode == fromNode && p.toNode == toNode;
      }
      return false;
    }
  }

  protected void visit(int nodeID) {

    if (nodeID >= 0) {
      seen(nodeID);
      visited[nodeID] = ++time;
      low[nodeID] = visited[nodeID];
    }
  }

  private void visitEdge(int explore, int neighbor) {

    edgesMap.get(explore).add(neighbor);
  }

  private boolean seenEdge(int explore, int neighbor) {

    return edgesMap.get(explore).contains(neighbor)
        || edgesMap.get(neighbor).contains(explore);
  }

  protected boolean setup(Nodes nodes, Edges edges, NodesVisitor visitor) {

    super.setup(nodes, edges, visitor);
    visited = new int[nodes.count()];
    low = new int[nodes.count()];
    time = 0;
    return true;
  }

  protected boolean traverseBCC(int nodeID, Nodes nodes, Edges edges,
      boolean all) {

    low[nodeID] = ++time;
    visited[nodeID] = low[nodeID];
    for (int i = 0; i < edges.getOutDegree(nodeID); i++) {
      NodePair lastEdge = null;
      if (!edgesStack.isEmpty()) {
        lastEdge = edgesStack.getLast();
      }
      int neighbor = edges.getOutNode(nodeID, i);
      if ((lastEdge != null) && (neighbor == lastEdge.fromNode)) {
        continue;
      }
      if (visited[neighbor] < visited[nodeID]) {
        edgesStack.addLast(new NodePair(nodeID, neighbor));
      }
      if (visited[neighbor] == 0) {
        if (!traverseBCC(neighbor, nodes, edges, all)) {
          return false;
        }
        low[nodeID] = Math.min(low[nodeID], low[neighbor]);
        // a) if nodeID is the root (visited at time 1), and one of
        // its neighbors is visited from nodeID after time 2, then
        // the root node must be an articulation point;
        if ((visited[nodeID] == 1) && (visited[neighbor] != 2)) {
          localVisit(GraphTraversal.ARTICULATION_POINT);
          localVisit(nodeID);
        }
        // b) for every non-root node, if the low of a neighbor of
        // nodeID is not smaller than nodeID, then the entire
        // subtree rooted at neighbor has no back edges to some
        // proper ancestor of nodeID; therefore, nodeID must be
        // an articulation point;
        if ((visited[nodeID] != 1) && (low[neighbor] >= visited[nodeID])) {
          localVisit(GraphTraversal.ARTICULATION_POINT);
          localVisit(nodeID);
        }
        // traverse the entire bicomponent
        if (low[neighbor] >= visited[nodeID]) {
          localVisit(GraphTraversal.VISIT_START);
          NodePair current = new NodePair(nodeID, neighbor);
          while (!edgesStack.getLast().equals(current)) {
            NodePair seenEdge = edgesStack.removeLast();
            localVisit(seenEdge.fromNode);
            seen(seenEdge.fromNode);
          }
          NodePair seenEdge = edgesStack.removeLast();
          localVisit(seenEdge.fromNode);
          localVisit(seenEdge.toNode);
          seen(seenEdge.fromNode);
          seen(seenEdge.toNode);
          localVisit(GraphTraversal.VISITED_BICOMPONENT);
          if (!all) {
            return false;
          }
        }
      }
      else {
        low[nodeID] = Math.min(low[nodeID], visited[neighbor]);
      }
    }
    return true;
  }

  public void traverseAll(Nodes nodes, Edges edges, NodesVisitor visitor) {

    if (setup(nodes, edges, visitor)) {
      while (!seenAll()) {
        int nodeID = BitmapUtils.leftmostBit(~getSeen());
        // traverses all bicomponents of a connected component
        traverseBCC(nodeID, nodes, edges, true);
      }
      visit(VISITED_ALL_COMPONENTS);
    }
  }

  protected void traverseComponent(int nodeID, Nodes nodes, Edges edges) {

    traverseBCC(nodeID, nodes, edges, false);
  }
}