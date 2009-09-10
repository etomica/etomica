package etomica.virial.cluster2.graph.impl;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.Nodes;
import etomica.virial.cluster2.graph.NodesVisitor;
import etomica.virial.cluster2.util.BitmapUtils;

public class BCTraversal extends AbstractGraphTraversal {

  private int time;
  private boolean[] articulation;
  private int[] visited;
  private int[] low;
  private Map<Integer,Set<Integer>> edgesMap = new HashMap<Integer,Set<Integer>>();

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

    return edgesMap.get(explore).contains(neighbor) || edgesMap.get(neighbor).contains(explore);
  }

  protected boolean setup(Nodes nodes, Edges edges, NodesVisitor visitor) {

    super.setup(nodes, edges, visitor);
    visited = new int[nodes.count()];
    low = new int[nodes.count()];
    articulation = new boolean[nodes.count()];
    edgesMap.clear();
    for (int nodeID = 0; nodeID < nodes.count(); nodeID++) {
      edgesMap.put(nodeID, new HashSet<Integer>());
    }
    time = 0;
    return true;
  }

  protected void traverseCC(int nodeID, Nodes nodes, Edges edges) {

    int localSeen = getSeen();
    LinkedList<Integer> toExploreQ = new LinkedList<Integer>();
    // visit the node and update the seen nodes
    visit(nodeID);
    // queue the visited node for exploration of its neighbors
    toExploreQ.add(nodeID);
    // done when:
    // (1) all nodes seen OR
    // (2) no new nodes to explore in this connected component
    while (!seenAll() && !toExploreQ.isEmpty()) {
      // retrieve the next node to explore from the queue (DO NOT REMOVE!!!)
      int explore = toExploreQ.peekLast();
      // does the node we are exploring have an unseen neighbor?
      boolean unseenNeighbor = false;
      // visit the first unseen neighbor (if any) and enqueue it for traversal
      for (int i = 0; i < edges.getOutDegree(explore); i++) {
        int neighbor = edges.getOutNode(explore, i);
        // each edge should be processed once
        if (seenEdge(explore, neighbor)) {
          continue;
        }
        if (unseenNeighbor(explore, neighbor)) {
          visit(neighbor);
          visitEdge(explore, neighbor);
          toExploreQ.add(neighbor);
          unseenNeighbor = true;
          // in DF, we follow the first unseen neighbor
          break;
        }
        else {
          low[explore] = Math.min(low[explore], visited[neighbor]);
        }
      }
      // dequeue a node only when all its neighbors are seen;
      // early dequeuing breaks backtracking of the DF traversal;
      if (!unseenNeighbor) {
        int removed = toExploreQ.removeLast();
        // when we remove the root, break
        if (toExploreQ.isEmpty()) {
          break;
        }
        int backNeighbor = toExploreQ.peekLast();
        low[backNeighbor] = Math.min(low[backNeighbor], low[removed]);
        // backNeighbor is the root and an articulation point
        if ((visited[backNeighbor] == 1) && (visited[removed] != 2)) {
          articulation[backNeighbor] = true;
        }
        // backNeighbor is a non-root articulation point
        if ((visited[backNeighbor] != 1)
            && (low[removed] >= visited[backNeighbor])) {
          articulation[backNeighbor] = true;
        }
      }
    }
    /**
     * 
     * 
     * 
     */
    // process the neighbors of the last node added to the queue
    if (seenAll() && !toExploreQ.isEmpty()) {
      // retrieve the next node to explore from the queue (DO NOT REMOVE!!!)
      int explore = toExploreQ.peekLast();
      // traverse all edges that have not been traversed yet
      for (int i = 0; i < edges.getOutDegree(explore); i++) {
        int neighbor = edges.getOutNode(explore, i);
        // each edge should be processed once
        if (seenEdge(explore, neighbor)) {
          continue;
        }
        low[explore] = Math.min(low[explore], visited[neighbor]);
      }
    }
    /**
     * 
     * 
     * 
     */
    // when we remove the root, break
    while (!toExploreQ.isEmpty()) {
      int removed = toExploreQ.removeLast();
      // when we remove the root, break
      if (toExploreQ.isEmpty()) {
        break;
      }
      int backNeighbor = toExploreQ.peekLast();
      low[backNeighbor] = Math.min(low[backNeighbor], low[removed]);
      // backNeighbor is the root and an articulation point
      if ((visited[backNeighbor] == 1) && (visited[removed] != 2)) {
        articulation[backNeighbor] = true;
      }
      // backNeighbor is a non-root articulation point
      if ((visited[backNeighbor] != 1)
          && (low[removed] >= visited[backNeighbor])) {
        articulation[backNeighbor] = true;
      }
    }
    // reset the nodes we've already seen
    setSeen(localSeen);
  }

  protected void traverseBCC(int nodeID, Nodes nodes, Edges edges, boolean all) {

    traverseCC(nodeID, nodes, edges);
    // traverse one bicomponent of this connected component
    for (int visitTime = visited[nodeID]; visitTime <= time; visitTime++) {
      for (int visitedID = 0; visitedID < visited.length; visitedID++) {
        // identify the node visited at time visitTime
        if (visited[visitedID] == visitTime) {
          // we've actually seen this node
          seen(visitedID);
          // visit the next biconnected component node
          localVisit(visitedID);
          // every articulation point is part of two or
          // more bicomponents, so we may need to visit
          // it more than once
          if (articulation[visitedID]) {
            localVisit(VISITED_BC_COMPONENT);
            if (!all) {
              return;
            }
            localVisit(visitedID);
          }
          // if this is the last visited node, then we
          // must finish visiting the current component
          else if (visitTime == time) {
            localVisit(VISITED_BC_COMPONENT);
            if (!all) {
              return;
            }
          }
          // short-circuit the traversal of the visited nodes
          break;
        }
      }
    }
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
