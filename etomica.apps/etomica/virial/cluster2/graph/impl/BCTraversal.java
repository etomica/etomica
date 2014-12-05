/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph.impl;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.GraphTraversal;
import etomica.virial.cluster2.graph.Nodes;
import etomica.virial.cluster2.graph.NodesVisitor;
import etomica.virial.cluster2.util.BitmapUtils;

public class BCTraversal extends AbstractGraphTraversal {

  protected int time;
  protected int[] visited;
  protected int[] low;
  protected LinkedList<NodePair> edgesStack = new LinkedList<NodePair>();

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

  @Override
  protected boolean setup(Nodes nodes, Edges edges, NodesVisitor visitor) {

    super.setup(nodes, edges, visitor);
    visited = new int[nodes.count()];
    low = new int[nodes.count()];
    time = 0;
    return true;
  }

  protected boolean traverseBCC(int nodeID, Nodes nodes, Edges edges, boolean all) {

    low[nodeID] = ++time;
    visited[nodeID] = low[nodeID];
    // this is a singleton component
    if (edges.getOutDegree(nodeID) == 0) {
      localVisit(GraphTraversal.START_BICOMPONENT);
      localVisit(nodeID);
      seen(nodeID);
      localVisit(GraphTraversal.VISITED_BICOMPONENT);
    }
    else {
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
            localVisit(GraphTraversal.START_BICOMPONENT);
            Set<Integer> bicomponent = new HashSet<Integer>();
            NodePair current = new NodePair(nodeID, neighbor);
            while (!edgesStack.getLast().equals(current)) {
              NodePair seenEdge = edgesStack.removeLast();
              bicomponent.add(seenEdge.fromNode);
              bicomponent.add(seenEdge.toNode);
            }
            NodePair seenEdge = edgesStack.removeLast();
            bicomponent.add(seenEdge.fromNode);
            bicomponent.add(seenEdge.toNode);
            for (int node : bicomponent) {
              localVisit(node);
              seen(node);
            }
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
    }
    return true;
  }

  @Override
  public void traverseAll(Nodes nodes, Edges edges, NodesVisitor visitor) {

    if (setup(nodes, edges, visitor)) {
      while (!seenAll()) {
        int nodeID = BitmapUtils.leftmostBit(~getSeen());
        // every connected component traversal starts at time 0;
        // this ensures that every connected component has its
        // own traversal root
        time = 0;
        // traverses all bicomponents of a connected component
        localVisit(GraphTraversal.START_COMPONENT);
        traverseBCC(nodeID, nodes, edges, true);
        localVisit(GraphTraversal.VISITED_COMPONENT);
      }
      visit(VISITED_ALL);
    }
  }

  @Override
  protected void traverseComponent(int nodeID, Nodes nodes, Edges edges) {

    traverseBCC(nodeID, nodes, edges, false);
  }
}