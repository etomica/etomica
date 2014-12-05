/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.traversal;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;

import etomica.graph.model.Graph;

public class Biconnected extends AbstractTraversal {

  protected int time;
  protected int[] visited;
  protected int[] low;
  protected LinkedList<NodePair> edgesStack = new LinkedList<NodePair>();

  private class NodePair {

    public byte fromNode;
    public byte toNode;

    public NodePair(byte x, byte y) {

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
  protected boolean setup(Graph graph, TraversalVisitor visitor) {

    super.setup(graph, visitor);
    visited = new int[graph.nodeCount()];
    low = new int[graph.nodeCount()];
    time = 0;
    edgesStack.clear();
    return true;
  }

  protected boolean traverseBCC(byte nodeID, Graph graph, boolean all) {

    low[nodeID] = ++time;
    visited[nodeID] = low[nodeID];
    // this is a singleton component
    byte od = graph.getOutDegree(nodeID);
    if (od == 0) {
      status(STATUS_START_BICOMPONENT);
      localVisit(nodeID);
      seen(nodeID);
      status(STATUS_VISITED_BICOMPONENT);
    }
    else {
      for (byte i = 0; i < od; i++) {
        NodePair lastEdge = null;
        if (!edgesStack.isEmpty()) {
          lastEdge = edgesStack.getLast();
        }
        byte neighbor = graph.getOutNode(nodeID, i);
        if ((lastEdge != null) && (neighbor == lastEdge.fromNode)) {
          continue;
        }
        if (visited[neighbor] < visited[nodeID]) {
          edgesStack.addLast(new NodePair(nodeID, neighbor));
        }
        if (visited[neighbor] == 0) {
          if (!traverseBCC(neighbor, graph, all)) {
            return false;
          }
          low[nodeID] = Math.min(low[nodeID], low[neighbor]);
          // a) if nodeID is the root (visited at time 1), and one of
          // its neighbors is visited from nodeID after time 2, then
          // the root node must be an articulation point;
          if ((visited[nodeID] == 1) && (visited[neighbor] != 2)) {
            status(STATUS_ARTICULATION_POINT);
            localVisit(nodeID);
          }
          // b) for every non-root node, if the low of a neighbor of
          // nodeID is not smaller than nodeID, then the entire
          // subtree rooted at neighbor has no back edges to some
          // proper ancestor of nodeID; therefore, nodeID must be
          // an articulation point;
          if ((visited[nodeID] != 1) && (low[neighbor] >= visited[nodeID])) {
            status(STATUS_ARTICULATION_POINT);
            localVisit(nodeID);
          }
          // traverse the entire bicomponent
          if (low[neighbor] >= visited[nodeID]) {
            status(STATUS_START_BICOMPONENT);
            Set<Byte> bicomponent = new HashSet<Byte>();
            NodePair current = new NodePair(nodeID, neighbor);
            while (!edgesStack.getLast().equals(current)) {
              NodePair seenEdge = edgesStack.removeLast();
              bicomponent.add(seenEdge.fromNode);
              bicomponent.add(seenEdge.toNode);
            }
            NodePair seenEdge = edgesStack.removeLast();
            bicomponent.add(seenEdge.fromNode);
            bicomponent.add(seenEdge.toNode);
            for (byte node : bicomponent) {
              localVisit(node);
              seen(node);
            }
            status(STATUS_VISITED_BICOMPONENT);
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
  public byte traverseAll(Graph graph, TraversalVisitor visitor) {

    byte result = 0;
    if (setup(graph, visitor)) {
      while (!seenAll()) {
        byte nodeID = BitmapUtils.leftmostBit(~getSeen());
        // every connected component traversal starts at time 0;
        // this ensures that every connected component has its
        // own traversal root
        time = 0;
        // traverses all bicomponents of a connected component
        status(STATUS_START_COMPONENT);
        traverseBCC(nodeID, graph, true);
        status(STATUS_VISITED_COMPONENT);
        result++;
      }
      status(STATUS_VISITED_ALL);
    }
    return result == 0 ? (byte) graph.nodes().size() : result;
  }

  @Override
  protected void traverseComponent(byte nodeID, Graph graph) {
    traverseBCC(nodeID, graph, false);
  }
}