/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.traversal;

import etomica.graph.model.Graph;

public abstract class AbstractTraversal implements Traversal {

  private int seenNodes;
  private int goalNodes;
  private TraversalVisitor localVisitor;

  protected TraversalVisitor getLocalVisitor() {

    return localVisitor;
  }

  public int getSeen() {

    return seenNodes;
  }

  public int getGoal() {

    return goalNodes;
  }

  protected void setLocalVisitor(TraversalVisitor visitor) {

    localVisitor = visitor;
  }

  public void setSeen(int bitmap) {

    seenNodes = bitmap;
  }

  public void setGoal(int bitmap) {

    goalNodes = bitmap;
  }

  protected void error(byte errorID) {

    if (getLocalVisitor() != null) {
      getLocalVisitor().visit(NODE_NULL, errorID);
    }
  }

  protected void seen(byte nodeID) {

    setSeen(getSeen() | BitmapUtils.bitOnMask(nodeID));
  }

  protected boolean seenAll() {

    return (getSeen() == getGoal());
  }

  protected boolean unseenNeighbor(byte nodeID, byte neighborID) {

    return ((getSeen() & BitmapUtils.bitOnMask(neighborID)) == 0);
  }

  protected void status(byte statusID) {

    if (getLocalVisitor() != null) {
      getLocalVisitor().visit(NODE_NULL, statusID);
    }
  }

  protected void localVisit(byte nodeID) {

    if (getLocalVisitor() != null) {
      getLocalVisitor().visit(nodeID, STATUS_VISITED_NODE);
    }
  }

  protected void visit(byte nodeID) {

    localVisit(nodeID);
    if (nodeID >= 0) {
      seen(nodeID);
    }
  }

  protected boolean setup(Graph graph, TraversalVisitor visitor) {

    setLocalVisitor(visitor);
    if ((graph == null) || (graph.nodeCount() == 0)) {
      error(ERROR_TRAVERSAL_NODES);
      return false;
    }
    // we have seen no nodes
    setSeen(0);
    // one bit per every node we may see
    setGoal(BitmapUtils.bitMask(graph.nodeCount()));
    return true;
  }

  protected boolean setup(byte nodeID, Graph graph, TraversalVisitor visitor) {

    if (!setup(graph, visitor)) {
      return false;
    }
    if ((nodeID < 0) || (nodeID >= graph.nodeCount())) {
      error(ERROR_TRAVERSAL_ROOT);
      return false;
    }
    return true;
  }

  protected abstract void traverseComponent(byte nodeID, Graph graph);

  /**
   * There should be no reason to override this method.
   */
  public final boolean traverseComponent(byte nodeID, Graph graph, TraversalVisitor visitor) {

    if (setup(nodeID, graph, visitor)) {
      traverseComponent(nodeID, graph);
      return seenAll();
    }
    return false;
  }

  /**
   * The only reason to override this method is to change the node choice before
   * traversing each component, or to traverse more complex components such as
   * biconnected, K-connected, etc.
   */
  public byte traverseAll(Graph graph, TraversalVisitor visitor) {

    byte result = 0;
    if (setup(graph, visitor)) {
      while (!seenAll()) {
        byte nodeID = BitmapUtils.leftmostBit(~getSeen());
        // traverse the component
        traverseComponent(nodeID, graph);
        result++;
      }
      status(STATUS_VISITED_ALL);
    }
    return result == 0 ? graph.nodeCount() : result;
  }
}
