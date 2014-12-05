/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.GraphTraversal;
import etomica.virial.cluster2.graph.Nodes;
import etomica.virial.cluster2.graph.NodesVisitor;
import etomica.virial.cluster2.util.BitmapUtils;

public abstract class AbstractGraphTraversal implements GraphTraversal {

  private int seenNodes;
  private int goalNodes;
  private NodesVisitor localVisitor;

  protected NodesVisitor getLocalVisitor() {

    return localVisitor;
  }

  public int getSeen() {

    return seenNodes;
  }

  public int getGoal() {

    return goalNodes;
  }

  protected void setLocalVisitor(NodesVisitor visitor) {

    localVisitor = visitor;
  }

  public void setSeen(int bitmap) {

    seenNodes = bitmap;
  }

  public void setGoal(int bitmap) {

    goalNodes = bitmap;
  }

  protected void error(int errorID) {

    if (getLocalVisitor() != null) {
      getLocalVisitor().visit(errorID);
    }
  }

  protected void seen(int nodeID) {

    setSeen(getSeen() | BitmapUtils.bitOnMask(nodeID));
  }

  protected boolean seenAll() {

    return (getSeen() == getGoal());
  }

  protected boolean unseenNeighbor(int nodeID, int neighborID) {

    return ((getSeen() & BitmapUtils.bitOnMask(neighborID)) == 0);
  }

  protected void localVisit(int nodeID) {

    if (getLocalVisitor() != null) {
      getLocalVisitor().visit(nodeID);
    }
  }

  protected void visit(int nodeID) {

    localVisit(nodeID);
    if (nodeID >= 0) {
      seen(nodeID);
    }
  }

  protected boolean setup(Nodes nodes, Edges edges, NodesVisitor visitor) {

    setLocalVisitor(visitor);
    if ((nodes == null) || (nodes.count() == 0)) {
      error(TRAVERSAL_NODES_ERROR);
      return false;
    }
    if (edges == null) {
      error(TRAVERSAL_EDGES_ERROR);
      return false;
    }
    // we have seen no nodes
    setSeen(0);
    // one bit per every node we may see
    setGoal(BitmapUtils.bitMask(nodes.count()));
    return true;
  }

  protected boolean setup(int nodeID, Nodes nodes, Edges edges,
      NodesVisitor visitor) {

    if (!setup(nodes, edges, visitor)) {
      return false;
    }
    if ((nodeID < 0) || (nodeID >= nodes.count())) {
      error(TRAVERSAL_ROOT_ERROR);
      return false;
    }
    return true;
  }

  protected abstract void traverseComponent(int nodeID, Nodes nodes, Edges edges);

  /**
   * There should be no reason to override this method.
   */
  public final boolean traverseComponent(int nodeID, Nodes nodes, Edges edges,
      NodesVisitor visitor) {

    if (setup(nodeID, nodes, edges, visitor)) {
      traverseComponent(nodeID, nodes, edges);
      return seenAll();
    }
    return false;
  }

  /**
   * The only reason to override this method is to change the node choice before
   * traversing each component, or to traverse more complex components such as
   * biconnected, triconnected, etc.
   */
  public void traverseAll(Nodes nodes, Edges edges, NodesVisitor visitor) {

    if (setup(nodes, edges, visitor)) {
      while (!seenAll()) {
        int nodeID = BitmapUtils.leftmostBit(~getSeen());
        // traverse the component
        traverseComponent(nodeID, nodes, edges);
      }
      visit(VISITED_ALL);
    }
  }
}