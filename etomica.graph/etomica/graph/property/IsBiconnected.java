/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.property;

import etomica.graph.model.Graph;
import etomica.graph.traversal.Biconnected;
import etomica.graph.traversal.Traversal;

public class IsBiconnected implements Property {

  private Traversal bcTraversal;

  public IsBiconnected() {

    this.bcTraversal = new Biconnected();
  }

  public boolean check(Graph graph) {
    byte n = graph.nodeCount();
    // by definition, a null graph is not biconnected
    if (n == 0) {
      return false;
    }
    // by definition, a singleton graph is biconnected
    if (n == 1) {
      return true;
    }
    byte edgeCount = graph.edgeCount();
    // a two-node graph having a single edge is biconnected
    if (n == 2 && edgeCount == 1) {
      return true;
    }
    // invariant: a biconnected graph with N > 2 nodes has at least N edges
    //
    // SKETCH: Assume a graph G has N-1 edges. We show that this graph cannot
    // be biconnected. By definition, G is connected. We claim that G has no
    // cycles. Assume G has a simple cycle C_K involving K nodes. C_K has K
    // edges, so the remaining N-K nodes in G must be connected to C_K through
    // the remaining (N-1-K) edges, which is not possible. So G has no cycles.
    // By definition, G is a tree (a tree is a connected, cycle-free graph).
    // Any two nodes in G are connected by exactly one simple path so removing
    // ANY node with more than one incident edge disconnects the tree. Since
    // G has at least one node with two incident edges, G can be split. Q.E.D.
    if (edgeCount < n) {
      return false;
    }
    // traverse the graph starting at nodeID and return true IFF all nodes in the graph
    // are traversed in the same connected component
    return bcTraversal.traverseComponent((byte) 0, graph, null);
  }
}