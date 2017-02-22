/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Edge;
import etomica.graph.model.Graph;

/**
 * This operation multiplies two graphs (or each pair of graphs from two sets)
 * by superimposing all points (with the idea that all nodes are
 * distinguishable).  The result is a graph with all bonds from both graphs.
 * 
 * @author Andrew Schultz
 */
public class MulUnion implements Binary {

  public Set<Graph> apply(Set<Graph> left, Set<Graph> right, Parameters params) {
    assert(params == null);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph lg : left) {
      for (Graph rg : right) {
        Graph graph = apply(lg, rg);
        if (graph != null) {
          result.add(graph);
        }
      }
    }
    return result;
  }

  public Graph apply(Graph left, Graph right) {
    assert(left.nodeCount() == right.nodeCount());
    Graph result = left.copy();
    for (Edge e : right.edges()) {
      byte eid = e.getId();
      if (result.hasEdge(eid)) {
        throw new RuntimeException("can't union-multiply two graphs with edges in the same position");
      }
      result.putEdge(eid);
      result.getEdge(eid).setColor(e.getColor());
    }
    return result;
  }

}
