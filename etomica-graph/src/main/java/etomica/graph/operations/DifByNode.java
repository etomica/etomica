/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Graph;
import static etomica.graph.model.Metadata.*;
import etomica.graph.model.Node;

public class DifByNode implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {

    assert (params instanceof DifParameters);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      Set<Graph> newSet = apply(g, (DifParameters) params);
      if (newSet != null) {
        result.addAll(newSet);
      }
    }
    Unary isoFree = new IsoFree();
    return isoFree.apply(result, null);
  }

  private Set<Graph> apply(Graph g, DifParameters params) {

    Set<Graph> result = new HashSet<Graph>();
    for (byte nodeId = 0; nodeId < g.nodeCount(); nodeId++) {
      Node node = g.getNode(nodeId);
      if (node.getType() == TYPE_NODE_FIELD && node.getColor() == params.color()) {
        Graph newGraph = g.copy();
        newGraph.getNode(nodeId).setType(TYPE_NODE_ROOT);
        result.add(newGraph);
      }
    }
    return result;
  }
}