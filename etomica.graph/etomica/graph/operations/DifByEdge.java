/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import static etomica.graph.model.Metadata.*;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Edge;
import etomica.graph.model.Graph;

public class DifByEdge implements Unary {

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
    return isoFree.apply(result, params);
  }

  public Set<Graph> apply(Graph g, DifParameters params) {

    Set<Graph> result = new HashSet<Graph>();
    for (Edge edge : g.edges()) {
      byte fromNode = g.getFromNode(edge.getId());
      byte toNode = g.getToNode(edge.getId());
      if (g.getNode(fromNode).getType() == TYPE_NODE_FIELD && g.getNode(toNode).getType() == TYPE_NODE_FIELD
          && edge.getColor() == params.color()) {
        Graph newGraph = g.copy();
        newGraph.getNode(fromNode).setType(TYPE_NODE_ROOT);
        newGraph.getNode(toNode).setType(TYPE_NODE_ROOT);
        result.add(newGraph);
      }
    }
    return result;
  }
}