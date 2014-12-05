/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import static etomica.graph.model.Metadata.*;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Graph;

public class Int implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {

    assert (params instanceof DifParameters);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      Graph newGraph = apply(g, (IntParameters) params);
      if (newGraph != null) {
        result.add(newGraph);
      }
    }
    Unary isoFree = new IsoFree();
    return isoFree.apply(result, params);
  }

  public Graph apply(Graph g, IntParameters params) {

    Graph result = null;
    if (params.nodeId() < g.nodeCount() && g.getNode(params.nodeId()).getType() == TYPE_NODE_ROOT) {
      result = g.copy();
      result.getNode(params.nodeId()).setType(TYPE_NODE_FIELD);
    }
    return result;
  }
}