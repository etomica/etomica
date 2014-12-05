/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import static etomica.graph.model.Metadata.TYPE_NODE_ROOT;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Graph;

public class Conv implements Binary {

  public Set<Graph> apply(Set<Graph> left, Set<Graph> right, Parameters params) {

    assert (params instanceof ConvParameters);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph lg : left) {
      for (Graph rg : right) {
        Graph newGraph = apply(lg, rg, (ConvParameters) params);
        result.add(newGraph);
      }
    }
    Unary isoFree = new IsoFree();
    return isoFree.apply(result, null);
  }

  public Graph apply(Graph lg, Graph rg, ConvParameters params) {

    Graph result = null;
    if (params.nodeId() < lg.nodeCount() && params.nodeId() < rg.nodeCount()
        && lg.getNode(params.nodeId()).getType() == TYPE_NODE_ROOT
        && rg.getNode(params.nodeId()).getType() == TYPE_NODE_ROOT) {
      Mul multiply = new Mul();
      Int integration = new Int();
      return integration.apply(multiply.apply(lg, rg, params.mulParameters()), new IntParameters(params.nodeId()));
    }
    return result;
  }
}