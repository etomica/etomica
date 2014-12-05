/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Coefficient;
import etomica.graph.model.Graph;

/**
 * Multiply graphs by a scalar (coefficient)
 * 
 * @author Andrew Schultz
 */
public class MulScalar implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {
    assert (params instanceof MulScalarParameters);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      Graph newGraph = apply(g, (MulScalarParameters) params);
      result.add(newGraph);
    }
    return result;
  }

  public Graph apply(Graph argument, MulScalarParameters params) {
    Graph result = argument.copy();
    Coefficient c = result.coefficient();
    c.multiply(params.factor());
    return result;
  }
}
