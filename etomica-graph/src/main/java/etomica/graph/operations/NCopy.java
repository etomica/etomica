/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Graph;

public class NCopy implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {

    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      Graph ncopy = g.copy();
      ncopy.coefficient().setNumerator(-g.coefficient().getNumerator());
      result.add(ncopy);
    }
    return result;
  }
}