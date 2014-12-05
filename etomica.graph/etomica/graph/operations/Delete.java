/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.isomorphism.Match;
import etomica.graph.model.Graph;

/*
 * Delete from left all graphs that have an isomorph in right.
 */
public class Delete implements Binary {

  public Set<Graph> apply(Set<Graph> left, Set<Graph> right, Parameters params) {

    Set<Graph> result = new HashSet<Graph>();
    for (Graph lg : left) {
      boolean deleteLeft = false;
      for (Graph rg : right) {
        // should we delete lg from the result?
        if (Match.match(lg, rg, true)) {
          deleteLeft = true;
          break;
        }
      }
      if (!deleteLeft) {
        result.add(lg.copy());
      }
    }
    return result;
  }
}