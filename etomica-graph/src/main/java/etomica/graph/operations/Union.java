/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.Set;

import etomica.graph.model.Graph;

public class Union implements Binary {

  public Set<Graph> apply(Set<Graph> left, Set<Graph> right, Parameters params) {

    PCopy pcopy = new PCopy();
    Set<Graph> result = pcopy.apply(left, params);
    result.addAll(pcopy.apply(right, params));
    return result;
  }
}