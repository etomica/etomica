/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.Set;

import etomica.graph.model.Graph;

public class Exp implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {

    assert (params instanceof ExpParameters);
    Unary power = new Pow();
    Binary sum = new Sum();
    Set<Graph> result = null;
    for (byte pow = ((ExpParameters) params).expLo(); pow < ((ExpParameters) params).expHi(); pow++) {
      Set<Graph> newSet = power.apply(argument, new PowParameters(pow));
      for (Graph g : newSet) {
        g.coefficient().setDenominator(g.coefficient().getDenominator() * fact(pow));
      }
      if (result == null) {
        result = newSet;
      }
      else {
        result = sum.apply(result, newSet, params);
      }
    }
    return result;
  }

  private int fact(byte k) {

    int result = 1;
    for (int n = 2; n <= k; n++) {
      result = result * n;
    }
    return result;
  }
}