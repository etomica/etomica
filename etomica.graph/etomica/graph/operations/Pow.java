/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.Set;

import etomica.graph.model.Graph;

public class Pow implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {

    assert (params instanceof PowParameters);
    Unary pcopy = new PCopy();
    Binary multiply = new Mul();
    Set<Graph> result = pcopy.apply(argument, params);
    for (byte pow = 2; pow < ((PowParameters) params).exponent(); pow++) {
      result = multiply.apply(result, argument, params);
    }
    return result;
  }
}