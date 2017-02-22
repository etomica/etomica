/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine.impl;

import etomica.graph.engine.Parser.Variable;

public class VariableImpl implements Variable {

  private String variable;

  public VariableImpl(String variable) {

    this.variable = variable;
  }

  public String getName() {

    return this.variable;
  }
}