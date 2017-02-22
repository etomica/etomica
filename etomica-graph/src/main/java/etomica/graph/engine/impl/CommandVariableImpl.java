/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine.impl;

import etomica.graph.engine.Parser.CommandVariable;
import etomica.graph.engine.Parser.Variable;

public class CommandVariableImpl extends CommandImpl implements CommandVariable {

  private Variable variable;

  public CommandVariableImpl(String command, Variable variable) {

    super(command);
    this.variable = variable;
  }

  public Variable getVariable() {

    return this.variable;
  }
}
