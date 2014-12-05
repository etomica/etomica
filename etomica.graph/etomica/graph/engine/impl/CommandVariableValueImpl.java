/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine.impl;

import etomica.graph.engine.Parser.CommandVariableValue;
import etomica.graph.engine.Parser.Value;
import etomica.graph.engine.Parser.Variable;

public class CommandVariableValueImpl extends CommandImpl implements CommandVariableValue {

  private Value value;
  private Variable variable;

  public CommandVariableValueImpl(String command, Variable variable, Value value) {

    super(command);
    this.variable = variable;
    this.value = value;
  }

  public Value getValue() {

    return this.value;
  }

  public Variable getVariable() {

    return this.variable;
  }
}
