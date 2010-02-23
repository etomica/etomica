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
