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
