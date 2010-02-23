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