package etomica.graph.engine.impl;

import etomica.graph.engine.Parser.Assignment;
import etomica.graph.engine.Parser.StrictExpression;
import etomica.graph.engine.Parser.Variable;

public class AssignmentImpl implements Assignment {

  private StrictExpression expression;
  private Variable variable;

  public AssignmentImpl(Variable variable, StrictExpression expression) {

    this.variable = variable;
    this.expression = expression;
  }

  public StrictExpression getExpression() {

    return expression;
  }

  public Variable getVariable() {

    return variable;
  }
}