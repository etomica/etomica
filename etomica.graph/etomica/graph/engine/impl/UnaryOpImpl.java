package etomica.graph.engine.impl;

import etomica.graph.engine.Parser.Expression;
import etomica.graph.engine.Parser.UnaryOp;

public class UnaryOpImpl implements UnaryOp {

  private Expression expression;
  private String operation;

  public UnaryOpImpl(String operation, Expression expression) {

    this.operation = operation;
    this.expression = expression;
  }

  public Expression getExpression() {

    return expression;
  }

  public String getOperation() {

    return operation;
  }
}