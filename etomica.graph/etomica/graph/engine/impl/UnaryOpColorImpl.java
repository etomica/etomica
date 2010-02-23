package etomica.graph.engine.impl;

import etomica.graph.engine.Parser.Expression;
import etomica.graph.engine.Parser.UnaryOpColor;

public class UnaryOpColorImpl extends UnaryOpImpl implements UnaryOpColor {

  private char color;

  public UnaryOpColorImpl(String operation, Expression expression, char color) {

    super(operation, expression);
    this.color = color;
  }

  public char getColor1() {

    return this.color;
  }
}