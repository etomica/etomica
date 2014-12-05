/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine.impl;

import etomica.graph.engine.Parser.BinaryOp;
import etomica.graph.engine.Parser.Expression;

public class BinaryOpImpl implements BinaryOp {

  private Expression expression1;
  private Expression expression2;
  private String operation;

  public BinaryOpImpl(String operation, Expression expression1, Expression expression2) {

    this.operation = operation;
    this.expression1 = expression1;
    this.expression2 = expression2;
  }

  public Expression getExpression1() {

    return expression1;
  }

  public Expression getExpression2() {

    return expression2;
  }

  public String getOperation() {

    return operation;
  }
}