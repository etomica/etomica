/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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