/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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