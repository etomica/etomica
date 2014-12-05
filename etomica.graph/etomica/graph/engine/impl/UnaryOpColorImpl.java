/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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