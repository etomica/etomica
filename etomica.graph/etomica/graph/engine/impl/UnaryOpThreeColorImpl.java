/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine.impl;

import etomica.graph.engine.Parser.Expression;
import etomica.graph.engine.Parser.UnaryOpThreeColor;

public class UnaryOpThreeColorImpl extends UnaryOpImpl implements UnaryOpThreeColor {

  private char color1;
  private char color2;
  private char color3;

  public UnaryOpThreeColorImpl(String operation, Expression expression, char color1, char color2, char color3) {

    super(operation, expression);
    this.color1 = color1;
    this.color2 = color2;
    this.color3 = color3;
  }

  public char getColor1() {

    return color1;
  }

  public char getColor2() {

    return color2;
  }

  public char getColor3() {

    return color3;
  }
}