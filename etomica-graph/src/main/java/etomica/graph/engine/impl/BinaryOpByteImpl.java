/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine.impl;

import etomica.graph.engine.Parser.BinaryOpByte;
import etomica.graph.engine.Parser.Expression;

public class BinaryOpByteImpl extends BinaryOpImpl implements BinaryOpByte {

  private byte value;

  public BinaryOpByteImpl(String operation, Expression expression1, Expression expression2, byte value) {

    super(operation, expression1, expression2);
    this.value = value;
  }

  public byte getByte1() {

    return this.value;
  }
}
