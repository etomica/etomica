/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine.impl;

import etomica.graph.engine.Parser.Expression;
import etomica.graph.engine.Parser.UnaryOpByte;

public class UnaryOpByteImpl extends UnaryOpImpl implements UnaryOpByte {

  private byte value;

  public UnaryOpByteImpl(String operation, Expression expression, byte value) {

    super(operation, expression);
    this.value = value;
  }

  public byte getByte1() {

    return value;
  }
}