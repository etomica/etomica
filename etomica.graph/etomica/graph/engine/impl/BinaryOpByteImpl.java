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
