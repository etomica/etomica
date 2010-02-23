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