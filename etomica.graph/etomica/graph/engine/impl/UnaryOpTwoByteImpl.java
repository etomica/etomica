package etomica.graph.engine.impl;

import etomica.graph.engine.Parser.Expression;
import etomica.graph.engine.Parser.UnaryOpTwoByte;

public class UnaryOpTwoByteImpl extends UnaryOpImpl implements UnaryOpTwoByte {

  private byte value1;
  private byte value2;

  public UnaryOpTwoByteImpl(String operation, Expression expression, byte value1, byte value2) {

    super(operation, expression);
    this.value1 = value1;
    this.value2 = value2;
  }

  public byte getByte1() {

    return value1;
  }

  public byte getByte2() {

    return value2;
  }
}