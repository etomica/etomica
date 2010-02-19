package etomica.graph.operations;

public class PowParameters implements Parameters {

  private byte exponent;

  public PowParameters(byte exponent) {

    this.exponent = exponent;
  }

  public byte exponent() {

    return exponent;
  }
}
