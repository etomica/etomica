package etomica.graph.operations;

public class ExpParameters implements Parameters {

  private byte expLo;
  private byte expHi;

  public ExpParameters(byte expLo, byte expHi) {

    this.expLo = expLo;
    this.expHi = expHi;
  }

  public byte expLo() {

    return expLo;
  }

  public byte expHi() {

    return expHi;
  }
}
