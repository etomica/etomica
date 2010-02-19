package etomica.graph.operations;

public final class RelabelParameters implements Parameters {

  private byte[] permutation;

  public RelabelParameters(byte[] permutation) {

    this.permutation = permutation;
  }

  public byte map(byte index) {

    return permutation[index];
  }

  public int size() {

    return permutation.length;
  }
}