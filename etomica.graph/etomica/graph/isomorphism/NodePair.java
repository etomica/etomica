package etomica.graph.isomorphism;

public final class NodePair {

  private byte first;
  private byte second;
  private static NodePair NODE_PAIR_NULL = new NodePair(SearchState.NULL_NODE, SearchState.NULL_NODE);

  public NodePair(Byte n1, Byte n2) {

    first = n1;
    second = n2;
  }

  public byte getN1() {

    return first;
  }

  public byte getN2() {

    return second;
  }

  public static NodePair nullPair() {

    return NODE_PAIR_NULL;
  }

  @Override
  public String toString() {

    return "(" + first + "," + second + ")";
  }
}