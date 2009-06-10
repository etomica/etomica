package etomica.virial.cluster2.graph.isomorphism;

public final class NodePair {

  private int firstNode;
  private int secondNode;

  public NodePair(Integer n1, Integer n2) {

    firstNode = n1;
    secondNode = n2;
  }

  public int getFirstNode() {

    return firstNode;
  }

  public int getSecondNode() {

    return secondNode;
  }
}