package etomica.graph.operations;

public class IntParameters implements Parameters {

  private byte nodeId;

  public IntParameters(byte nodeId) {

    this.nodeId = nodeId;
  }

  public byte nodeId() {

    return nodeId;
  }
}
