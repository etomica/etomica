package etomica.graph.operations;

public class ConvParameters implements Parameters {

  private byte nodeId;

  public ConvParameters(byte nodeId) {

    this.nodeId = nodeId;
  }

  public byte nodeId() {

    return nodeId;
  }
}
