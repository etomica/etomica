package etomica.graph.operations;

import etomica.graph.operations.Mul.MulParameters;

public class ConvParameters implements Parameters {

  private byte nodeId;
  private MulParameters mulParameters;

  public ConvParameters(byte nodeId, MulParameters mulParameters) {
    this.nodeId = nodeId;
    this.mulParameters = mulParameters;
  }

  public byte nodeId() {
    return nodeId;
  }
  
  public MulParameters mulParameters() {
    return mulParameters;
  }
}
