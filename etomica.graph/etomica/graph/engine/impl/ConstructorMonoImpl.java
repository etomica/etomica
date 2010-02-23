package etomica.graph.engine.impl;

import etomica.graph.engine.Parser.ConstructorMono;

public class ConstructorMonoImpl extends ConstructorImpl implements ConstructorMono {

  private byte fieldNodes;
  private boolean isoFree;
  private byte rootNodes;

  public ConstructorMonoImpl(byte rootNodes, byte fieldNodes, boolean isoFree) {

    this.rootNodes = rootNodes;
    this.fieldNodes = fieldNodes;
    this.isoFree = isoFree;
  }

  public byte getFieldNodes() {

    return fieldNodes;
  }

  public byte getRootNodes() {

    return rootNodes;
  }

  public boolean isIsoFree() {

    return isoFree;
  }
}