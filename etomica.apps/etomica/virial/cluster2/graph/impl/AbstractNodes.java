package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.graph.Nodes;


public abstract class AbstractNodes implements Nodes {

  private byte nodeCount = 0;

  public AbstractNodes(byte numNodes) {

    nodeCount = numNodes;
  }

  @Override
  public byte count() {

    return nodeCount;
  }
}