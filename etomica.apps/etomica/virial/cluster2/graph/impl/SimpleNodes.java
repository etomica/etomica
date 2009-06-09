package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.graph.Nodes;

public class SimpleNodes implements Nodes {

  private byte nodeCount = 0;

  public SimpleNodes(byte numNodes) {

    nodeCount = numNodes;
  }

  @Override
  public byte count() {

    return nodeCount;
  }

  @Override
  public byte count(byte color) {

    return (color == Nodes.NODE_COLOR_DEFAULT ? count() : 0);
  }

  @Override
  public boolean hasColoredNode(byte color) {

    return (color == Nodes.NODE_COLOR_DEFAULT);
  }
}