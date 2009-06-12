package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.bitmap.Bitmap;
import etomica.virial.cluster2.graph.EdgesRepresentation;

public abstract class AbstractBitmapRepresentation implements
    EdgesRepresentation {

  private byte nodeCount;
  private Bitmap edges;

  @SuppressWarnings("unused")
  private AbstractBitmapRepresentation() {

    // block instantiation
  }

  public AbstractBitmapRepresentation(byte numNodes) {

    nodeCount = numNodes;
  }

  protected Bitmap getEdges() {

    return edges;
  }

  public byte getNodeCount() {

    return nodeCount;
  }

  public boolean hasEdge(int fromNodeID, int toNodeID) {

    return getEdges().testBit(getEdgeID(fromNodeID, toNodeID));
  }

  public void setEdgesBitmap(Bitmap edgesStore) {

    edges = edgesStore;
  }

  public String toString(int edgeID) {

    return "(" + getFromNodeID(edgeID) + "," + getToNodeID(edgeID) + ")";
  }

  @Override
  public String toString() {

    String result = "";
    boolean first = true;
    for (int i = 0; i < getNodeCount(); i++) {
      for (int j = i + 1; j < getNodeCount(); j++) {
        if (hasEdge(i, j)) {
          if (!first) {
            result += ", ";
          }
          result += toString(getEdgeID(i, j));
          first = false;
        }
      }
    }
    return "<" + result + ">";
  }
}