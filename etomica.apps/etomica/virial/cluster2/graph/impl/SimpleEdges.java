package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.bitmap.Bitmap;
import etomica.virial.cluster2.graph.EdgeAttributes;
import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.EdgesRepresentation;
import etomica.virial.cluster2.graph.EdgesMetadata;
import etomica.virial.cluster2.graph.GraphFactory;

public class SimpleEdges implements Edges {

  private Bitmap edges = null;
  private EdgesMetadata metadata;
  private EdgesRepresentation representation;

  public SimpleEdges(Bitmap edges, EdgesRepresentation representation,
      EdgesMetadata metadata) {

    this.edges = edges;
    this.representation = representation;
    this.metadata = metadata;
  }

  @Override
  public int count() {

    return getEdges().bitCount();
  }

  protected EdgesRepresentation getRepresentation() {

    return representation;
  }

  @Override
  public EdgesMetadata getMetadata() {

    return metadata;
  }

  @Override
  public boolean hasEdge(int fromNodeID, int toNodeID) {

    return getEdges().testBit(getRepresentation().getEdgeID(fromNodeID, toNodeID));
  }

  @Override
  public String toString() {

    String result = "";
    Bitmap printSet = getEdges().copy();
    while (printSet.bitCount() > 0) {
      int bit = printSet.hsb();
      result += getRepresentation().toString(bit);
      printSet.clearBit(bit);
      if (printSet.bitCount() > 0) {
        result += ", ";
      }
    }
    return "<" + result + ">";
  }

  @Override
  public EdgeAttributes getAttributes(int fromNodeID, int toNodeID) {

    return GraphFactory.defaultEdgeAttributes();
  }

  protected Bitmap getEdges() {

    return edges;
  }

  @Override
  public int getInDegree(int nodeID) {

    int result = 0;
    for (int i = 0; i < getRepresentation().getNodeCount(); i++) {
      if (i != nodeID) {
        result += getEdges().testBit(getRepresentation().getEdgeID(i, nodeID)) ? 1 : 0;
      }
    }
    return result;
  }

  @Override
  public int getInNode(int nodeID, int index) {

    int found = 0;
    for (int i = 0; i < getRepresentation().getNodeCount(); i++) {
      if (i != nodeID) {
        found += getEdges().testBit(getRepresentation().getEdgeID(i, nodeID)) ? 1 : 0;
      }
      if (index == found - 1) {
        return i;
      }
    }
    return -1;
  }

  @Override
  public int getOutDegree(int nodeID) {

    int result = 0;
    for (int i = 0; i < getRepresentation().getNodeCount(); i++) {
      if (i != nodeID) {
        result += getEdges().testBit(getRepresentation().getEdgeID(nodeID, i)) ? 1 : 0;
      }
    }
    return result;
  }

  @Override
  public int getOutNode(int nodeID, int index) {

    int found = 0;
    for (int i = 0; i < getRepresentation().getNodeCount(); i++) {
      if (i != nodeID) {
        found += getEdges().testBit(getRepresentation().getEdgeID(nodeID, i)) ? 1 : 0;
      }
      if (index == found - 1) {
        return i;
      }
    }
    return -1;
  }
}