package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.graph.EdgeAttributes;
import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.EdgesRepresentation;
import etomica.virial.cluster2.graph.EdgesMetadata;
import etomica.virial.cluster2.graph.GraphFactory;

public class SimpleEdges implements Edges {

  private EdgesMetadata metadata;
  private EdgesRepresentation representation;

  public SimpleEdges(EdgesRepresentation representation, EdgesMetadata metadata) {

    this.representation = representation;
    this.metadata = metadata;
  }

  public int count() {

    return getRepresentation().getEdgeCount();
  }

  protected EdgesRepresentation getRepresentation() {

    return representation;
  }

  public EdgesMetadata getMetadata() {

    return metadata;
  }

  public boolean hasEdge(int fromNodeID, int toNodeID) {

    return getRepresentation().hasEdge(fromNodeID, toNodeID);
  }

  @Override
  public String toString() {

    return getRepresentation().toString();
  }

  public EdgeAttributes getAttributes(int fromNodeID, int toNodeID) {

    return GraphFactory.defaultEdgeAttributes();
  }

  public int getInDegree(int nodeID) {

    int result = 0;
    for (int i = 0; i < getRepresentation().getNodeCount(); i++) {
      if (i != nodeID) {
        result += hasEdge(i, nodeID) ? 1 : 0;
      }
    }
    return result;
  }

  public int getInNode(int nodeID, int index) {

    int found = 0;
    for (int i = 0; i < getRepresentation().getNodeCount(); i++) {
      if (i != nodeID) {
        found += hasEdge(i, nodeID) ? 1 : 0;
      }
      if (index == found - 1) {
        return i;
      }
    }
    return -1;
  }

  public int getOutDegree(int nodeID) {

    int result = 0;
    for (int i = 0; i < getRepresentation().getNodeCount(); i++) {
      if (i != nodeID) {
        result += hasEdge(nodeID, i) ? 1 : 0;
      }
    }
    return result;
  }

  public int getOutNode(int nodeID, int index) {

    int found = 0;
    for (int i = 0; i < getRepresentation().getNodeCount(); i++) {
      if (i != nodeID) {
        found += hasEdge(nodeID, i) ? 1 : 0;
      }
      if (index == found - 1) {
        return i;
      }
    }
    return -1;
  }
}