package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.graph.EdgeAttributes;
import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.EdgesRepresentation;
import etomica.virial.cluster2.graph.EdgesMetadata;
import etomica.virial.cluster2.graph.GraphFactory;

public class SimpleEdges implements Edges {

  private EdgesMetadata metadata;
  private EdgesRepresentation representation;

  public SimpleEdges(EdgesRepresentation rep, EdgesMetadata meta) {

    representation = rep;
    metadata = meta;
  }

  public Edges canonical() {

    return GraphFactory.canonicalEdges(this);
  }

  public Edges copy() {
    
    assert(representation instanceof AbstractBitmapRepresentation);
    AbstractBitmapRepresentation rep = (AbstractBitmapRepresentation) representation;
    return new SimpleEdges(rep.copy(), metadata.copy());
  }
  
  public Edges ncopy() {
    
    assert(representation instanceof AbstractBitmapRepresentation);
    AbstractBitmapRepresentation rep = (AbstractBitmapRepresentation) representation;
    return new SimpleEdges(rep.copy(), metadata.ncopy());
  }
  
  public Edges complement() {

    return GraphFactory.complementEdges(this);
  }

  public int count() {

    return getRepresentation().getEdgeCount();
  }

  public EdgesRepresentation getRepresentation() {

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

  @Override
  public int hashCode() {

    return getRepresentation().hashCode();
  }

  @Override
  public boolean equals(Object obj) {

    if (obj instanceof Edges) {
      Edges other = (Edges) obj;
      return getRepresentation().equals(other.getRepresentation());
    }
    return false;
  }

  public EdgeAttributes getAttributes(int fromNodeID, int toNodeID) {

    return GraphFactory.defaultEdgeAttributes();
  }

  public int getInDegree(int nodeID) {

    return getRepresentation().getInDegree(nodeID);
  }

  public int getInNode(int nodeID, int index) {

    return getRepresentation().getInNode(nodeID, index);
  }

  public int getOutDegree(int nodeID) {

    return getRepresentation().getOutDegree(nodeID);
  }

  public int getOutNode(int nodeID, int index) {

    return getRepresentation().getOutNode(nodeID, index);
  }
}