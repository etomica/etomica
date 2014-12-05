/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph.impl;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import etomica.virial.cluster2.graph.EdgeAttributes;
import etomica.virial.cluster2.graph.EdgeVisitor;
import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.EdgesMetadata;
import etomica.virial.cluster2.graph.EdgesRepresentation;
import etomica.virial.cluster2.graph.GraphFactory;

public class SimpleEdges implements Edges {

  private EdgesMetadata metadata;
  private EdgesRepresentation representation;

  public SimpleEdges(EdgesRepresentation rep, EdgesMetadata meta) {

    representation = rep;
    metadata = meta;
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

    if (getRepresentation().hasEdge(fromNodeID, toNodeID)) {
      return SimpleEdgeAttributes.getAttributes(Edges.EDGE_COLOR_DEFAULT);
    }
    else {
      return SimpleEdgeAttributes.getAttributes(Edges.EDGE_COLOR_EMPTY);
    }
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

  public Set<Character> getColors() {

    ColorVisitor v = new ColorVisitor(this);
    representation.visitEdges(v);
    return v.getColors();
  }

  public List<Integer> getPartition(char edgeColor) {

    ColorPartitionVisitor v = new ColorPartitionVisitor(this, edgeColor);
    representation.visitEdges(v);
    return v.getPartition();
  }
}

class ColorVisitor implements EdgeVisitor {

  Edges owner;
  Set<Character> colors = new HashSet<Character>();

  public ColorVisitor(Edges owner) {
    
    this.owner = owner;
  }

  public Set<Character> getColors() {
    return colors;
  }
  
  public boolean visit(int node1, int node2) {

    colors.add(owner.getAttributes(node1, node2).getColor());
    return true;
  }
  
}

class ColorPartitionVisitor implements EdgeVisitor {

  Edges owner;
  List<Integer> partition = new ArrayList<Integer>();
  private char color;

  public ColorPartitionVisitor(Edges owner, char color) {
    
    this.owner = owner;
    this.color = color;
  }

  public List<Integer> getPartition() {
    return partition;
  }
  
  public boolean visit(int node1, int node2) {

    if (owner.getAttributes(node1, node2).getColor() == color) {
      partition.add(owner.getRepresentation().getEdgeID(node1, node2));
    }
    return true;
  } 
}