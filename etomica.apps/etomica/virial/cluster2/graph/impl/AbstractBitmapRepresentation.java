/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.bitmap.Bitmap;
import etomica.virial.cluster2.graph.EdgeVisitor;
import etomica.virial.cluster2.graph.EdgesRepresentation;

/**
 * Representation of an undirected graph.
 * 
 * @author Demian Lessa
 */
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

    return getEdges().testBit(getEdgeID(fromNodeID, toNodeID))
        || getEdges().testBit(getEdgeID(toNodeID, fromNodeID));
  }

  public int getInDegree(int nodeID) {

    int result = 0;
    for (int i = 0; i < getNodeCount(); i++) {
      if (i != nodeID) {
        result += hasEdge(i, nodeID) ? 1 : 0;
      }
    }
    return result;
  }

  public int getInNode(int nodeID, int index) {

    int found = 0;
    for (int i = 0; i < getNodeCount(); i++) {
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
    for (int i = 0; i < getNodeCount(); i++) {
      if (i != nodeID) {
        result += hasEdge(nodeID, i) ? 1 : 0;
      }
    }
    return result;
  }

  public int getOutNode(int nodeID, int index) {

    int found = 0;
    for (int i = 0; i < getNodeCount(); i++) {
      if (i != nodeID) {
        found += hasEdge(nodeID, i) ? 1 : 0;
      }
      if (index == found - 1) {
        return i;
      }
    }
    return -1;
  }
  
  public String getUpperTriangle() {
    
    String result = "";
    for (int n1 = 0; n1 < getNodeCount(); n1++) {
      for (int n2 = n1+1; n2 < getNodeCount(); n2++) {
        result += getEdges().testBit(getEdgeID(n1, n2)) ? '1' : '0';
      }
    }
    return result;
  }

  public void setEdgesBitmap(Bitmap edgesStore) {

    edges = edgesStore;
  }

  public String toString(int edgeID) {

    return "(" + getFromNodeID(edgeID) + "," + getToNodeID(edgeID) + ")";
  }

  @Override
  public boolean equals(Object obj) {
  
    if (obj instanceof EdgesRepresentation) {
      EdgesRepresentation other = (EdgesRepresentation) obj;
      return getUpperTriangle().equals(other.getUpperTriangle());
    }
    return false;
  }
  
  @Override
  public int hashCode() {
  
    return getUpperTriangle().hashCode();
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

  public void visitEdges(EdgeVisitor visitor) {

    for (int i = 0; i < getNodeCount(); i++) {
      for (int j = i + 1; j < getNodeCount(); j++) {
        if (hasEdge(i, j)) {
          if (!visitor.visit(i, j)) {
            return;
          }
        }
      }
    }
  }
}