/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.bitmap.Bitmap;
import etomica.virial.cluster2.graph.EdgesRepresentation;

/*
 * This class maps the edges of an N by N matrix onto a Bitmap using N(N-1)/2 
 * bits. The linearization that maps each edge (n1,n2) such that n1 > n2 onto 
 * the bitmap is illustrated by example below and it has an average cost 
 * proportional to N. Since this class reserves storage for one of the edges
 * (n1,n2) and (n2,n1), it can only be used to represent edges of undirected 
 * graphs.
 * 
 * edgeID (from,to) 
 * ------ ---------
 *     0     (0,1) 
 *     1     (0,2) 
 *     2     (0,3) 
 *     3     (1,2) 
 *     4     (1,3) 
 *     5     (2,3)
 * 
 * @author Demian Lessa
 */
public class UpperTriangleRepresentation extends AbstractBitmapRepresentation {

  public UpperTriangleRepresentation(byte numNodes) {

    super(numNodes);
  }

  public EdgesRepresentation copy() {
    
    AbstractBitmapRepresentation er = new UpperTriangleRepresentation(getNodeCount());
    er.setEdgesBitmap(getEdges().copy());
    return er;
  }
  
  public EdgesRepresentation complement() {

    UpperTriangleRepresentation c = new UpperTriangleRepresentation(getNodeCount());
    Bitmap store = getEdges().copy();
    store.not();
    c.setEdgesBitmap(store);
    return c;
  }

  public int getCapacity() {

    return getNodeCount() * (getNodeCount() - 1) / 2;
  }

  public int getEdgeCount() {

    return getEdges().bitCount();
  }

  public int getEdgeID(int fromNodeID, int toNodeID) {

    assert (fromNodeID != toNodeID);
    if (fromNodeID > toNodeID) {
      return getEdgeID(toNodeID, fromNodeID);
    }
    return (toNodeID - fromNodeID - 1) + sumMaxEdges(0, fromNodeID - 1);
  }

  public int getFromNodeID(int edgeID) {

    int fromNodeID = 0;
    int offset = maxEdges(fromNodeID) - 1;
    while (edgeID > offset) {
      fromNodeID++;
      offset += maxEdges(fromNodeID);
    }
    return fromNodeID;
  }

  public int getToNodeID(int edgeID) {

    int fromNodeID = getFromNodeID(edgeID);
    int fromEdgeID = getEdgeID(fromNodeID, fromNodeID + 1);
    return (fromNodeID + 1) + (edgeID - fromEdgeID);
  }

  public String getUpperTriangle() {
    
    return getEdges().toString();
  }
  
  public boolean hasEdge(int fromNodeID, int toNodeID) {

    return getEdges().testBit(getEdgeID(fromNodeID, toNodeID))
        || getEdges().testBit(getEdgeID(toNodeID, fromNodeID));
  }
  
// @Override
// public String toString() {
//
// String result = "";
// Bitmap printSet = getEdges().copy();
// while (printSet.bitCount() > 0) {
// int bit = printSet.hsb();
// result += toString(bit);
// printSet.clearBit(bit);
// if (printSet.bitCount() > 0) {
// result += ", ";
// }
// }
// return "<" + result + ">";
// }
  /*
   * Returns the largest number of edges encoded with the fromNodeID as the
   * first node in the edge pair (fromNodeID, toNodeID).
   */
  protected int maxEdges(int fromNodeID) {

    return (getNodeCount() - fromNodeID - 1);
  }

  protected int sumMaxEdges(int firstNodeID, int lastNodeID) {

    int result = 0;
    for (int nodeID = firstNodeID; nodeID <= lastNodeID; nodeID++) {
      result += maxEdges(nodeID);
    }
    return result;
  }
}