/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.bitmap.Bitmap;
import etomica.virial.cluster2.bitmap.BitmapFactory;
import etomica.virial.cluster2.graph.EdgeVisitor;
import etomica.virial.cluster2.graph.EdgesRepresentation;
import etomica.virial.cluster2.graph.EdgesRepresentationFactory;

/**
 * This class maps every edge of an N by N matrix onto a Bitmap using N^2 bits.
 * The edge (n1,n2) exists in the matrix iff the bit (n1*N + n2) is set. Since
 * this class assumes storage for both edges (n1,n2) and (n2,n1), it can be used
 * to represent edges for both directed and undirected graphs.
 * 
 * @author Demian Lessa
 */
public class AdjacencyMatrixRepresentation extends AbstractBitmapRepresentation {

  public AdjacencyMatrixRepresentation(byte numNodes) {

    super(numNodes);
  }

  public EdgesRepresentation copy() {
    
    AbstractBitmapRepresentation er = new AdjacencyMatrixRepresentation(getNodeCount());
    er.setEdgesBitmap(getEdges().copy());
    return er;
  }

  public EdgesRepresentation complement() {

    AdjacencyMatrixRepresentation c = new AdjacencyMatrixRepresentation(getNodeCount());
    Bitmap store = getEdges().copy();
    store.not();
    for (int i = 0; i < getNodeCount(); i++) {
      store.clearBit(getEdgeID(i,i));
    }
    c.setEdgesBitmap(store);
    return c;
  }

  public int getCapacity() {

    return getNodeCount() * getNodeCount();
  }

  public int getEdgeCount() {

    return getEdges().bitCount() / 2;
  }

  public int getEdgeID(int fromNodeID, int toNodeID) {

    return fromNodeID * getNodeCount() + toNodeID;
  }

  public int getFromNodeID(int edgeID) {

    return edgeID / getNodeCount();
  }

  public int getToNodeID(int edgeID) {

    return edgeID % getNodeCount();
  }

  public void setEdgesBitmap(Bitmap edgesStore) {

    EdgesRepresentationFactory utf = EdgesRepresentationFactory.getFactory(
        true, getNodeCount());
    EdgesRepresentation ut = utf.getRepresentation(edgesStore);
    Bitmap newStore = BitmapFactory.getBitmap(getCapacity(), false);
    ut.visitEdges(new AdjacencyBuilder(newStore, this));
    super.setEdgesBitmap(newStore);
  }
}

class AdjacencyBuilder implements EdgeVisitor {

  private Bitmap store;
  private EdgesRepresentation rep;

  public AdjacencyBuilder(Bitmap store, EdgesRepresentation rep) {

    this.store = store;
    this.rep = rep;
  }

  public boolean visit(int node1, int node2) {

    store.setBit(rep.getEdgeID(node1, node2));
    store.setBit(rep.getEdgeID(node2, node1));
    return true;
  }
}