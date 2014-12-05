/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph;

/*
 * This interface provides a simple mechanism to map pairs of nodes onto a
 * linear scale of edge identifiers. Hence, it provides a means to access
 * edge information from some particular storage that is indexed on a single
 * integer, such as a vector or a bitmap.
 * 
 * @author Demian Lessa
 * 
 */
public interface EdgesRepresentation {

  public EdgesRepresentation complement();

  public EdgesRepresentation copy();

  public int getCapacity();

  public int getEdgeCount();

  public int getEdgeID(int fromNodeID, int toNodeID);

  public int getFromNodeID(int edgeID);

  public byte getNodeCount();

  public int getToNodeID(int edgeID);

  public int getInDegree(int nodeID);

  public int getInNode(int nodeID, int index);

  public int getOutDegree(int nodeID);

  public int getOutNode(int nodeID, int index);

  public boolean hasEdge(int fromNodeID, int toNodeID);

  public String toString(int edgeID);

  public void visitEdges(EdgeVisitor visitor);

  public String getUpperTriangle();
}