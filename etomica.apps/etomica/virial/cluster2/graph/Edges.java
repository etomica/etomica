/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph;

import java.util.List;
import java.util.Set;

public interface Edges {

  public static final char EDGE_COLOR_EMPTY = '$';
  public static final char EDGE_COLOR_0 = 'a';
  public static final char EDGE_COLOR_1 = 'b';
  public static final char EDGE_COLOR_2 = 'c';
  public static final char EDGE_COLOR_3 = 'd';
  public static final char EDGE_COLOR_4 = 'e';
  public static final char EDGE_COLOR_5 = 'f';
  public static final char EDGE_COLOR_6 = 'g';
  public static final char EDGE_COLOR_7 = 'h';
  public static final char EDGE_COLOR_8 = 'i';
  public static final char EDGE_COLOR_9 = 'j';
  public static final char EDGE_COLOR_DEFAULT = EDGE_COLOR_0;

  // returns the set of edges corresponding to the complement of this set of edges
  public Edges complement();

  // returns an independent copy of this set of edges
  public Edges copy();

  // returns an independent copy of this set of edges with the inverse of the coefficient
  public Edges ncopy();

  // returns the number of edges
  public int count();

  // get the attributes of the edge (node1, node2)
  public EdgeAttributes getAttributes(int fromNodeID, int toNodeID);

  // returns the number of head nodes adjacent to nodeID
  public int getInDegree(int nodeID);

  // returns the i-th head node adjacent to nodeID
  public int getInNode(int nodeID, int index);

  public EdgesMetadata getMetadata();

  // returns the number of tail nodes adjacent to nodeID
  public int getOutDegree(int nodeID);

  // returns the i-th tail node adjacent to nodeID
  public int getOutNode(int nodeID, int index);

  public boolean hasEdge(int fromNodeID, int toNodeID);

  // returns the object with the internal representation of the edges
  public EdgesRepresentation getRepresentation();
  
  // a set with the distinct edge colors
  public Set<Character> getColors();

  // a list with all edge IDs of a given color
  public List<Integer> getPartition(char edgeColor);
}