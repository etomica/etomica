/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.traversal;

import etomica.graph.model.Graph;

/**
 * This interface generalizes algorithms that traverse the nodes of a graph. The traversal
 * allows a visitor to visit each of the nodes as they are traversed. If the graph is not
 * connected, only one connected component is traversed.
 */
public interface Traversal {

  public static byte ERROR = 96;
  public static byte ERROR_TRAVERSAL_EDGES = 125;
  public static byte ERROR_TRAVERSAL_NODES = 126;
  public static byte ERROR_TRAVERSAL_ROOT = 127;

  public static byte NODE_NULL = (byte) 0xFF;

  public static byte STATUS_VISITED_NODE = 0;
  // public static byte VISITED_NONE = 1;
  public static byte STATUS_START_COMPONENT = 1;
  public static byte STATUS_VISITED_COMPONENT = 2;
  public static byte STATUS_START_BICOMPONENT = 4;
  public static byte STATUS_VISITED_BICOMPONENT = 8;
  public static byte STATUS_VISITED_ALL = 16;
  public static byte STATUS_ARTICULATION_POINT = 32;

  /**
   * Traverses all components of the graph, starting at an arbitrary node for each
   * component, and returns the total number of components traversed.
   */
  public byte traverseAll(Graph graph, TraversalVisitor visitor);

  /**
   * Traverses a single component of the graph starting at the designated node.
   *
   * @return true if all nodes in the graph were seen during traversal
   */
  public boolean traverseComponent(byte nodeID, Graph graph, TraversalVisitor visitor);
}