/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.property;

import java.util.HashSet;
import java.util.Set;

import static etomica.graph.model.Metadata.*;
import static etomica.graph.traversal.Traversal.*;


import etomica.graph.model.Graph;
import etomica.graph.traversal.Biconnected;
import etomica.graph.traversal.TraversalVisitor;


public class HasArticulationPoint implements Property {

  public boolean check(Graph graph) {

    if (graph == null) {
      return false;
    }
    Set<Byte> rootNodes = new HashSet<Byte>();
    for (byte nodeID = 0; nodeID < graph.nodeCount(); nodeID++) {
      if (graph.nodes().get(nodeID).getType() == TYPE_NODE_ROOT) {
        rootNodes.add(nodeID);
      }
    }
    // invoke a BC traversal starting and return true IFF the visitor
    // detected that the graph has an articulation point
    APVisitor v = new APVisitor(graph, rootNodes);
    new Biconnected().traverseAll(graph, v);
    return v.hasArticulationPoint();
  }
}

class APVisitor implements TraversalVisitor {

  private Graph graph;
  private boolean isArticulated = false;
  private boolean isArticulation = false;
  private Set<Byte> rootNodes = new HashSet<Byte>();
  private boolean isConnected = false;
  private Set<Byte> visitedNodes = new HashSet<Byte>();
  private Set<Byte> bicomponent;
  private Set<Set<Byte>> bicomponents = new HashSet<Set<Byte>>();
  private Set<Byte> articulationPoints = new HashSet<Byte>();

  public APVisitor(Graph graph, Set<Byte> rootNodes) {

    this.graph = graph;
    this.rootNodes = rootNodes;
  }

  public boolean visit(byte nodeID, byte status) {

    if (status == STATUS_START_BICOMPONENT) {
      bicomponent = new HashSet<Byte>();
    }
    if (status == STATUS_VISITED_BICOMPONENT) {
      bicomponents.add(bicomponent);
      bicomponent = null;
    }
    else if (status == STATUS_VISITED_COMPONENT) {
      isConnected = isConnected || (visitedNodes.size() == graph.nodeCount());
      visitedNodes.clear();
    }
    // the next node is an articulation point and should not be processed
    else if (status == STATUS_ARTICULATION_POINT) {
      isArticulation = true;
      isArticulated = true;
    }
    // visiting a node in the current biconnected component
    else if (status == STATUS_VISITED_NODE) {
      // if it is an articulation point, ignore it
      if (isArticulation) {
        isArticulation = !isArticulation;
        articulationPoints.add(nodeID);
      }
      else {
        visitedNodes.add(nodeID);
      }
    }
    return true;
  }

  public boolean hasArticulationPoint() {

    // this graph is not connected
    if (!isConnected) {
      return true;
    }
    // this graph has no (graph-theoretical) articulation points
    if (!isArticulated) {
      return false;
    }
    // if we don't have at least 2 root nodes, then we just care about the graph
    if (rootNodes.size() < 2) {
      return isArticulated;
    }
    // some bicomponent has no root node, or has a single root node that is an
    // articulation point
    for (Set<Byte> bc : bicomponents) {
      int rootPoints = 0;
      int singletonRootNode = -1;
      for (Byte node : bc) {
        if (rootNodes.contains(node)) {
          rootPoints++;
          singletonRootNode = node;
        }
      }
      // no root point is found in this bicomponent
      if (rootPoints == 0) {
        return true;
      }
      // the singleton root point in this bicomponent is an articulation point
      else if ((rootPoints == 1) && (articulationPoints.contains(singletonRootNode))) {
        return true;
      }
    }
    // if we have not found bicomponents matching the articulation point
    // criteria, it is because this graph has no articulation points
    return false;
  }
}
