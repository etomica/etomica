/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.property;

import static etomica.graph.traversal.Traversal.STATUS_ARTICULATION_POINT;
import static etomica.graph.traversal.Traversal.STATUS_VISITED_COMPONENT;
import static etomica.graph.traversal.Traversal.STATUS_VISITED_NODE;

import java.util.ArrayList;
import java.util.List;

import etomica.graph.model.Graph;
import etomica.graph.traversal.Biconnected;
import etomica.graph.traversal.TraversalVisitor;

/**
 * Checks for "simple" articulation points -- ignoring root vs. field nodes.
 * Articulation points are points which divide any component into two or more
 * components.
 * 
 * This class also keeps track of the articulation points seen and checks for
 * connectedness.
 * 
 * @author Andrew Schultz
 */
public class HasSimpleArticulationPoint implements Property {

  public boolean check(Graph graph) {

    if (graph == null || graph.nodeCount() < 2) {
      isConnected = true;
      return false;
    }
    // invoke a BC traversal starting and return true IFF the visitor
    // detected that the graph has an articulation point
    articulationPoints.clear();
    SAPVisitor v = new SAPVisitor(articulationPoints);
    new Biconnected().traverseAll(graph, v);
    isConnected = v.isConnected();
    return v.hasArticulationPoint();
  }

  public List<Byte> getArticulationPoints() {
    return articulationPoints;
  }
  
  public boolean isConnected() {
    return isConnected;
  }

  protected final List<Byte> articulationPoints = new ArrayList<Byte>();
  protected boolean isConnected;
}

class SAPVisitor implements TraversalVisitor {

  private boolean isArticulated = false;
  private boolean isArticulation = false;
  private List<Byte> articulationPoints;
  private int components;

  public SAPVisitor(List<Byte> articulationPoints) {
    this.articulationPoints = articulationPoints;
  }

  public boolean visit(byte nodeID, byte status) {

    if (status == STATUS_VISITED_COMPONENT) {
      components++;
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
    }
    return true;
  }

  public boolean isConnected() {
    return components == 1;
  }
  
  public boolean hasArticulationPoint() {
    return isArticulated;
  }
}
