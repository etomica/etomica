/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.traversal;

import static etomica.graph.traversal.Traversal.STATUS_ARTICULATION_POINT;
import static etomica.graph.traversal.Traversal.STATUS_START_BICOMPONENT;
import static etomica.graph.traversal.Traversal.STATUS_VISITED_BICOMPONENT;
import static etomica.graph.traversal.Traversal.STATUS_VISITED_NODE;

import java.util.ArrayList;
import java.util.List;

import etomica.graph.model.Graph;

/**
 * This TraversalVisitor visits nodes and keeps track of biconnected
 * components.
 */
public class BCVisitor implements TraversalVisitor {

  private List<List<Byte>> biComponents;
  private List<Byte> biComponent;
  private boolean isArticulation = false;

  public BCVisitor(List<List<Byte>> biComponents) {
    this.biComponents = biComponents;
  }

  public boolean visit(byte nodeID, byte status) {

    // the next node is an articulation point and should not be processed
    if (status == STATUS_START_BICOMPONENT) {
      biComponent = new ArrayList<Byte>();
      biComponents.add(biComponent);
    }
    else if (status == STATUS_VISITED_BICOMPONENT) {
      biComponent = null;
    }
    else if  (status == STATUS_ARTICULATION_POINT) {
      isArticulation = true;
    }
    // visiting a node in the current biconnected component
    else if (status == STATUS_VISITED_NODE) {
      // if it is an articulation point, ignore it
      if (isArticulation) {
        isArticulation = false;
      }
      else {
        biComponent.add(nodeID);
      }
    }
    return true;
  }

  /**
   * Returns the biconnected components of the given graph.
   */
  public static List<List<Byte>> getBiComponents(Graph g) {
    List<List<Byte>> biComponents = new ArrayList<List<Byte>>();
    BCVisitor v = new BCVisitor(biComponents);
    new Biconnected().traverseAll(g, v);
    return biComponents;
  }
}