/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.ArrayList;
import java.util.List;

import etomica.graph.model.Graph;
import etomica.graph.model.Metadata;
import etomica.graph.model.Node;
import etomica.graph.traversal.CVisitor;

/**
 * GraphOp class that adjusts root points in each component so that
 * + the number of root points is maintained
 * + no component has more than one root point
 * + all root points are included in all the "earliest" components
 * + if a component has a root point, it is the first node in that component
 *
 * This class only makes sense if the specific identity of each root
 * point is unimportant (Metadata.ROOT_POINTS_SPECIAL=false).
 *
 * @author Andrew Schultz
 */
public class GraphOpMaxRoot implements GraphOp {
  public Graph apply(Graph g) {
    List<Byte> roots = new ArrayList<Byte>();
    for (Node node : g.nodes()) {
      if (node.getType() == Metadata.TYPE_NODE_ROOT) {
        roots.add(node.getId());
      }
    }
    if (roots.size() == 0) return g;
    List<List<Byte>> comps = CVisitor.getComponents(g);
    if (roots.size() > comps.size()) {
      // more root points than components.  perhaps some single-point root
      // nodes could be combined (not our problem!).
      throw new RuntimeException("Something is terribly wrong!");
    }
    int iRoot = 0;
    g = g.copy();
    // we'll go out on a limb here and assume that the components are in order
    // (for some definition of "in order") and that the nodes within the
    // components are in order (with the same caveat).  The second is actually
    // probably not true, but it is still true that the first node listed is
    // the first node; the rest are probably out of order, but we don't care.
    for (List<Byte> comp : comps) {
      if (iRoot < roots.size()) {
        g.getNode(comp.get(0)).setType(Metadata.TYPE_NODE_ROOT);
        for (byte j=1; j<comp.size(); j++) {
          g.getNode(comp.get(j)).setType(Metadata.TYPE_NODE_FIELD);
        }
      }
      else {
        for (byte j=0; j<comp.size(); j++) {
          g.getNode(comp.get(j)).setType(Metadata.TYPE_NODE_FIELD);
        }
      }
      iRoot++;
    }
    return g;
  }
}