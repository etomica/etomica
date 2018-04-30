/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;


import java.util.List;
import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.GraphList;
import etomica.graph.property.IsConnected;
import etomica.graph.traversal.CVisitor;

/**
 * Split a disconnected graph into multiple graphs.
 * 
 * @author Andrew Schultz
 */
public class SplitGraph {

  protected IsConnected isCon = new IsConnected();
    
  public Set<Graph> apply(Graph g) {
    Set<Graph> result = new GraphList(null); // we want to return these in order
    if (isCon.check(g)) {
      result.add(g.copy());
      return result;
    }
    List<List<Byte>> components = CVisitor.getComponents(g);

    for (int iComp = 0; iComp<components.size(); iComp++) {
      List<Byte> nodes = components.get(iComp);
      Graph newG = GraphFactory.createGraph((byte)nodes.size());
      result.add(newG);
      for (byte id = 0; id<nodes.size(); id++) {
        byte nodeID = nodes.get(id);
        newG.getNode(id).setColor(g.getNode(nodeID).getColor());
        newG.getNode(id).setType(g.getNode(nodeID).getType());
        for (byte id2 = (byte)(id+1); id2<nodes.size(); id2++) {
          byte nodeID2 = nodes.get(id2);
          if (g.hasEdge(nodeID, nodeID2)) {
            newG.putEdge(id, id2);
            newG.getEdge(id, id2).setColor(g.getEdge(nodeID, nodeID2).getColor());
          }
        }
      }
      if (iComp == 0) {
        newG.coefficient().multiply(g.coefficient());
      }
    } 
    return result;
  }
}
