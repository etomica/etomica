/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.List;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;

/**
 * This class takes the given nodes from the given graph and returns a new
 * graph consisting only of those nodes (and the bonds between them).
 */
public class GraphExtract {

  public Graph apply(Graph argument, List<Byte> nodeIDs) {
    byte n = (byte)nodeIDs.size();
    Graph g = GraphFactory.createGraph(n);
    for (byte i=0; i<n; i++) {
      byte ii = nodeIDs.get(i);
      g.getNode(i).setColor(argument.getNode(ii).getColor());
      g.getNode(i).setType(argument.getNode(ii).getType());
      for (byte j=0; j<i; j++) {
        byte jj = nodeIDs.get(j);
        if (argument.hasEdge(jj,ii)) {
          g.putEdge(j,i);
          g.getEdge(j,i).setColor(argument.getEdge(jj,ii).getColor());
        }
      }
    }
    return g;
  }
}
