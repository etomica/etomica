/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Edge;
import etomica.graph.model.Graph;

/**
 * Deletes all edges of a given color.  This is useful for e=f+1.  After using
 * substitution with an explicit 1-bond, you can use this class to delete the
 * resulting 1-bonds.
 * 
 * @author Andrew Schultz
 */
public class DeleteEdge implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {

    assert (params instanceof DeleteEdgeParameters);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      result.add(apply(g, (DeleteEdgeParameters) params));
    }
    return result;
  }

  public Graph apply(Graph g, DeleteEdgeParameters params) {

    Graph newGraph = g.copy();
    for (Edge edge : g.edges()) {
      if (edge.getColor() == params.color()) {
        byte fromNode = g.getFromNode(edge.getId());
        byte toNode = g.getToNode(edge.getId());
        newGraph.deleteEdge(fromNode, toNode);
      }
    }
    return newGraph;
  }
}