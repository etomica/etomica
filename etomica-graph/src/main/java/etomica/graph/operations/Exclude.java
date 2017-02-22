/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Edge;
import etomica.graph.model.Graph;
import etomica.graph.model.Node;

/**
 * exclude diagrams having multiple association bond
 */
public class Exclude implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {

    assert (params instanceof ExcludeParameters);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      Set<Graph> newSet = apply(g, (ExcludeParameters) params);
      if (newSet != null) {
        result.addAll(newSet);
      }
    }
    Unary isoFree = new IsoFree();
    return isoFree.apply(result, params);
  }

  public Set<Graph> apply(Graph graph, ExcludeParameters params) {

    assert (params instanceof ExcludeParameters);
    Set<Graph> result = new HashSet<Graph>();

    for (Node node1 : graph.nodes()) {
      for (Node node2 : graph.nodes()) {
        if (node1 == node2 || !graph.hasEdge(node1.getId(), node2.getId())) {
          continue;
        }
        Edge edge12 = graph.getEdge(node1.getId(), node2.getId());
        char[] elementEdges1 = params.bondMap.get(edge12.getColor());
        if (elementEdges1 == null) {
          continue;
        }
        char elementEdge12 = elementEdges1[node2.getId() > node1.getId() ? 0 : 1];
        for (Node node3 : graph.nodes()) {
          if (node1 == node3 || node2 == node3 || node3.getId() < node2.getId()
              || !graph.hasEdge(node1.getId(), node3.getId())) {
            continue;
          }
          Edge edge13 = graph.getEdge(node1.getId(), node3.getId());
          char[] elementEdges2 = params.bondMap.get(edge13.getColor());
          if (elementEdges2 == null) {
            continue;
          }
          char elementEdge13 = elementEdges2[node3.getId() > node1.getId() ? 0 : 1];
          if (elementEdge12 == elementEdge13) {
            return null;
          }
        }
      }
    }
    result.add(graph.copy());
    return result;
  }
}