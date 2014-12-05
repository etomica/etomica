/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Edge;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.Node;

public class Relabel implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {

    assert (params instanceof RelabelParameters);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      Graph newGraph = apply(g, (RelabelParameters) params);
      result.add(newGraph);
    }
    return result;
  }

  public Graph apply(Graph argument, RelabelParameters params) {

    Node[] nodes = new Node[argument.nodeCount()];
    for (byte nodeId = 0; nodeId < nodes.length; nodeId++) {
      byte newNodeId = params.map(nodeId);
      nodes[newNodeId] = argument.getNode(newNodeId).copy();
      nodes[newNodeId].setColor(argument.getNode(nodeId).getColor());
      nodes[newNodeId].setType(argument.getNode(nodeId).getType());
//      nodes[nodeId] = argument.getNode(params.map(nodeId)).copy();
    }
    // copy the mapped nodes
    Graph result = GraphFactory.createGraph(nodes);
    // copy the coefficient
    result.coefficient().multiply(argument.coefficient());
    result.setNumFactors(argument.factors().length);
    result.addFactors(argument.factors());
    // copy the edges from and to mapped nodes
    for (Edge edge : argument.edges()) {
      byte fromNode = params.map(argument.getFromNode(edge.getId()));
      byte toNode = params.map(argument.getToNode(edge.getId()));
      result.putEdge(fromNode, toNode);
      result.getEdge(fromNode, toNode).setColor(edge.getColor());
    }
    return result;
  }
}
