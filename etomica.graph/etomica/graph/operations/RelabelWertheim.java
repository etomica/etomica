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
import etomica.graph.model.impl.MetadataImpl;

public class RelabelWertheim implements Unary {

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
      nodes[nodeId] = argument.getNode(nodeId).copy();
      nodes[nodeId].setColor(argument.getNode(params.map(nodeId)).getColor());
      nodes[nodeId].setType(argument.getNode(params.map(nodeId)).getType());
    }
    // copy the mapped nodes
    Graph result = GraphFactory.createGraph(nodes);
    // copy the coefficient
    result.coefficient().multiply(argument.coefficient());
    result.setNumFactors(argument.factors().length);
    result.addFactors(argument.factors());
    // copy the edges from and to mapped nodes
    for (Edge edge : argument.edges()) {
      byte oldFromNode = argument.getFromNode(edge.getId());
      byte oldToNode = argument.getToNode(edge.getId());
      byte fromNode = params.map(oldFromNode);
      byte toNode = params.map(oldToNode);
      result.putEdge(fromNode, toNode);
      char oldColor = edge.getColor();
      if ((oldFromNode-oldToNode)*(fromNode-toNode) < 0){
    	  oldColor = MetadataImpl.getReverseEdgeColor(oldColor);
      }
      result.getEdge(fromNode, toNode).setColor(oldColor);
    }
    return result;
  }
}
