/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.util.GraphNumber;

public class DropOrphanNodes {

  public Graph apply(Graph g) {
    List<Byte> orphanNodes = new ArrayList<Byte>();
    for (byte i=0; i<g.nodeCount(); i++) {
      if (g.getOutDegree(i) == 0) {
        orphanNodes.add(i);
      }
    }
    Map<Byte,Byte> nodeMap = new HashMap<Byte,Byte>();
    byte newNodes = 0;
    for (byte i=0; i<g.nodeCount(); i++) {
      if (orphanNodes.contains(i)) continue;
      nodeMap.put(i, newNodes);
      newNodes++;
    }
    Graph result = GraphFactory.createGraph(newNodes);
    for (byte i=0; i<g.nodeCount(); i++) {
      if (orphanNodes.contains(i)) continue;
      byte newNode = nodeMap.get(i);
      result.getNode(newNode).setColor(g.getNode(i).getColor());
      result.getNode(newNode).setType(g.getNode(i).getType());
      for (byte j=0; j<g.getOutDegree(i); j++) {
        byte jNode = g.getOutNode(i, j);
        byte newJ = nodeMap.get(jNode);
        result.putEdge(newNode, newJ);
        result.getEdge(newNode, newJ).setColor(g.getEdge(i,jNode).getColor());
        result.getEdge(newNode, newJ).setType(g.getEdge(i,jNode).getType());
      }
    }
    return result;
  }
  
  public static void main(String[] args) {
    Graph g = GraphNumber.makeGraph(12);
    System.out.println(g);
    DropOrphanNodes dropper = new DropOrphanNodes();
    System.out.println(dropper.apply(g));
  }
}
