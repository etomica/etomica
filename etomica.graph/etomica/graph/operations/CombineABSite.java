/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;


import java.util.HashSet;
import java.util.List;
import java.util.Set;

import etomica.graph.model.Edge;
import etomica.graph.model.Graph;
import etomica.graph.model.impl.MetadataImpl;

/**
 * Operation that attempts to change bonds in a graph so that equivalent graphs
 * in a set can be combined.  This class switch bonds corresponding to A and B
 * site with a preference for A.  A and B are considered distinguishable but 
 * equivalent. 
 */
public class CombineABSite implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {
    assert (params instanceof CombineABSiteParameters);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      Graph newGraph = apply(g, (CombineABSiteParameters)params);
      result.add(newGraph);
    }
    return result;
  }

  public Graph apply(Graph g, CombineABSiteParameters params) {
    g = g.copy();
    for (byte j = 0;j<g.nodeCount();j++){//looping over nodes
      int bondedACount = 0;//from 0 to ..
      int bondedBCount = 0;//from 0 to ..
      boolean siteAFirst = false;
      for (byte i=0; i <g.nodeCount();i++){
        if(i==j||!g.hasEdge(j, i))continue;
        Edge edge = g.getEdge(j, i); 
        if (params.bondsA.contains(edge.getColor())) {
          if(bondedBCount == 0){
            siteAFirst = true;
          }
          bondedACount++;
        }
        else if (params.bondsB.contains(edge.getColor())) {
          bondedBCount++;
        }
      }
      if((bondedBCount>bondedACount)||(bondedACount>0&&bondedACount==bondedBCount&&!siteAFirst)){
        for(byte i=(byte)0; i <g.nodeCount();i++){
          if(i==j||!g.hasEdge(j, i))continue;
          Edge edge = g.getEdge(j, i);
          Edge reverseEdge = g.getEdge(i, j);
          char ec = edge.getColor();
          char newC = ' ';
          int index = params.bondsA.indexOf(ec);
          if (index > -1) {
            // edge is an "A" edge.  turn it into the equivalent "B"
            newC = params.bondsB.get(index);
          }
          else {
            index = params.bondsB.indexOf(ec);
            if (index > -1) {
              // edge is an "B" edge.  turn it into the equivalent "A"
              newC = params.bondsA.get(index);
            }
            else {
              continue;
            }
          }
          edge.setColor(newC);
          reverseEdge.setColor(MetadataImpl.getReverseEdgeColor(newC));
        }
      }
    }
    return g;
  }

  public static class CombineABSiteParameters implements Parameters {
    protected final List<Character> bondsA, bondsB;
    public CombineABSiteParameters(List<Character> bonds1, List<Character> bonds2) {
      this.bondsA = bonds1;
      this.bondsB = bonds2;
    }
  }
}
