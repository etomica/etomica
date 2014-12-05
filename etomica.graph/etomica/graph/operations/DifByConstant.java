/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Edge;
import etomica.graph.model.Graph;

public class DifByConstant implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {

    assert (params instanceof DifByConstantParameters);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      Set<Graph> newSet = apply(g, (DifByConstantParameters) params);
      if (newSet != null) {
        result.addAll(newSet);
      }
    }
    Unary isoFree = new IsoFree();
    return isoFree.apply(result, params);
  }

  public Set<Graph> apply(Graph g, DifByConstantParameters params) {

    Set<Graph> result = new HashSet<Graph>();
    boolean deriv = false;
    for (Edge edge : g.edges()) {
      if (edge.getColor() == params.color()) {
        Graph newGraph = g.copy();
        newGraph.getEdge(edge.getId()).setColor(params.colorNew);
        result.add(newGraph);
        deriv = true;
      }
    }
    if (deriv) {
    	return result; 
    } else {
    	result.add(g);
    	return result;
    }
  }
  
  public static class DifByConstantParameters implements Parameters {

	  private char color;
	  private char colorNew;

	  public DifByConstantParameters(char color, char colorNew) {

	    this.color = color;
	    this.colorNew = colorNew;
	  }

	  public char color() {

	    return color;
	  }
	  
	  public char colorNew() {

		    return colorNew;
	  }
	}

}