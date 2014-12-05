/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import etomica.graph.iterators.RangePermutator;
import etomica.graph.model.Graph;
import etomica.graph.model.Permutator;

/**
 * Performs the substitution 1=a+b.  Can be used to perform
 * 1 = e + (-f)
 * (subsequently, you'll need to replace (-f) with (f) and
 * multiply the diagram by -1 if it had an odd number of f bonds).
 *
 * @author Andrew Schultz
 */
public class SplitOne implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {

    assert (params instanceof SplitParameters);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      Set<Graph> newSet = apply(g, (SplitOneParameters) params);
      if (newSet != null) {
        result.addAll(newSet);
      }
    }
    Unary isoFree = new IsoFree();
    result = isoFree.apply(result, null);
    return result;
  }

  public Set<Graph> apply(Graph graph, SplitOneParameters params) {

    // collect the Ids of all edges we must replace
    List<Byte> edges = new ArrayList<Byte>();
    byte nodeCount = graph.nodeCount();
    for (byte edgeId = 0; edgeId < nodeCount*(nodeCount-1)/2; edgeId++) {
      if (!graph.hasEdge(edgeId)) {
        edges.add(edgeId);
      }
    }

    Set<Graph> result = new HashSet<Graph>();
    // compute all possible permutations of length |edges| consisting of two colors
    // such that each of the colors can appear any number of times from 0 to |edges|
    Permutator permutations = new RangePermutator(edges.size(), 0, edges.size());
    // each permutation is a color assignment for the edges
    while (permutations.hasNext()) {
      // copy the graph
      Graph newGraph = graph.copy();
      byte[] permutation = permutations.next();
      // modify edge colors: partition 0 => newColor0, partition 1 => newColor1
      for (byte edgePtr = 0; edgePtr < edges.size(); edgePtr++) {
        char newColor = permutation[edgePtr] == 0 ? params.newColor0 : params.newColor1;
        byte edgeId = edges.get(edgePtr);
        newGraph.putEdge(edgeId);
        newGraph.getEdge(edgeId).setColor(newColor);
      }
      result.add(newGraph);
    }
    return result;
  }
  
  public static class SplitOneParameters implements Parameters {
    
    public final char newColor0;
    public final char newColor1;
  
    public SplitOneParameters(char newColor0, char newColor1) {
      this.newColor0 = newColor0;
      this.newColor1 = newColor1;
    }
  }
}