/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import etomica.graph.iterators.RangePermutator;
import etomica.graph.model.Edge;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphList;
import etomica.graph.model.Permutator;
import etomica.graph.property.Property;

public class Split implements Unary {

  protected final IsoFree isoFree = new IsoFree();

  
  public Set<Graph> apply(Set<Graph> argument, Parameters params) {

    assert (params instanceof SplitParameters);
    Set<Graph> result = new GraphList(null);
    int count = 0;
    long t1 = System.currentTimeMillis();
    int interval = 10000;
    int totalCount = 0;
    // We assume here that the set we're receiving is isofree.  We can then
    // assume that graphs that result from splitting can be isomers only if
    // they come from the same original graph.
    for (Graph g : argument) {
      Set<Graph> newSet = apply(g, (SplitParameters) params);
              
      if (newSet.size() == 1) {
        // split had no effect (g did not contain the bond of interest)
        result.addAll(newSet);
      }
      else {
        result.addAll(isoFree.apply(newSet, null));
      }
      if (++count==interval) {
        if (totalCount == 0) System.out.print("Split =>");
        System.out.print(" "+result.size());
        totalCount += count;
        if (totalCount == interval*10) interval *= 10;
        count = 0;
      }
    }
    long t2 = System.currentTimeMillis();
    if (totalCount>0) System.out.println("; "+(t2-t1)/1000+"s");
    return result;
  }

  public Set<Graph> apply(Graph graph, SplitParameters params) {
    Property dp = params.getDiscardProperty();

    // collect the Ids of all edges we must replace
    List<Byte> edges = new ArrayList<Byte>();
    for (Edge edge : graph.edges()) {
      if (edge.getColor() == params.edgeColor()) {
        edges.add(edge.getId());
      }
    }

    Set<Graph> result = new GraphList(null);
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
        char newColor = permutation[edgePtr] == 0 ? params.newColor0() : params.newColor1();
        newGraph.getEdge(edges.get(edgePtr)).setColor(newColor);
      }
      if (dp == null || !dp.check(newGraph)) {
        result.add(newGraph);
      }
    }
    return result;
  }
}