package etomica.graph.operations;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import etomica.graph.iterators.RangePermutator;
import etomica.graph.model.Edge;
import etomica.graph.model.Graph;
import etomica.graph.model.Permutator;

public class Split implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {

    assert (params instanceof SplitParameters);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      Set<Graph> newSet = apply(g, (SplitParameters) params);
      if (newSet != null) {
        result.addAll(newSet);
      }
    }
    Unary isoFree = new IsoFree();
    return isoFree.apply(result, params);
  }

  public Set<Graph> apply(Graph graph, SplitParameters params) {

    // collect the Ids of all edges we must replace
    List<Byte> edges = new ArrayList<Byte>();
    for (Edge edge : graph.edges()) {
      if (edge.getColor() == params.edgeColor()) {
        edges.add(edge.getId());
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
      for (int edgePtr = 0; edgePtr < edges.size(); edgePtr++) {
        char newColor = permutation[edgePtr] == 0 ? params.newColor0() : params.newColor1();
        newGraph.getEdge(edges.get(edgePtr)).setColor(newColor);
      }
      result.add(newGraph);
    }
    return result;
  }
}