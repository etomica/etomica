package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.iterators.IteratorWrapper;
import etomica.graph.iterators.filters.FieldNodeCount;
import etomica.graph.model.Graph;
import etomica.graph.model.Node;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;

/**
 * Performs decoration of diagrams by points.
 *
 * @author Andrew Schultz
 */
public class Decorate implements Binary {

  public Set<Graph> apply(Set<Graph> argument, Set<Graph> argument2, Parameters params) {

    assert (params instanceof DecorateParameters);
    char color = ((DecorateParameters)params).color;
    MulFlexibleParameters mfp = ((DecorateParameters)params).mfp;
    Set<Graph> result = new HashSet<Graph>();
    int maxNodes = 0;
    int maxPow = 0;
    for (Graph g : argument) {
      if (g.nodeCount() > maxNodes) {
        maxNodes = g.nodeCount();
      }
      int colorCount = 0;
      for (Node node : g.nodes()) {
        if (node.getColor() == color) {
          colorCount++;
        }
      }
      if (colorCount > maxPow) {
        maxPow = colorCount;
      }
    }
    Set<Graph>[] allSet1 = new Set[maxPow+1];
    for (int i=1; i<allSet1.length; i++) {
      allSet1[i] = new HashSet<Graph>();
    }
    for (Graph g : argument) {
      int colorCount = 0;
      for (Node node : g.nodes()) {
        if (node.getColor() == color) {
          colorCount++;
        }
      }
      if (colorCount == 0) {
        result.add(g);
      }
      else {
        allSet1[colorCount].add(g);
      }
    }
    // set2Pow holds (set2)^i
    Set<Graph> set2Pow = argument2;
    MulFlexible mulFlex = new MulFlexible();
    IsoFree isoFree = new IsoFree();
    for (int i=1; i<allSet1.length; i++) {
      result.addAll(mulFlex.apply(allSet1[i], set2Pow, mfp));

      FieldNodeCount truncater = new FieldNodeCount(new IteratorWrapper(result.iterator()), maxNodes);
      Set<Graph> truncatedSet = new HashSet<Graph>();
      while (truncater.hasNext()) {
        truncatedSet.add(truncater.next());
      }
      result = isoFree.apply(truncatedSet, null);
      
      if (i+1<allSet1.length) {
        // we're going to make another pass.  calculate (set2)^(i+1)
        set2Pow = mulFlex.apply(set2Pow, argument2, mfp);

        // truncate to avoid diagram explosion
        truncater = new FieldNodeCount(new IteratorWrapper(set2Pow.iterator()), maxNodes-i);
        truncatedSet = new HashSet<Graph>();
        while (truncater.hasNext()) {
          truncatedSet.add(truncater.next());
        }
        set2Pow = isoFree.apply(truncatedSet, null);
      }
    }
    return result;
  }

  public static class DecorateParameters implements Parameters {
    public final char color;
    public final MulFlexibleParameters mfp;
    public DecorateParameters(char color, MulFlexibleParameters mfp) {
      this.color = color;
      this.mfp = mfp;
    }
  }
}
