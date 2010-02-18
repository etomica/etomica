package etomica.graph.iterators.filters;

import java.util.Set;


import etomica.graph.isomorphism.Match;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;

public class IsomorphismFilter extends GlobalFilter {

  public static boolean DEBUG_MODE = true;
  private static int DEBUG_FREQUENCY = 1000;

  private int countSeen = 0;
  private int countUnique = 0;
  private int countDiscarded = 0;
  private long debugStart = System.nanoTime();

  public IsomorphismFilter(GraphIterator iterator) {

    super(iterator);
  }

  protected boolean accept(Graph g1, Set<Graph> set) {

    boolean result = true;
    if (!set.isEmpty()) {
      for (Graph isoGraph : set) {
        // test for isomorphism and, if they don't match, keep the graph lower in the
        // graph order; update the graph coefficients;
        if (Match.match(isoGraph, g1)) {
          result = false;
          if (isoGraph.compareTo(g1) <= 0) {
            isoGraph.coefficient().add(g1.coefficient());
          }
          else {
            set.remove(isoGraph);
            g1.coefficient().add(isoGraph.coefficient());
            // replace the graph in the set with an isomorph with lower score
            result = true;
          }
          break;
        }
      }
    }
    if (!result) {
      countDiscarded++;
    }
    else {
      countUnique++;
    }
    countSeen++;
    debugReport();
    return result;
  }

  private void debugReport() {

    if (!DEBUG_MODE || (countSeen % DEBUG_FREQUENCY != 0)) {
      return;
    }
    long debugDuration = (System.nanoTime() - debugStart) / 1000000000;
    System.out.println(String.format("unique: %d; discarded: %d; total: %d; time: %d sec (%d min)",
        countUnique, countDiscarded, countSeen, debugDuration, (debugDuration / 60)));
  }
}