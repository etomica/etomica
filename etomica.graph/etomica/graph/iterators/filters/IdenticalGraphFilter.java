package etomica.graph.iterators.filters;

import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;

/**
 * Filter that combines graphs that are identical
 *
 * @author Andrew Schultz
 */
public class IdenticalGraphFilter extends GlobalFilter {

  public IdenticalGraphFilter(GraphIterator iterator) {

    super(iterator);
  }

  protected boolean accept(Graph g1, Set<Graph> set) {

    boolean result = true;
    if (!set.isEmpty()) {
      for (Graph isoGraph : set) {
        // test for identical graphs and, if they do match, update the graph
        // coefficient and reject the new graph
        if (isoGraph.compareTo(g1) == 0) {
          isoGraph.coefficient().add(g1.coefficient());
          result = false;
          if (isoGraph.coefficient().getNumerator() == 0) {
            set.remove(isoGraph);
            break;
          }
        }
      }
    }
    return result;
  }
}