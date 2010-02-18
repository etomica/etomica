package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.iterators.IteratorWrapper;
import etomica.graph.iterators.filters.IsomorphismFilter;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;

/*
 * IsoFree (Isomorphism Free GraphSet): GraphSet --> GraphSet
 *
 * Given a graph set S, IsoFree(S) returns a graph set S' such that no two graphs g1, g2
 * in S' are isomorphic to each other. For every graph g' in S', the coefficient of g'
 * is the sum of the coefficients of every graph g in S such that g is isomorphic to g'.
 */
public class IsoFree implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {

    PCopy pcopy = new PCopy();
    IteratorWrapper wrapper = new IteratorWrapper(pcopy.apply(argument, null).iterator());
    GraphIterator isomorphs = new IsomorphismFilter(wrapper);
    Set<Graph> result = new HashSet<Graph>();
    while (isomorphs.hasNext()) {
      result.add(isomorphs.next());
    }
    return result;
  }
}