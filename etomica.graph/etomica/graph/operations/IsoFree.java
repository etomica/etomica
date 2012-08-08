package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.iterators.IteratorWrapper;
import etomica.graph.iterators.filters.GlobalFilter;
import etomica.graph.iterators.filters.IsomorphismFilter;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;

public class IsoFree implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {

    IteratorWrapper wrapper = new IteratorWrapper(argument.iterator(), true);
    GraphIterator isomorphs = new IsomorphismFilter(wrapper, (IsoFreeParams)params);
    Set<Graph> result = new HashSet<Graph>();
    while (isomorphs.hasNext()) {
      result.add(isomorphs.next());
    }
    return result;
  }

  public static class IsoFreeParams extends GlobalFilter.SignatureMaker implements Parameters {
  }
}