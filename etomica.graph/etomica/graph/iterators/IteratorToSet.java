package etomica.graph.iterators;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;

public class IteratorToSet {

  public Set<Graph> getSet(GraphIterator iterator) {

    Set<Graph> result = new HashSet<Graph>();
    while (iterator.hasNext()) {
      result.add(iterator.next());
    }
    return result;
  }
}