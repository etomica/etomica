package etomica.graph.iterators;

import java.util.Iterator;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;

public class IteratorWrapper implements GraphIterator {

  private Iterator<Graph> iterator;

  public IteratorWrapper(Iterator<Graph> iterator) {

    this.iterator = iterator;
  }

  public boolean hasNext() {

    return iterator.hasNext();
  }

  public Graph next() {

    return iterator.next();
  }

  public void remove() {

    // no-op
  }
}
