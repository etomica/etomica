package etomica.graph.iterators.filters;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;


import etomica.graph.iterators.ChainedIterator;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;

/**
 * A GlobalFilter is one that must see the entire graph stream before returning any graph.
 * These filters block the outer iteration until the inner iteration is complete.
 */
public abstract class GlobalFilter implements GraphIterator {

  private Map<String, Set<Graph>>  blockingMap = new HashMap<String, Set<Graph>> ();
  private GraphIterator iterator;
  private ChainedIterator blockingIterator = null;

  public GlobalFilter(GraphIterator iterator) {

    this.iterator = iterator;
  }

  protected Set<Graph> getBlockingSet(Graph g) {

    String key = g.getSignature();
    Set<Graph> set = blockingMap.get(key);
    if (set == null) {
      set = new HashSet<Graph>();
      blockingMap.put(key, set);
    }
    return set;
  }

  protected abstract boolean accept(Graph g, Set<Graph> set);

  // hasNext is blocking: it only returns after processing the input iterator completely
  public boolean hasNext() {

    if (blockingIterator == null) {
      while (iterator.hasNext()) {
        Graph next = iterator.next();
        Set<Graph> set = getBlockingSet(next);
        if (accept(next, set)) {
          set.add(next);
        }
      }
      blockingIterator = new ChainedIterator();
      for (String key : blockingMap.keySet()) {
        Set<Graph> set = blockingMap.get(key);
        if (set.size() > 0) {
          blockingIterator.chainIterator(set.iterator());
        }
      }
      blockingIterator.start();
    }
    return blockingIterator.hasNext();
  }

  public Graph next() {

    if (blockingIterator != null && blockingIterator.hasNext()) {
      return blockingIterator.next();
    }
    return null;
  }

  public void remove() {

    // no-op
  }
}