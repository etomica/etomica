/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
  protected final SignatureMaker signatureMaker;

  public GlobalFilter(GraphIterator iterator) {
    this(iterator, null);
  }
  
  public GlobalFilter(GraphIterator iterator, SignatureMaker signatureMaker) {
    this.iterator = iterator;
    this.signatureMaker = signatureMaker == null ? new SignatureMaker() : signatureMaker;
  }

  protected Set<Graph> getBlockingSet(Graph g) {

    String key = signatureMaker.getSignature(g);
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

  public static class SignatureMaker {
    public String getSignature(Graph g) {
      return g.getSignature();
    }
  }
}