/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.iterators;

import java.util.Iterator;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;

public class IteratorWrapper implements GraphIterator {

  private Iterator<Graph> iterator;

  public IteratorWrapper(Iterator<Graph> iterator) {
    this(iterator, false);
  }

  public IteratorWrapper(Iterator<Graph> iterator, boolean doCopy) {
    this.doCopy = doCopy;
    this.iterator = iterator;
  }

  public boolean hasNext() {

    return iterator.hasNext();
  }

  public Graph next() {
    Graph g = iterator.next();
    return doCopy ? g.copy() : g;
  }

  public void remove() {

    // no-op
  }

  protected final boolean doCopy;
}
