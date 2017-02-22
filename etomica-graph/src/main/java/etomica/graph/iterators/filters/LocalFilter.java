/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.iterators.filters;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;

public abstract class LocalFilter implements GraphIterator {

  private Graph next;
  private GraphIterator iterator;

  public LocalFilter(GraphIterator iterator) {

    this.next = null;
    this.iterator = iterator;
  }

  protected abstract boolean accept(Graph g);

  public boolean hasNext() {

    while (next == null && iterator.hasNext()) {
      next = iterator.next();
      if (!accept(next)) {
        next = null;
      }
    }
    return (next != null);
  }

  public Graph next() {

    Graph result = next;
    next = null;
    return result;
  }

  public void remove() {

    // no-op
  }
}