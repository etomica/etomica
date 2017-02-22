/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.iterators;

import etomica.graph.model.Bitmap;
import etomica.graph.model.BitmapFactory;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.GraphIterator;

public class DefaultIterator implements GraphIterator {

  private Bitmap current;
  private Bitmap target;
  private byte nodeCount;
  private byte rootNodeCount;

  public DefaultIterator(byte nodeCount) {

    this(nodeCount, (byte) 0);
  }

  public DefaultIterator(byte nodeCount, byte rootNodeCount) {

    this.nodeCount = nodeCount;
    this.rootNodeCount = rootNodeCount;
    if (nodeCount == 0) {
      current = BitmapFactory.EMPTY;
      target = BitmapFactory.EMPTY;
    }
    else if (nodeCount == 1) {
      current = BitmapFactory.ZERO;
      target = BitmapFactory.ZERO;
    }
    else {
      current = BitmapFactory.createBitmap(nodeCount, false);
      target = BitmapFactory.createBitmap(nodeCount, true);
    }
  }

  public boolean hasNext() {

    return (!current.equals(BitmapFactory.EMPTY));
  }

  public Graph next() {

    if (hasNext()) {
      Graph g = GraphFactory.createGraph(nodeCount, rootNodeCount, current);
      if (current.equals(target)) {
        current = BitmapFactory.EMPTY;
      }
      else {
        current = current.copy();
        current.inc();
      }
      return g;
    }
    return null;
  }

  public void remove() {

    // TODO: no-op
  }
}