/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.iterators;

import etomica.graph.model.BitmapFactory;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.GraphIterator;

public class FixedEdgeCountIterator implements GraphIterator {

  private byte nodeCount;
  private DefaultPermutator permutator;

  public FixedEdgeCountIterator(byte nodeCount, byte edgeCount) {

    this.nodeCount = nodeCount;
    int[] partition = new int[2];
    // number of absent edges forms partition 0
    partition[0] = (nodeCount * (nodeCount - 1) / 2) - edgeCount;
    // number of present edges forms partition 1
    partition[1] = edgeCount;
    this.permutator = new DefaultPermutator(partition);
  }

  public boolean hasNext() {

    return permutator.hasNext();
  }

  public Graph next() {

    if (hasNext()) {
      return GraphFactory.createGraph(nodeCount, BitmapFactory.createBitmap(permutator.next()));
    }
    return null;
  }

  public void remove() {

    // no-op
  }
}