/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.util;

import java.util.Iterator;
import java.util.NoSuchElementException;

public class PermutationRangeIterator implements Iterator<int[]> {

  private int cur;
  private PermutationIterator[] permutations;

  public PermutationRangeIterator(int setSize, int minP1, int maxP1) {

    assert ((minP1 >= 0) && (maxP1 >= minP1) && (setSize >= maxP1));
    permutations = new PermutationIterator[maxP1 - minP1 + 1];
    for (int i = minP1; i <= maxP1; i++) {
      permutations[i - minP1] = new PermutationIterator(new int[] { i,
          setSize - i });
    }
    cur = 0;
  }

  public boolean hasNext() {

    return (cur < permutations.length) && (permutations[cur].hasNext());
  }

  // pop the current permutation then push the next one
  public int[] next() {

    if (!hasNext()) {
      throw new NoSuchElementException();
    }
    int[] result = permutations[cur].next();
    if (!permutations[cur].hasNext()) {
      cur++;
    }
    return result;
  }

  // this iterator is read-only
  public void remove() {

    throw new UnsupportedOperationException();
  }
}