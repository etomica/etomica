/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.iterators;

import java.util.NoSuchElementException;

import etomica.graph.model.Permutator;

/**
 * This is a facility class that iterates through (maxP1 - minP1 + 1) permutations one at
 * a time. When the last element of the K-th permutation is produced, the next element is
 * the first element of the (K+1)-th permutation.
 *
 * Each permutation iteration consists of two classes with a fixed total number of elements.
 * For example, an instance of this class created as new PermutationRangeIterator(8, 0, 8)
 * would generate all combinations of 8 bits.
 *
 */
public class RangePermutator implements Permutator {

  private int cur;
  private DefaultPermutator[] permutations;

  public RangePermutator(int setSize, int minP1, int maxP1) {

    assert ((minP1 >= 0) && (maxP1 >= minP1) && (setSize >= maxP1));
    permutations = new DefaultPermutator[maxP1 - minP1 + 1];
    for (int i = minP1; i <= maxP1; i++) {
      permutations[i - minP1] = new DefaultPermutator(new int[] { i, setSize - i });
    }
    cur = 0;
  }

  public boolean hasNext() {

    return (cur < permutations.length) && (permutations[cur].hasNext());
  }

  // pop the current permutation and push the next one
  public byte[] next() {

    if (!hasNext()) {
      throw new NoSuchElementException();
    }
    byte[] result = permutations[cur].next();
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