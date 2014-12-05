/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.test;

import java.util.Iterator;

import etomica.virial.cluster2.util.PermutationRangeIterator;

public class PermutatorTest extends CustomTestCase {

  private void println(int[] vector, int base) {

    String result = "";
    for (int i = 0; i < vector.length; i++) {
      result += base - vector[i] - 1;
    }
    System.out.println(result);
  }

  private void println(int[] vector) {

    String result = "";
    for (int i = 0; i < vector.length; i++) {
      result += vector[i];
    }
    System.out.println(result);
  }

  public void testPermutator1() {

// Iterator<int[]> iter = new PermutationIterator(new int[] { 6, 3 });
// while (iter.hasNext()) {
// println(iter.next(), 2);
// }
  }

  public void testPermutator2() {

// Iterator<int[]> iter = new PermutationIterator(new int[] { 4, 3, 3 });
// while (iter.hasNext()) {
// println(iter.next(), 3);
// }
  }

  public void testPermutator3() {

// Iterator<int[]> iter = new PermutationIterator(new int[] { 1, 1, 1, 1, 1, 1
    // });
// while (iter.hasNext()) {
// println(iter.next(), 6);
// }
  }

  public void testPermutator4() {

    int numNodes = 6;
    int numGraphs = (numNodes) * (numNodes - 1) / 2;
    expected = 1 << numGraphs;
    enumerated = 0;
    Iterator<int[]> iter = new PermutationRangeIterator(numGraphs, 0, numGraphs);
    long time1 = System.nanoTime();
    while (iter.hasNext()) {
      iter.next();
      // println(iter.next(), 2);
      enumerated++;
    }
    elapsed = (System.nanoTime() - time1) / K; // µs
    printRuntime();
  }
}