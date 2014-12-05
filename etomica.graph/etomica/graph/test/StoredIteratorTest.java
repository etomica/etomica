/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.test;

import etomica.graph.iterators.StoredIterator;

public class StoredIteratorTest extends GraphIteratorTest {

  protected void testNaive(byte nodeCount) {

    testTemplate(nodeCount, new StoredIterator(nodeCount));
  }

  public void reset() {

    super.reset();
    checkAssertion = false;
    printPermutations = true;
    printMemory = true;
  }

  public void testGraphs() {

    reset();
    byte rangeBegin = 5;
    byte rangeEnd = 5;
    for (byte i = rangeBegin; i <= rangeEnd; i++) {
      testNaive(i);
    }
  }
}