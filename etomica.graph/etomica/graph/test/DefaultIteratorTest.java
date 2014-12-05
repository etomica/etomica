/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.test;


import etomica.graph.iterators.DefaultIterator;
import etomica.graph.model.GraphIterator;

public class DefaultIteratorTest extends GraphIteratorTest {

  protected void testNaive(byte nodeCount, GraphIterator iterator) {

    testTemplate(nodeCount, iterator);
  }

  public void reset() {

    super.reset();
    printPermutations = true;
    printMemory = true;
    checkAssertion = true;
  }

  public void testGraphs() {

    reset();
    // nodes = 7: total of 2097152 graphs after 2 secs (0 min)
    byte rangeBegin = 5;
    byte rangeEnd = 5;
    for (byte i = rangeBegin; i <= rangeEnd; i++) {
      expected = i == 0 ? 0 : 1 << (i) * (i - 1) / 2;
      testNaive(i, new DefaultIterator(i));
    }
  }

  public void testGraphsWithRootNodes() {

    reset();
    checkAssertion = false;
    byte rootNodes = 2;
    byte rangeBegin = 4;
    byte rangeEnd = 6;
    for (byte i = rangeBegin; i <= rangeEnd; i++) {
      expected = i == 0 ? 0 : 1 << (i) * (i - 1) / 2;
//      testNaive(i, new DefaultIterator(i, rootNodes));
    }
  }
}