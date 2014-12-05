/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.test;


import etomica.graph.iterators.DefaultIterator;
import etomica.graph.iterators.filters.PropertyFilter;
import etomica.graph.model.GraphIterator;
import etomica.graph.property.HasEdgeCount;

public class EdgeCountFilterTest extends GraphIteratorTest {

  protected void testNaive(byte nodeCount, GraphIterator iterator) {

    testTemplate(nodeCount, iterator);
  }

  public void reset() {

    super.reset();
    printPermutations = false;
    printMemory = true;
    checkAssertion = true;
  }

  public void testEdgeFilteredGraphs() {

    reset();
    // nodes = 7: total of 2097152 graphs after 56 secs (0 min)
    byte nodes = 7;
    checkAssertion = false;
    int count = 0;
    long time = System.nanoTime();
    for (byte edges = 0; edges <= nodes * (nodes - 1) / 2; edges++) {
      System.out.println("edges....: " + edges);
      System.out.println("partial..: " + count);
      testNaive(nodes, new PropertyFilter(new DefaultIterator(nodes), new HasEdgeCount(edges)));
      count += enumerated;
    }
    long elapsed = (System.nanoTime() - time) / 1000000000;
    System.out.println("total of " + count + " graphs after " + elapsed + " secs (" + elapsed / 60 + " min)");
  }
}