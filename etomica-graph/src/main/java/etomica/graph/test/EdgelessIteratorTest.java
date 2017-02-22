/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.test;

import java.util.HashMap;
import java.util.Map;

import static etomica.graph.model.Metadata.*;


import etomica.graph.iterators.EdgelessIterator;
import etomica.graph.iterators.filters.IsomorphismFilter;
import etomica.graph.model.GraphIterator;

public class EdgelessIteratorTest extends GraphIteratorTest {

  private Map<Character, Byte> fieldMap = new HashMap<Character, Byte>();
  private Map<Character, Byte> rootMap = new HashMap<Character, Byte>();

  protected void testNaive(byte nodeCount, GraphIterator iterator) {

    testTemplate(nodeCount, iterator);
  }

  public void reset() {

    super.reset();
    printPermutations = true;
    printMemory = true;
    checkAssertion = false;
    rootMap.clear();
    rootMap.put(COLOR_CODE_0, (byte) 2);
    fieldMap.clear();
    fieldMap.put(COLOR_CODE_0, (byte) 2);
    fieldMap.put(COLOR_CODE_1, (byte) 1);
  }

  public void testGraphs() {

    reset();
    // nodes = 5: total of 6 graphs after < 50ms
    testNaive((byte) 5, new EdgelessIterator(rootMap, fieldMap));
  }

  public void testIsomorphFreeGraphs() {

    reset();
    // nodes = 5: total of 1 graphs after < 50ms
    testNaive((byte) 5, new IsomorphismFilter(new EdgelessIterator(rootMap, fieldMap)));
  }
}