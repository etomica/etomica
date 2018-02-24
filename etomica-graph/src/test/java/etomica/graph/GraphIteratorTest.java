/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph;

import java.util.Iterator;

import etomica.graph.model.Graph;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.fail;

public class GraphIteratorTest extends CustomTestCase {

  protected void testTemplate(byte nodeCount, Iterator<Graph> iterator) {

    long time1;
    permutations = "";
    enumerated = 0;
    try {
      runGC();
      time1 = System.nanoTime();
      // iterate
      while (iterator.hasNext()) {
        enumerated++;
        if (printPermutations) {
          System.out.println(iterator.next());
        }
        else {
          iterator.next();
        }
      }
      elapsed = (System.nanoTime() - time1) / K; // µs
      memoryUse();
      runGC();
      printEnumerated(nodeCount);
      printRuntime();
      printMemory();
      if (checkAssertion) {
        assertEquals(expected, enumerated, 0.0001);
      }
    }
    catch (RuntimeException ge) {
      ge.printStackTrace();
      fail("Unexpected exception: " + nodeCount + ": " + ge.getStackTrace());
    }
  }

  public void reset() {

    printMemory = false;
    printRuntime = true;
    printPermutations = true;
  }
}
