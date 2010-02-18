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
    byte rangeBegin = 6;
    byte rangeEnd = 6;
    for (byte i = rangeBegin; i <= rangeEnd; i++) {
      testNaive(i);
    }
  }
}