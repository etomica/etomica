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

}