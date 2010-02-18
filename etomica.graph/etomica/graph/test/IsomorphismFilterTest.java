package etomica.graph.test;


import etomica.graph.isomorphism.Match;
import etomica.graph.iterators.DefaultIterator;
import etomica.graph.iterators.filters.IsomorphismFilter;
import etomica.graph.model.GraphIterator;

public class IsomorphismFilterTest extends GraphIteratorTest {

  protected void testNaive(byte nodeCount, GraphIterator iterator) {

    testTemplate(nodeCount, iterator);
  }

  public void reset() {

    super.reset();
    printPermutations = true;
    printMemory = true;
    checkAssertion = true;
  }

  public void testIsoFreeGraphs() {

    reset();
    byte rangeBegin = 5;
    byte rangeEnd = 5;
    for (byte i = rangeBegin; i <= rangeEnd; i++) {
      expected = i == 0 ? 0 : Match.ISMORPHS_COUNT[i - 1];
      testNaive(i, new IsomorphismFilter(new DefaultIterator(i)));
    }
  }
}