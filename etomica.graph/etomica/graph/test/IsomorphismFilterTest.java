package etomica.graph.test;


import etomica.graph.isomorphism.Match;
import etomica.graph.iterators.DefaultIterator;
import etomica.graph.iterators.filters.IsomorphismFilter;
import etomica.graph.iterators.filters.PropertyFilter;
import etomica.graph.model.GraphIterator;
import etomica.graph.property.HasNoRootEdge;

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
    byte rangeBegin = 2;
    byte rangeEnd = 6;
    // nodes = 7: total of 1044 out of 2097152 graphs after 216 secs
    for (byte i = rangeBegin; i <= rangeEnd; i++) {
      expected = i == 0 ? 0 : Match.ISMORPHS_COUNT[i - 1];
      testNaive(i, new IsomorphismFilter(new DefaultIterator(i)));
    }
  }

  public void testIsoFreeGraphsWithRootNodes() {

    reset();
    checkAssertion = false;
    byte rootNodes = 2;
    byte rangeBegin = 3;
    byte rangeEnd = 6;
//    for (byte i = rangeBegin; i <= rangeEnd; i++) {
//      expected = i == 0 ? 0 : Match.ISMORPHS_COUNT[i - 1];
//      testNaive(i, new IsomorphismFilter(new PropertyFilter(new DefaultIterator(i, rootNodes), new HasNoRootEdge())));
//    }
  }
}