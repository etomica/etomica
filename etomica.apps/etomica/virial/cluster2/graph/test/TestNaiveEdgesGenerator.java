package etomica.virial.cluster2.graph.test;

import etomica.virial.cluster2.graph.EdgesFilter;
import etomica.virial.cluster2.graph.EdgesGenerator;
import etomica.virial.cluster2.graph.GraphException;
import etomica.virial.cluster2.graph.impl.ConnectedEdgesFilter;
import etomica.virial.cluster2.graph.impl.NaiveEdgesGenerator;
import etomica.virial.cluster2.graph.impl.SimpleGraphSet;
import etomica.virial.cluster2.graph.impl.NullEdgesFilter;
import etomica.virial.cluster2.test.CustomTestCase;

public class TestNaiveEdgesGenerator extends CustomTestCase {

  private boolean nullFiltered = false;
  private boolean connectedFiltered = true;
  private SimpleGraphSet family;
  private EdgesGenerator generator;

  private void testCompleteFamilyN(int numNodes) {

    long time1;
    printMemory = true;
    printPermutations = false;
    permutations = "";
    expected = 1 << (numNodes) * (numNodes - 1) / 2;
    enumerated = 1 << (numNodes) * (numNodes - 1) / 2;
    if (nullFiltered) {
      expected = 0;
    }
    try {
      family = new SimpleGraphSet((byte) numNodes, enumerated < M);
      EdgesFilter filter = null;
      if (connectedFiltered && nullFiltered) {
        filter = new ConnectedEdgesFilter(family.getNodes(), new NullEdgesFilter());
      } else if (connectedFiltered) {
        filter = new ConnectedEdgesFilter(family.getNodes());
      } else if (nullFiltered) {
        filter = new NullEdgesFilter();
      } else {
        generator = new NaiveEdgesGenerator(family);
      }
      if (filter == null) {
        generator = new NaiveEdgesGenerator(family);
      } else {
        generator = new NaiveEdgesGenerator(family, filter);
      }
      runGC();
      time1 = System.nanoTime();
      family.generateEdgesSet(generator);
      elapsed = (System.nanoTime() - time1) / K; // µs
      runGC();
      memoryUse();
      enumerated = family.getEdgesSetSize();
      printEnumerated(numNodes);
      printRuntime();
      printMemory();
      if (printPermutations) {
        permutations = family.toString();
        printPermutations();
      }
      System.out.println();
      assertEquals(expected, family.getEdgesSetSize());
    } catch (GraphException ge) {
      fail("Unexpected exception: " + ge.getStackTrace());
    }
  }

  public void testCompleteFamily0() {

    // try {
    // family = new SimpleGraphSet((byte) 0);
    // fail("Invalid node count");
    // } catch (GraphException ge) {
    // }
  }

  public void testCompleteFamily1() {

    testCompleteFamilyN(1);
  }

  public void testCompleteFamily2() {

    testCompleteFamilyN(2);
  }

  public void testCompleteFamily3() {

    testCompleteFamilyN(3);
  }

  public void testCompleteFamily4() {

    testCompleteFamilyN(4);
  }

  public void testCompleteFamily5() {

    testCompleteFamilyN(5);
  }

  public void testCompleteFamily6() {

    testCompleteFamilyN(6);
  }

  public void testCompleteFamily7() {

    testCompleteFamilyN(7);
  }

  public void testCompleteFamily8() {

    // with pre-allocation: out of Java heap space IMMEDIATELY
    // 280M graphs @44 bytes/graph (asymptotically, with BitmapOfLong)
    // requires approximately 12GB of heap space
    // with dynamic allocation, 14sec to process all graphs
    // if (nullFiltered ) {
    // testCompleteFamilyN(8 );
    // }
  }

  public void testCompleteNullFamily9() {
    // requires approximately 6.5h

    // if (nullFiltered ) {
    // testCompleteFamilyN(9 );
    // }
  }
}