package etomica.virial.cluster2.graph.test;

import etomica.virial.cluster2.graph.EdgesFilter;
import etomica.virial.cluster2.graph.GraphSet;
import etomica.virial.cluster2.graph.GraphFactory;
import etomica.virial.cluster2.graph.Nodes;
import etomica.virial.cluster2.test.CustomTestCase;

public class TestNaiveEdgesGenerator extends CustomTestCase {

  private boolean nullFiltered = false;
  private boolean connectedFiltered = true;
  private GraphSet family;

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
    // TODO: factory for nodes, factory methods for generators and filters in
    // the GraphFamilyFactory
    try {
      // create nodes
      Nodes nodes = GraphFactory.defaultNodes((byte) numNodes);
      // create generator
      EdgesFilter filter = null;
      if (connectedFiltered && nullFiltered) {
        filter = GraphFactory.connectedFilter(nodes, GraphFactory
            .nullFilter(null));
      }
      else if (connectedFiltered) {
        filter = GraphFactory.connectedFilter(nodes, null);
      }
      else if (nullFiltered) {
        filter = GraphFactory.nullFilter(null);
      }
      runGC();
      time1 = System.nanoTime();
      // create family
      GraphFactory.DYNAMIC_ALLOCATION = enumerated < M;
      family = GraphFactory.naiveGraphSet(nodes, filter);
      elapsed = (System.nanoTime() - time1) / K; // µs
      runGC();
      memoryUse();
      enumerated = family.getSize();
      printEnumerated(numNodes);
      printRuntime();
      printMemory();
      if (printPermutations) {
        permutations = family.toString();
        printPermutations();
      }
      System.out.println();
      assertEquals(expected, family.getSize());
    }
    catch (RuntimeException ge) {
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