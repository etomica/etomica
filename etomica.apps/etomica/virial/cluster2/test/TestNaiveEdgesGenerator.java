package etomica.virial.cluster2.test;

import etomica.virial.cluster2.graph.*;

public class TestNaiveEdgesGenerator extends CustomTestCase {

  private boolean enumConnected = false;
// private boolean enumBiconnected = false;
  private boolean nullFiltered = false;
  private int rangeEnd = 0;
  private int minEdges = 0;
  private int maxEdges = 0;
  private GraphSet family;

  private void testTemplate(int numNodes) {

    long time1;
    permutations = "";
    expected = 1 << (numNodes) * (numNodes - 1) / 2;
    enumerated = 1 << (numNodes) * (numNodes - 1) / 2;
    if (nullFiltered) {
      expected = 0;
    }
    try {
      // create nodes
      Nodes nodes = GraphFactory.defaultNodes((byte) numNodes);
      // pass-through filter
      EdgesFilter filter = GraphFactory.trueFilter();
      if (minEdges != maxEdges) {
        filter.chain(GraphFactory.rangeFilter(minEdges, maxEdges));
      }
      if (enumConnected) {
        filter.chain(GraphFactory.connectedFilter(nodes));
      }
      if (nullFiltered) {
        filter.chain(GraphFactory.falseFilter());
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
      System.out.println(family.getTags());
      printEnumerated(numNodes);
      printRuntime();
      printMemory();
      if (printPermutations) {
        permutations = family.toString();
        printPermutations();
      }
      if (!enumConnected && minEdges == maxEdges) {
        assertEquals(expected, family.getSize());
      }
    }
    catch (RuntimeException ge) {
      fail("Unexpected exception: " + ge.getStackTrace());
    }
  }

  @Override
  public void setUp() {

    super.setUp();
    enumConnected = false;
// enumBiconnected = false;
    nullFiltered = false;
    printMemory = true;
    printPermutations = true;
    GraphFactory.USE_UPPER_TRIANGLE = true;
    rangeEnd = 5;
  }

  public void testEmptyFamily() {

    try {
      family = GraphFactory.naiveGraphSet(GraphFactory.defaultNodes((byte) 0),
          null);
      fail("Invalid node count");
    }
    catch (Exception ge) {
// ge.printStackTrace();
    }
  }

  public void testGeneral() {

    for (int i = 1; i <= rangeEnd; i++) {
      testTemplate(i);
    }
  }

  public void testConnected() {

    enumConnected = true;
    testGeneral();
  }

// public void testBiconnected() {
//
// enumBiconnected = true;
// testGeneral();
// }
//
  public void testNullFilteredGeneral() {

    nullFiltered = true;
    testGeneral();
  }

  public void testNullFilteredConnected() {

    nullFiltered = true;
    enumConnected = true;
    testGeneral();
  }

  public void testRangedEdgesGeneral() {

    minEdges = 3;
    maxEdges = 4;
    testGeneral();
  }
// public void testNullFilteredBiconnected() {
//
// nullFiltered = true;
// enumBiconnected = true;
// testGeneral();
// }
//
  // N=8
  // with pre-allocation: out of Java heap space IMMEDIATELY
  // 280M graphs @44 bytes/graph (asymptotically, with BitmapOfLong)
  // requires approximately 12GB of heap space
  // with dynamic allocation, 14sec to process all graphs
  //
  // N=9 requires approximately 6.5h
}