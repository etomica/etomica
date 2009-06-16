package etomica.virial.cluster2.test;

import etomica.virial.cluster2.graph.*;

public class TestNaiveEdgesGenerator extends CustomTestCase {

  private boolean enumConnected = false;
  private boolean enumBiconnected = false;
  private boolean nullFiltered = false;
  private boolean isomorphFree = false;
  private int rangeBegin = 0;
  private int rangeEnd = 0;
  private int minEdges = 0;
  private int maxEdges = 0;
  private int index = 0;
  private int isomorphCount = 0;
  private long factor = 1;
  private GraphSet family;
  private FilterFactory ffactory = new FilterFactory();
  private EdgesSetVisitor printVisitor = new EdgesSetVisitor() {

    public boolean visit(Edges edges) {

      enumerated++;
      isomorphCount += Math
          .round(factor / edges.getMetadata().getCoefficient());
      if (printPermutations) {
        System.out.println(index++ + ": ("
            + Math.round(factor / edges.getMetadata().getCoefficient()) + ") "
            + edges.toString());
      }
      return true;
    }
  };

  private void testTemplate(int numNodes) {

    long time1;
    index = 0;
    factor = 1;
    isomorphCount = 0;
    permutations = "";
    enumerated = 0;
    expected = 1 << (numNodes) * (numNodes - 1) / 2;
    for (int i = 1; i <= numNodes; i++) {
      factor *= i;
    }
    if (nullFiltered) {
      expected = 0;
    }
    try {
      // create nodes
      Nodes nodes = GraphFactory.defaultNodes((byte) numNodes);
      // pass-through filter
// EdgesFilter filter = ffactory.trueFilter();
      EdgesFilter filter = null;
      // ignore-all filter
      if (nullFiltered) {
        filter = ffactory.falseFilter();
      }
      // cheap filter for number of edges within a range
      if (minEdges >= 0 && maxEdges >= 0) {
        EdgesFilter chained = ffactory.rangeFilter(minEdges, maxEdges);
        if (filter != null) {
          filter.chain(chained);
        }
        else {
          filter = chained;
        }
      }
      // cheap filter for connected graph
      if (enumConnected) {
        EdgesFilter chained = ffactory.connectedFilter(nodes);
        if (filter != null) {
          filter.chain(chained);
        }
        else {
          filter = chained;
        }
      }
      runGC();
      time1 = System.nanoTime();
      // create family
      GraphFactory.DYNAMIC_ALLOCATION = enumerated < M;
      family = GraphFactory.naiveGraphSet(nodes, filter, isomorphFree);
      elapsed = (System.nanoTime() - time1) / K; // µs
      runGC();
      memoryUse();
      if (printPermutations) {
        System.out.println();
      }
      family.visitEdgesSet(printVisitor);
      System.out.println(family.getTags());
      printEnumerated(numNodes);
      printRuntime();
      printMemory();
      if (!enumConnected && minEdges < 0 && maxEdges < 0 && !isomorphFree) {
        assertEquals(expected, enumerated, 0.0001);
      }
    }
    catch (RuntimeException ge) {
      ge.printStackTrace();
      fail("Unexpected exception: " + numNodes + ": " + ge.getStackTrace());
    }
  }

  @Override
  public void setUp() {

    super.setUp();
    minEdges = -1;
    maxEdges = -1;
    rangeBegin = 1;
    rangeEnd = 5;
    enumConnected = false;
    enumBiconnected = false;
    nullFiltered = false;
    printMemory = true;
    printPermutations = false;
    isomorphFree = false;
    GraphFactory.USE_UPPER_TRIANGLE = true;
  }

  public void testEmptyFamily() {

    try {
      family = GraphFactory.naiveGraphSet(GraphFactory.emptyNodes(), null);
      fail("Invalid node count");
    }
    catch (Exception ge) {
// ge.printStackTrace();
    }
  }

  protected void baseTest() {

    // OK: 09-06-13
    for (int i = rangeBegin; i <= rangeEnd; i++) {
      testTemplate(i);
    }
  }

  public void testGeneral() {

    printTest("General");
    // OK: 09-06-13
    // baseTest();
  }

  public void testConnected() {

    printTest("Connected");
    // OK: 09-06-13
    enumConnected = true;
    // baseTest();
  }

  public void testBiconnected() {

    printTest("Biconnected");
    enumBiconnected = true;
    // baseTest();
  }

  public void testNullFilteredGeneral() {

    printTest("Null Filtered");
// OK: 09-06-13
    nullFiltered = true;
// baseTest();
  }

  public void testNullFilteredConnected() {

    printTest("Connected, Null Filtered");
    nullFiltered = true;
    enumConnected = true;
// baseTest();
  }

  public void testRangedEdgesGeneral() {

    printTest("Ranged Edges");
    // OK: 09-06-13
    minEdges = 3;
    maxEdges = 4;
// baseTest();
  }

  public void testNullFilteredBiconnected() {

    printTest("Biconnected, Null Filtered");
    nullFiltered = true;
    enumBiconnected = true;
// baseTest();
  }

  public void testIsomorphFree() {

    printTest("Isomorph-Free");
// for N=7, about 12 minutes using VF2
// rangeBegin = 7;
// rangeEnd = 7;
    isomorphFree = true;
    for (int i = rangeBegin; i <= rangeEnd; i++) {
      testTemplate(i);
    }
  }

  public void testIsomorphFreeConnected() {

    printTest("Isomorph-Free, Connected");
    // for N=7, about 12 minutes using VF2
    // rangeBegin = 7;
    // rangeEnd = 7;
    enumConnected = true;
    isomorphFree = true;
    for (int i = rangeBegin; i <= rangeEnd; i++) {
      testTemplate(i);
    }
  }

  public void testOptimizedIsomorphFree() {

    printTest("Isomorph-Free (Optimized)");
// rangeBegin = 8;
// rangeEnd = 8;
    isomorphFree = true;
    for (int i = rangeBegin; i <= rangeEnd; i++) {
      int numNodes = i * (i - 1) / 2;
      minEdges = 0;
      maxEdges = numNodes / 2;
      testTemplate(i);
    }
  }
  // N=8
  // with pre-allocation: out of Java heap space IMMEDIATELY
  // 280M graphs @44 bytes/graph (asymptotically, with BitmapOfLong)
  // requires approximately 12GB of heap space
  // with dynamic allocation, 14sec to process all graphs
  //
  // N=9 requires approximately 6.5h
}