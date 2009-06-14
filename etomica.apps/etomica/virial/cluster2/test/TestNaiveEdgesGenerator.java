package etomica.virial.cluster2.test;

import etomica.virial.cluster2.graph.*;

public class TestNaiveEdgesGenerator extends CustomTestCase {

  private boolean enumConnected = false;
  private boolean enumBiconnected = false;
  private boolean nullFiltered = false;
  private boolean isomorphFree = false;
  private boolean isomorphOptimal = false;
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
      EdgesFilter filter = ffactory.trueFilter();
      // ignore-all filter
      if (nullFiltered) {
        filter.chain(ffactory.falseFilter());
      }
      // cheap filter for number of edges within a range
      if (minEdges >= 0 && maxEdges >= 0) {
        filter.chain(ffactory.rangeFilter(minEdges, maxEdges));
      }
      // cheap filter for connected graph
      if (enumConnected) {
        filter.chain(ffactory.connectedFilter(nodes));
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
    isomorphOptimal = false;
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

    // OK: 09-06-13
    // templateTest();
  }

  public void testConnected() {

    // OK: 09-06-13
    enumConnected = true;
    // templateTest();
  }

  public void testBiconnected() {

    enumBiconnected = true;
    // templateTest();
  }

  public void testNullFilteredGeneral() {

// OK: 09-06-13
    nullFiltered = true;
// templateTest();
  }

  public void testNullFilteredConnected() {

    nullFiltered = true;
    enumConnected = true;
// templateTest();
  }

  public void testRangedEdgesGeneral() {

    // OK: 09-06-13
    minEdges = 3;
    maxEdges = 4;
// templateTest();
  }

  public void testNullFilteredBiconnected() {

    nullFiltered = true;
    enumBiconnected = true;
// templateTest();
  }

  public void testIsomorphFree() {
// for N=7, about 12 minutes using VF2
    
//    rangeBegin = 7;
//    rangeEnd = 7;
//    isomorphFree = true;
//    isomorphOptimal = false;
//    for (int i = rangeBegin; i <= rangeEnd; i++) {
//      testTemplate(i);
//    }
  }

  public void testOptimizedIsomorphFree() {

    rangeBegin = 8;
    rangeEnd = 8;
    isomorphFree = true;
    isomorphOptimal = false;
    for (int i = rangeBegin; i <= rangeEnd; i++) {
      int numNodes = i * (i - 1) / 2;
      if (isomorphOptimal) {
        for (int j = 0; j <= (numNodes / 2); j++) {
          minEdges = j;
          maxEdges = j;
          testTemplate(i);
        }
      }
      else {
        minEdges = 0;
        maxEdges = numNodes / 2;
        testTemplate(i);
      }
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