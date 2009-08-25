package etomica.virial.cluster2.test;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.EdgesFilter;
import etomica.virial.cluster2.graph.EdgesSetVisitor;
import etomica.virial.cluster2.graph.FilterFactory;
import etomica.virial.cluster2.graph.GraphFactory;
import etomica.virial.cluster2.graph.GraphSet;
import etomica.virial.cluster2.graph.Nodes;
import etomica.virial.cluster2.graph.isomorphism.Match;

public class TestNaiveEdgesGenerator extends CustomTestCase {

  private boolean enumConnected = false;
  private boolean enumBiconnected = false;
  private boolean dropRootEdges = false;
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
          .round(factor / edges.getMetadata().getCoefficient().getValue1());
      if (printPermutations) {
        System.out.println(index++ + ": ("
            + Math.round(factor / edges.getMetadata().getCoefficient().getValue1()) + ") "
            + edges.toString());
      }
      return true;
    }
  };

  private void testTemplate(Nodes nodes) {

    long time1;
    byte numNodes = nodes.count();
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
      // pass-through filter
// EdgesFilter filter = ffactory.trueFilter();
      EdgesFilter filter = null;
      // ignore-all filter
      if (nullFiltered) {
        filter = ffactory.falseFilter();
      }
      // cheap filter for dropping root edges
      if (dropRootEdges) {
        EdgesFilter chained = ffactory.rootEdgesFilter(nodes);
        if (filter != null) {
          filter.chain(chained);
        }
        else {
          filter = chained;
        }
      }
      // cheap filter for number of edges within a range
      if (minEdges >= 0 && maxEdges >= 0) {
        EdgesFilter chained = null;
        if (isomorphFree) {
//          chained = ffactory.rangedFieldEdgesFilter(nodes, 0, 0);
//          chained.chain(ffactory.rangedRootEdgesFilter(nodes, minEdges,
//              maxEdges));
// chained = ffactory.rangedFieldEdgesFilter(nodes, minEdges, maxEdges);
 chained = ffactory.rangeFilter(minEdges, maxEdges);
// chained.chain(ffactory.rangedRootEdgesFilter(nodes, 0, 0));
        }
        else {
          chained = ffactory.rangeFilter(minEdges, maxEdges);
        }
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
      if (!enumConnected && minEdges < 0 && maxEdges < 0 && !isomorphFree
          && !dropRootEdges) {
        assertEquals(expected, enumerated, 0.0001);
      }
    }
    catch (RuntimeException ge) {
      ge.printStackTrace();
      fail("Unexpected exception: " + numNodes + ": " + ge.getStackTrace());
    }
  }

  @Override
  public void reset() {

    minEdges = -1;
    maxEdges = -1;
    rangeBegin = 2;
    rangeEnd = 5;
    enumConnected = false;
    enumBiconnected = false;
    nullFiltered = false;
    printMemory = true;
    printPermutations = true;
    isomorphFree = false;
    GraphFactory.USE_UPPER_TRIANGLE = false;
  }

  public void testEmptyFamily() {

    reset();
    try {
      family = GraphFactory.naiveGraphSet(GraphFactory.emptyNodes(), null);
      fail("Invalid node count");
    }
    catch (Exception ge) {
// ge.printStackTrace();
    }
  }

  protected Nodes getNodes(byte numNodes) {

    return GraphFactory.defaultNodes(numNodes);
  }

  protected Nodes getNodes(byte fieldNodes, byte rootNodes) {

    return GraphFactory.defaultNodes(fieldNodes, rootNodes);
  }

  protected void baseTest() {

    // OK: 09-06-13
    for (int i = rangeBegin; i <= rangeEnd; i++) {
      testTemplate(getNodes((byte) i));
    }
  }

  public void testGeneral() {

    reset();
    printTest("General");
    // OK: 09-06-13
//    baseTest();
  }

  public void testConnected() {

    reset();
    printTest("Connected");
    // OK: 09-06-13
    enumConnected = true;
    // baseTest();
  }

  public void testBiconnected() {

    reset();
    printTest("Biconnected");
    enumBiconnected = true;
    // baseTest();
  }

  public void testNullFilteredGeneral() {

    reset();
    printTest("Null Filtered");
// OK: 09-06-13
    nullFiltered = true;
// baseTest();
  }

  public void testNullFilteredConnected() {

    reset();
    printTest("Connected, Null Filtered");
    nullFiltered = true;
    enumConnected = true;
// baseTest();
  }

  public void testRangedEdgesGeneral() {

    reset();
    printTest("Ranged Edges");
    // OK: 09-06-13
    minEdges = 3;
    maxEdges = 4;
// baseTest();
  }

  public void testNullFilteredBiconnected() {

    reset();
    printTest("Biconnected, Null Filtered");
    nullFiltered = true;
    enumBiconnected = true;
// baseTest();
  }

  public void testIsomorphFree() {

    reset();
    printTest("Isomorph-Free");
// for N=7, about 12 minutes using VF2
// rangeBegin = 4;
// rangeEnd = 4;
    isomorphFree = true;
    for (int i = rangeBegin; i <= rangeEnd; i++) {
 testTemplate(getNodes((byte) i));
    }
  }

  /**
   * General graphs with 2 root nodes
   */
  public void testGeneralR2() {

    reset();
    printTest("2 root nodes");
// for N=7, about 12 minutes using VF2
    rangeBegin = 4;
    rangeEnd = 5;
    dropRootEdges = true;
    for (int i = rangeBegin; i <= rangeEnd; i++) {
//      testTemplate(getNodes((byte) (i - 2), (byte) 2));
    }
  }

  /**
   * Isomorph-Free graphs with 2 root nodes
   */
  public void testIsomorphFreeR2() {

    reset();
    printTest("Isomorph-Free, 2 root nodes");
// for N=7, about 12 minutes using VF2
    int rootNodes = 3;
    rangeBegin = 4;
    rangeEnd = 6;
    dropRootEdges = true;
    isomorphFree = true;
    for (int i = rangeBegin; i <= rangeEnd; i++) {
     // minEdges = 0;
// maxEdges = (i - rootNodes) * rootNodes;
     // maxEdges = (i - rootNodes) * (i - rootNodes - 1) / 4 - (((i - rootNodes) * (i - rootNodes - 1) % 4 == 0) ? 1 : 0);
  //    testTemplate(getNodes((byte) (i - rootNodes), (byte) rootNodes));
    }
    System.out.println("Called: " + Match.called);
    System.out.println("Pre-Filtered: " + (Match.total - Match.called));
  }

  public void testIsomorphFreeConnected() {

    reset();
    printTest("Isomorph-Free, Connected");
    // for N=7, about 12 minutes using VF2
    // rangeBegin = 7;
    // rangeEnd = 7;
    enumConnected = true;
    isomorphFree = true;
    for (int i = rangeBegin; i <= rangeEnd; i++) {
// testTemplate(getNodes((byte) i));
    }
  }

  public void testOptimizedIsomorphFree() {

    reset();
    printTest("Isomorph-Free (Optimized)");
 rangeBegin = 2;
 rangeEnd = 6;
    isomorphFree = true;
    for (int i = rangeBegin; i <= rangeEnd; i++) {
//      int numEdges = i * (i - 1) / 2;
//      minEdges = 0;
//      maxEdges = numEdges / 2;
// testTemplate(getNodes((byte) i));
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