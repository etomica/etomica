package etomica.virial.cluster2.graph.test;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.EdgesFilter;
import etomica.virial.cluster2.graph.EdgesSetVisitor;
import etomica.virial.cluster2.graph.GraphFactory;
import etomica.virial.cluster2.graph.GraphSet;
import etomica.virial.cluster2.graph.Nodes;
import etomica.virial.cluster2.nauty.NautyInfo;
import etomica.virial.cluster2.nauty.ProcessWrapper;
import etomica.virial.cluster2.nauty.impl.SimpleProcessWrapper;
import etomica.virial.cluster2.test.CustomTestCase;

public class TestNautyEdgesGenerator extends CustomTestCase {

  private static final String NAUTY_PATH = "/opt/nauty/nauty24b7/test/";
  private boolean enumConnected = false;
  private boolean enumBiconnected = false;
  private boolean nullFiltered = false;
  private GraphSet family;
  private double isomorphCount = 0;
  private EdgesSetVisitor nautyVisitor = new EdgesSetVisitor() {

    @Override
    public boolean visit(Edges edges) {

      isomorphCount += edges.getMetadata().getCoefficient();
      return true;
    }
  };

  private EdgesSetVisitor printVisitor = new EdgesSetVisitor() {

    private int index = 0;
    
    @Override
    public boolean visit(Edges edges) {

      System.out.println(index++ + ": " + edges.toString());
      return true;
    }
  };

  private void testNautyFamilyN(int numNodes, boolean nullFiltered) {

    long time1;
    isomorphCount = 0;
    printMemory = true;
    printPermutations = true;
    permutations = "";
    expected = 1 << (numNodes) * (numNodes - 1) / 2;
    enumerated = 1 << (numNodes) * (numNodes - 1) / 2;
    if (nullFiltered) {
      expected = 0;
    }
    try {
      // create nodes
      Nodes nodes = GraphFactory.defaultNodes((byte) numNodes);
      // create generator
      EdgesFilter filter = null;
      if (nullFiltered) {
        filter = GraphFactory.nullFilter(null);
      }
      runGC();
      time1 = System.nanoTime();
      // create family
      GraphFactory.DYNAMIC_ALLOCATION = enumerated < M;
      family = GraphFactory.nautyGraphSet(nodes, filter, getNautyProcess(numNodes));
      elapsed = (System.nanoTime() - time1) / K; // µs
      runGC();
      memoryUse();
      System.out.println(getNautyInfo(numNodes).toString());
      printEnumerated(numNodes);
      printRuntime();
      printMemory();
      if (printPermutations) {
        System.out.println();
        family.visitEdgesSet(printVisitor);
//        permutations = family.toString();
//        printPermutations();
      }
      System.out.println();
      family.visitEdgesSet(nautyVisitor);
      assertEquals(expected, isomorphCount, 0.0001);
    } catch (RuntimeException ge) {
      fail("Unexpected exception: " + ge.getStackTrace());
    }
  }

  private ProcessWrapper getNautyProcess(final int numNodes) {

    return new SimpleProcessWrapper(getNautyInfo(numNodes));
  }

  protected NautyInfo getNautyInfo(final int numNodes) {

    NautyInfo result = new NautyInfo(NAUTY_PATH, numNodes);
    result.setConnected(enumConnected && !enumBiconnected);
    result.setBiconnected(enumBiconnected);
    printPermutations = true;
    return result;
  }

  public void testNautyFamily0() {

//    try {
//      family = new SimpleGraphSet((byte) 0);
//      fail("Invalid node count");
//    } catch (GraphException ge) {
//    }
  }

  public void testNautyFamily1() {

  //  testNautyFamilyN(1, nullFiltered);
  }

  public void testNautyFamily2() {

  //  testNautyFamilyN(2, nullFiltered);
  }

  public void testNautyFamily3() {

  //  testNautyFamilyN(3, nullFiltered);
  }

  public void testNautyFamily4() {

  //  testNautyFamilyN(4, nullFiltered);
  }

  public void testNautyFamily5() {

  //  testNautyFamilyN(5, nullFiltered);
  }

  public void testNautyFamily6() {

  //  testNautyFamilyN(6, nullFiltered);
  }

  public void testNautyFamily7() {

 //   testNautyFamilyN(7, nullFiltered);
  }

  public void testNautyFamily8() {

 //   testNautyFamilyN(8, nullFiltered);
  }

  public void testNautyFamily9() {

    testNautyFamilyN(9, nullFiltered);
  }
}