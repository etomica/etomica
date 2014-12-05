/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.test;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.EdgesFilter;
import etomica.virial.cluster2.graph.EdgesSetVisitor;
import etomica.virial.cluster2.graph.FilterFactory;
import etomica.virial.cluster2.graph.GraphFactory;
import etomica.virial.cluster2.graph.GraphSet;
import etomica.virial.cluster2.graph.Nodes;
import etomica.virial.cluster2.nauty.NautyInfo;
import etomica.virial.cluster2.nauty.ProcessWrapper;
import etomica.virial.cluster2.nauty.impl.SimpleProcessWrapper;

public class TestNautyEdgesGenerator extends CustomTestCase {

  private static final String NAUTY_PATH = "/opt/nauty/nauty24b7/virial-customized/";
  private boolean enumConnected = false;
  private boolean enumBiconnected = false;
  private boolean nullFiltered = false;
  private boolean useUpperTriangle = true;
  private int rangeEnd = 0;
  private GraphSet family;
  private int index = 0;
  private double isomorphCount = 0;
  private FilterFactory ffactory = new FilterFactory();
  private EdgesSetVisitor nautyVisitor = new EdgesSetVisitor() {

    public boolean visit(Edges edges) {

      enumerated++;
      isomorphCount += edges.getMetadata().getCoefficient().getValue1();
      return true;
    }
  };
  private EdgesSetVisitor printVisitor = new EdgesSetVisitor() {

    public boolean visit(Edges edges) {

      System.out.println(index++ + ": " + edges.toString());
      return true;
    }
  };

  private void testTemplate(int numNodes) {

    long time1;
    index = 0;
    isomorphCount = 0;
    permutations = "";
    expected = 1 << (numNodes) * (numNodes - 1) / 2;
    enumerated = 0;
    if (nullFiltered) {
      expected = 0;
    }
    try {
      // create nodes
      Nodes nodes = GraphFactory.defaultNodes((byte) numNodes);
      // create generator
      // pass-through filter
      EdgesFilter filter = ffactory.trueFilter();
      if (nullFiltered) {
        filter.chain(ffactory.falseFilter());
      }
      runGC();
      time1 = System.nanoTime();
      // create family
      GraphFactory.DYNAMIC_ALLOCATION = enumerated < M;
      family = GraphFactory.nautyGraphSet(nodes, filter,
          getNautyProcess(numNodes));
      elapsed = (System.nanoTime() - time1) / K; // µs
      runGC();
      memoryUse();
      family.visitEdgesSet(nautyVisitor);
      System.out.println(family.getTags());
      System.out.println(getNautyInfo(numNodes).toString());
      printEnumerated(numNodes);
      printRuntime();
      printMemory();
      if (printPermutations) {
        System.out.println();
      }
      family.visitEdgesSet(printVisitor);
      System.out.println();
      if (!enumConnected && !enumBiconnected) {
        assertEquals(expected, isomorphCount, 0.0001);
      }
    }
    catch (RuntimeException ge) {
      fail("Unexpected exception: " + ge.getStackTrace());
    }
  }

  private ProcessWrapper getNautyProcess(int numNodes) {

    return new SimpleProcessWrapper(getNautyInfo(numNodes));
  }

  protected NautyInfo getNautyInfo(int numNodes) {

    NautyInfo result = new NautyInfo(NAUTY_PATH, numNodes);
    result.setConnected(enumConnected && !enumBiconnected);
    result.setBiconnected(enumBiconnected);
    result.setUpperTriangle(useUpperTriangle);
    return result;
  }

  @Override
  public void setUp() {

    super.setUp();
    enumConnected = false;
    enumBiconnected = false;
    nullFiltered = false;
    printMemory = true;
    printPermutations = true;
    useUpperTriangle = true;
    rangeEnd = 5;
  }

  public void testEmptyFamily() {

    try {
      family = GraphFactory.nautyGraphSet(GraphFactory.defaultNodes((byte) 0),
          null, getNautyProcess(0));
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

  public void testBiconnected() {

    enumBiconnected = true;
    testGeneral();
  }

  public void testNullFilteredGeneral() {

    nullFiltered = true;
    testGeneral();
  }

  public void testNullFilteredConnected() {

    nullFiltered = true;
    enumConnected = true;
    testGeneral();
  }

  public void testNullFilteredBiconnected() {

    nullFiltered = true;
    enumBiconnected = true;
    testGeneral();
  }
}