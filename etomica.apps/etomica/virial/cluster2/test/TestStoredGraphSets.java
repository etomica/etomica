/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.EdgesFilter;
import etomica.virial.cluster2.graph.EdgesSetVisitor;
import etomica.virial.cluster2.graph.FilterFactory;
import etomica.virial.cluster2.graph.GraphFactory;
import etomica.virial.cluster2.graph.GraphSet;
import etomica.virial.cluster2.graph.GraphSetIO;
import etomica.virial.cluster2.graph.Nodes;

public class TestStoredGraphSets extends CustomTestCase {

  private boolean enumConnected = false;
  private boolean enumBiconnected = false;
  private boolean nullFiltered = false;
  private int rangeEnd = 0;
  private int rangeStart = 0;
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
      family = GraphFactory.storedGraphSet(nodes, filter);
      elapsed = (System.nanoTime() - time1) / K; // µs
      runGC();
      memoryUse();
      enumerated = family.getEdgesSet().size();
//      family.visitEdgesSet(nautyVisitor);
      System.out.println(family.getTags());
      printEnumerated(numNodes);
      printRuntime();
      printMemory();
      if (printPermutations) {
        System.out.println();
        family.visitEdgesSet(printVisitor);
      }
      System.out.println();
    }
    catch (RuntimeException ge) {
      ge.printStackTrace();
      fail("Unexpected exception: " + ge.getStackTrace());
    }
  }

  @Override
  public void setUp() {

    super.setUp();
    enumConnected = false;
    enumBiconnected = false;
    nullFiltered = false;
    printMemory = true;
    printPermutations = true;
    rangeStart = 2;
    rangeEnd = 6;
  }

  private void write() {

    rangeStart = 2;
    rangeEnd = 6;
    for (int i = rangeStart; i <= rangeEnd; i++) {
      testTemplate(i);
      try {
        File f = new File("/home/demian/MyDownloads/graphset-" + i);
        if (f.exists()) {
          f.delete();
        }
        f.createNewFile();
        BufferedWriter w = new BufferedWriter(new FileWriter(f));
        GraphSetIO.writeAsText(family, w);
      }
      catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
    }
  }

  private void read() {

    rangeStart = 2;
    rangeEnd = 6;
    for (int i = rangeStart; i <= rangeEnd; i++) {
//      testTemplate(i);
      try {
        File f = new File("/home/demian/MyDownloads/graphset-" + i);
        if (!f.exists()) {
          continue;
        }
        index = 0;
        BufferedReader r = new BufferedReader(new FileReader(f));
        runGC();
        long time1 = System.nanoTime();
        // create family
        GraphFactory.DYNAMIC_ALLOCATION = enumerated < M;
        family = GraphSetIO.readAsText(r);
        elapsed = (System.nanoTime() - time1) / K; // µs
        runGC();
        memoryUse();
        enumerated = family.getEdgesSet().size();
//        family.visitEdgesSet(nautyVisitor);
        System.out.println(family.getTags());
        printEnumerated(family.getNodes().count());
        printRuntime();
        printMemory();
        if (printPermutations) {
          System.out.println();
          family.visitEdgesSet(printVisitor);
        }
        System.out.println();
      }
      catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
    }
  }

  public void testWrite() {
//    write();
  }
  
  public void testRead() {
    read();
  }
  //  public void testConnected() {
//
//    enumConnected = true;
//    testGeneral();
//  }
//
//  public void testBiconnected() {
//
//    enumBiconnected = true;
//    testGeneral();
//  }
//
//  public void testNullFilteredGeneral() {
//
//    nullFiltered = true;
//    testGeneral();
//  }
//
//  public void testNullFilteredConnected() {
//
//    nullFiltered = true;
//    enumConnected = true;
//    testGeneral();
//  }
//
//  public void testNullFilteredBiconnected() {
//
//    nullFiltered = true;
//    enumBiconnected = true;
//    testGeneral();
//  }
}
