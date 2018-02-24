/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph;

import junit.framework.TestCase;

public class CustomTestCase {

  private Runtime rt = Runtime.getRuntime();
  protected static int K = 1000;
  protected static int M = K * K;
  protected static int G = M * K;
  protected static int T = G * K;
  protected static int KB = 1024;
  protected static int MB = KB * KB;
  protected long mem1 = 0;
  protected long mem2 = 0;
  protected long elapsed = 0;
  protected long enumerated = 0;
  protected long expected = 0;
  protected String permutations = "";
  protected String memoryUse = "";
  protected boolean printMemory = false;
  protected boolean printPermutations = false;
  protected boolean printRuntime = true;
  protected boolean checkAssertion = true;

  private long usedMemory() {

    return rt.totalMemory() - rt.freeMemory();
  }

  protected void runGC() throws RuntimeException {

    // It helps to call Runtime.gc() using several method calls:
    for (int r = 0; r < 6; ++r)
      _runGC();
  }

  private void _runGC() throws RuntimeException {

    long usedMem1 = usedMemory(), usedMem2 = Long.MAX_VALUE;
    for (int i = 0; (usedMem1 < usedMem2) && (i < 500); ++i) {
      rt.runFinalization();
      rt.gc();
      Thread.yield();
      usedMem2 = usedMem1;
      usedMem1 = usedMemory();
    }
  }

  protected void printEnumerated(int numNodes) {

    System.out.println("nodes....: " + numNodes);
    if (enumerated > 100*M) {
      System.out.println("graphs...: " + (enumerated / M) + "M");
    }
    else if (enumerated > 100*K) {
      System.out.println("graphs...: " + (enumerated / K) + "K");
    }
    else {
      System.out.println("graphs...: " + enumerated);
    }
  }

  protected void printPermutations() {

    if (printPermutations) {
      System.out.println(permutations);
    }
  }

  protected void printRuntime() {

    if (printRuntime) {
      if (elapsed < K) {
        System.out.print("runtime..: " + (elapsed) + "µs (");
      }
      else if (elapsed < M) {
        System.out.print("runtime..: " + (elapsed / K) + "ms (");
      }
      else if (elapsed < G) {
        System.out.print("runtime..: " + (elapsed / M) + "sec (");
      }
      else if (elapsed < (60 * G)) {
        System.out.print("runtime..: " + (elapsed / (60 * M)) + "min (");
      }
      else if (elapsed < (60 * 60 * G)) {
        System.out.print("runtime..: " + (elapsed / (60 * 60 * M)) + "h (");
      }
      if (enumerated > 0) {
        if ((elapsed / enumerated) < K) {
          System.out.print((elapsed / enumerated) + "µs/graph)");
        }
        else {
          System.out.print(((elapsed / K) / enumerated) + "ms/graph)");
        }
        if (enumerated < 10*M) {
          System.out.println(" [" + enumerated / K + "K graphs]");
        }
        else {
          System.out.println(" [" + enumerated / M + "M graphs]");
        }
      }
      else {
        System.out.println(")");
      }
    }
  }

  protected void memoryUse() {

    if (printMemory) {
      memoryUse = "";
      mem1 = usedMemory();
      runGC();
      mem2 = usedMemory();
      if (mem1 > MB) {
        memoryUse += "memory...: " + (mem1 / MB) + "MB";
      }
      else {
        memoryUse += "memory...: " + (mem1 / KB) + "KB";
      }
      if (expected != 0) {
        if (expected > (mem1 / KB)) {
          memoryUse += " (" + (mem1 / expected) + "B/graph)\n";
        }
        else {
          memoryUse += " (" + (mem1 / expected / KB) + "KB/graph)\n";
        }
      }
      else {
        memoryUse += "\n";
      }
      if ((mem1 - mem2) > MB) {
        memoryUse += "reclaim..: " + ((mem1 - mem2) / MB) + "MB\n";
      }
      else {
        memoryUse += "reclaim..: " + ((mem1 - mem2) / KB) + "KB\n";
      }
    }
  }

  public void setUp() {

    reset();
  }

  protected void printMemory() {

    if (printMemory) {
      System.out.println(memoryUse);
    }
  }

  protected void printTest(String testName) {

    System.out.println();
    System.out.println("************************");
    System.out.println(testName);
    System.out.println("========================");
  }

  public void reset() {

    printMemory = false;
    printRuntime = true;
    printPermutations = false;
    checkAssertion = true;
  }

  // Asymptotic Object overhead: 4 bytes/object.
  //
  // public void testObjectMemoryFootprint() {
  //
  // Object[] foo = new Object[100 * M];
  // for (int i = 0; i < 2 * M; i++) {
  // foo[i] = new Object();
  // }
  // expected = 100 * M;
  // printMemory = true;
  // memoryUse();
  // printMemory();
  // }
  // Asymptotic SimpleEdges overhead: 4 bytes/object.
  //
  // public void testObjectMemoryFootprint() {
  //
  // Object[] foo = new Bitmap[10 * M];
  // for (int i = 0; i < 2 * M; i++) {
  // foo[i] = BitmapFactory.getBitmap(50,false);
  // }
  // if (Math.random() > 15) {
  // foo[15] = null;
  // }
  // expected = 100 * M;
  // printMemory = true;
  // memoryUse();
  // printMemory();
  // }
}
