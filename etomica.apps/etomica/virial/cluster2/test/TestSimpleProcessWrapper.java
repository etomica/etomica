/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.test;

import java.io.BufferedReader;
import java.io.IOException;

import etomica.virial.cluster2.nauty.NautyInfo;
import etomica.virial.cluster2.nauty.ProcessWrapper;
import etomica.virial.cluster2.nauty.impl.SimpleProcessWrapper;

public class TestSimpleProcessWrapper extends CustomTestCase {

  private static final String NAUTY_PATH = "/opt/nauty/nauty24b7/virial-customized/";

  public void testTemplate(final NautyInfo info) {

    System.out.println("============");
    System.out.println(info.toString());
    ProcessWrapper pw = new SimpleProcessWrapper(info);
    try {
      runGC();
      long time1 = System.nanoTime();
      pw.run();
      enumerated = 0;
      if (printPermutations) {
        System.out.println("============");
      }
      String line;
      BufferedReader br = new BufferedReader(pw.getProcessOutput());
      while ((line = br.readLine()) != null) {
        if (printPermutations) {
          System.out.println(line);
        }
        enumerated++;
      }
      if (printPermutations) {
        System.out.println("============");
      }
      enumerated /= 2;
      elapsed = (System.nanoTime() - time1) / K;
      printEnumerated(info.getNumBlackNodes());
      printRuntime();
      printMemory();
// printPermutations();
      System.out.println("============");
      System.out.println();
    }
    catch (IOException e) {
      e.printStackTrace();
      fail("IOException");
    }
  }

  boolean enumConnected = false;
  boolean enumBiconnected = true;
  boolean useUpperTriangle = true;

  @Override
  public void setUp() {

    super.setUp();
    enumConnected = false;
    enumBiconnected = false;
    useUpperTriangle = false;
  }
  
  protected NautyInfo getNautyInfo(int numNodes) {

    NautyInfo result = new NautyInfo(NAUTY_PATH, numNodes);
    result.setConnected(enumConnected);
    result.setBiconnected(enumBiconnected);
    result.setUpperTriangle(useUpperTriangle);
    return result;
  }

  public void testGeneral() {

    testTemplate(getNautyInfo(4));
    testTemplate(getNautyInfo(5));
    testTemplate(getNautyInfo(6));
    testTemplate(getNautyInfo(7));
  }

  public void testConnected() {

    enumConnected = true;
    testGeneral();
  }

  public void testBiconnected() {

    enumBiconnected = true;
    testGeneral();
  }
}