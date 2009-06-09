package etomica.virial.cluster2.nauty.test;

import java.io.BufferedReader;
import java.io.IOException;

import etomica.virial.cluster2.nauty.NautyInfo;
import etomica.virial.cluster2.nauty.ProcessWrapper;
import etomica.virial.cluster2.nauty.impl.SimpleProcessWrapper;
import etomica.virial.cluster2.test.CustomTestCase;

public class TestSimpleProcessWrapper extends CustomTestCase {

  private static final String NAUTY_PATH = "/opt/nauty/nauty24b7/test/";

  public void testTemplate(final NautyInfo info) {

    System.out.println(info.toString());
    ProcessWrapper pw = new SimpleProcessWrapper(info);
    try {
      runGC();
      long time1 = System.nanoTime();
      pw.run();
      elapsed = (System.nanoTime() - time1) / K;
//      if (enumGraphs) {
//        enumerated = pw.getStdout().size() / (info.getNumBlackNodes() + 1);
//      }
//      else {
//        StringTokenizer st = new StringTokenizer(pw.getStdout().get(1));
//        st.nextToken();
//        enumerated = new Integer(st.nextToken());
//      }
      printEnumerated(info.getNumBlackNodes());
      printRuntime();
      printMemory();
      if (printPermutations) {
//        permutations = pw.getStdout().toString();
        String line;
        BufferedReader br = new BufferedReader(pw.getProcessOutput());
        while ((line = br.readLine()) != null) {
        System.out.println(line);
        }
      }
//      printPermutations();
      System.out.println();
    } 
    catch (IOException e) {
      e.printStackTrace();
      fail("IOException");
    }
  }

  boolean enumConnected = false; 
  boolean enumBiconnected = true; 

  protected NautyInfo getNautyInfo(final int numNodes) {
    
    NautyInfo result = new NautyInfo(NAUTY_PATH, numNodes);
    result.setConnected(enumConnected);
    result.setBiconnected(enumBiconnected);
    printPermutations = true;
    return result;
  }
  
  public void testGetOutput1() {

    testTemplate(getNautyInfo(4));
  }
  
  public void testGetOutput2() {

//    testTemplate(getNautyInfo(6));
  }
  
  public void testGetOutput3() {

//    testTemplate(getNautyInfo(7));
  }
  
  public void testGetOutput4() {

//    testTemplate(getNautyInfo(8));
  }

  public void testGetOutput5() {

//    testTemplate(getNautyInfo(9));
  }

  public void testGetOutput6() {

//    testTemplate(getNautyInfo(10));
  }

  public void testGetOutput7() {

//    testTemplate(getNautyInfo(11));
  }
}