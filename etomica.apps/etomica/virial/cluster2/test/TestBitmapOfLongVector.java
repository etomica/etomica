/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.test;

import junit.framework.TestCase;
import etomica.virial.cluster2.bitmap.impl.BitmapOfLongVector;

public class TestBitmapOfLongVector extends TestCase {

  private int bmSize = 80;
  private BitmapOfLongVector bm1 = null;
  private BitmapOfLongVector bm2 = null;
  private BitmapOfLongVector bm3 = null;
  private BitmapOfLongVector bm4 = null;
  private BitmapOfLongVector bm5 = null;
  private BitmapOfLongVector bm6 = null;
  private String bmStr1 = "1110011001010100000011001101001100111001110000000000110100001111110100110011101";
  private String bmStr2 = "1111111111111111";
  private String bmStr3 = "0000000000000000";
  private String bmStr4 = "0000000000000001";
  private String bmStr5 = "00110010"; // 50
  private String bmStr6 = "00000000"; // 00

  @Override
  public void setUp() {

    bm1 = new BitmapOfLongVector(bmStr1);
    bm2 = new BitmapOfLongVector(bmStr2);
    bm3 = new BitmapOfLongVector(bmStr3);
    bm4 = new BitmapOfLongVector(bmStr4);
    bm5 = new BitmapOfLongVector(bmStr5);
    bm6 = new BitmapOfLongVector(bmStr6);
  }

  private void println(String val) {

    System.out.println(val);
  }

  public void testBitmapInt() {

    BitmapOfLongVector bm = new BitmapOfLongVector(bmSize);
    assertEquals(bmSize, bm.toString().length());
  }

  public void testBitmapIntBoolean1() {

    BitmapOfLongVector bm = new BitmapOfLongVector(bmSize, false);
    assertEquals(bmSize, bm.toString().length());
  }

  public void testBitmapIntBoolean2() {

    BitmapOfLongVector bm = new BitmapOfLongVector(bmSize, true);
    assertEquals(bmSize, bm.toString().length());
  }

  public void testBitmapString() {

    assertEquals(bmStr1, bm1.toString());
  }

  public void testClone() {

    BitmapOfLongVector bm = (BitmapOfLongVector) bm1.copy();
    assertEquals(bm1, bm);
  }

  public void testInc() {

    bm2.inc();
    println(bm2.toString());
    assertEquals(bm3, bm2);
    bm2.inc();
    println(bm2.toString());
    assertEquals(bm4, bm2);
    bm5.not();
    println(bm5.toString());
    for (int i = 0; i < 51; i++) {
      bm5.inc();
      println(bm5.toString());
    }
    assertEquals(bm6, bm5);
  }

  public void testDec() {

    bm4.dec();
    println(bm4.toString());
    assertEquals(bm3, bm4);
    bm4.dec();
    println(bm4.toString());
    assertEquals(bm2, bm4);
    println(bm5.toString());
    for (int i = 0; i < 50; i++) {
      bm5.dec();
      println(bm5.toString());
    }
    assertEquals(bm6, bm5);
  }

  public void testNot() {

    BitmapOfLongVector bmc = (BitmapOfLongVector) bm2.copy();
    bm2.not();
    println(bm2.toString());
    assertEquals(bm3, bm2);
    bm2.not();
    println(bm2.toString());
    assertEquals(bmc, bm2);

    bmc = (BitmapOfLongVector) bm1.copy();
    bm1.not();
    println(bm1.toString());
    bm1.not();
    println(bm1.toString());
    assertEquals(bmc, bm1);
  }

  public void testAnd() {
    fail("Not yet implemented");
  }

  public void testNand() {
    fail("Not yet implemented");
  }

  public void testOr() {
    fail("Not yet implemented");
  }

  public void testXor() {
    fail("Not yet implemented");
  }

  public void testCompareTo() {
    fail("Not yet implemented");
  }

  public void testTestBit() {
    fail("Not yet implemented");
  }

  public void testBitSize() {
    fail("Not yet implemented");
  }

  public void testHighestBit() {
    fail("Not yet implemented");
  }

  public void testLowestBit() {
    fail("Not yet implemented");
  }

  public void testSize() {
    fail("Not yet implemented");
  }
}