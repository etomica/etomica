/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.test;

import junit.framework.TestCase;
import etomica.virial.cluster2.bitmap.impl.BitmapOfLong;

public class TestBitmapOfLong extends TestCase {

  private int bmSize = 80;
  private BitmapOfLong bm1 = null;
  private BitmapOfLong bm2 = null;
  private BitmapOfLong bm3 = null;
  private BitmapOfLong bm4 = null;
  private BitmapOfLong bm5 = null;
  private BitmapOfLong bm6 = null;
  private String bmStr1 = "1110011001010100000011001101001100111001110000000000110100";
  private String bmStr2 = "1111111111111111";
  private String bmStr3 = "0000000000000000";
  private String bmStr4 = "0000000000000001";
  private String bmStr5 = "00110010"; // 50
  private String bmStr6 = "00000000"; // 00

  @Override
  public void setUp() {

    bm1 = new BitmapOfLong(bmStr1);
    bm2 = new BitmapOfLong(bmStr2);
    bm3 = new BitmapOfLong(bmStr3);
    bm4 = new BitmapOfLong(bmStr4);
    bm5 = new BitmapOfLong(bmStr5);
    bm6 = new BitmapOfLong(bmStr6);
  }

  private void println(String val) {

    System.out.println(val);
  }

  public void testBitmapInt() {

    BitmapOfLong bm = new BitmapOfLong(bmSize);
    assertEquals(bmSize, bm.toString().length());
  }

  public void testBitmapIntBoolean1() {

    BitmapOfLong bm = new BitmapOfLong(bmSize, false);
    assertEquals(bmSize, bm.toString().length());
  }

  public void testBitmapIntBoolean2() {

    BitmapOfLong bm = new BitmapOfLong(bmSize, true);
    assertEquals(bmSize, bm.toString().length());
  }

  public void testBitmapString() {

    assertEquals(bmStr1, bm1.toString());
  }

  public void testClone() {

    BitmapOfLong bm = (BitmapOfLong) bm1.copy();
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

    BitmapOfLong bmc = (BitmapOfLong) bm2.copy();
    bm2.not();
    println(bm2.toString());
    assertEquals(bm3, bm2);
    bm2.not();
    println(bm2.toString());
    assertEquals(bmc, bm2);

    bmc = (BitmapOfLong) bm1.copy();
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
