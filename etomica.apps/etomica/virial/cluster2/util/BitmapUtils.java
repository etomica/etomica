/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.util;

/*
 * BitmapUtils are used across all algorithm implementations below.
 * @author Demian Lessa
 */
public class BitmapUtils {

  protected static final int SZ_INT = 32;

  public static int bitMask(int bitSize) {

    return (0xFFFFFFFF >>> (SZ_INT - bitSize));
  }

  public static int bitOnMask(int bit) {

    return (0x00000001 << bit);
  }

  public static int bitOffMask(int bit) {

    return ~bitOnMask(bit);
  }

  public static int leftmostBit(int map) {

    for (int i = 0; i < SZ_INT; i++) {
      if ((map & bitOnMask((byte) i)) > 0) {
        return i;
      }
    }
    return -1;
  }
}
