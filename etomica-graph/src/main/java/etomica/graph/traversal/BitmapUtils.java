/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.traversal;

/*
 * BitmapUtils are used by the traversal algorithms.
 *
 * @author Demian Lessa
 */
public class BitmapUtils {

  protected static final byte SZ_INT = 32;

  public static int bitMask(byte bitSize) {

    return (0xFFFFFFFF >>> (SZ_INT - bitSize));
  }

  public static int bitOnMask(byte bit) {

    return (0x00000001 << bit);
  }

  public static int bitOffMask(byte bit) {

    return ~bitOnMask(bit);
  }

  public static byte leftmostBit(int map) {

    for (byte i = 0; i < SZ_INT; i++) {
      if ((map & bitOnMask(i)) > 0) {
        return i;
      }
    }
    return -1;
  }
}
