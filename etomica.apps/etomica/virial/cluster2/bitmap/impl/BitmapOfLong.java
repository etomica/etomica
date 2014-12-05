/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.bitmap.impl;

import etomica.virial.cluster2.bitmap.AbstractBitwiseBitmap;
import etomica.virial.cluster2.bitmap.Bitmap;

/**
 * This class is the most space efficient encoding of bitmaps, however it can
 * only encode up to 64 bits. Assume N is the number of nodes in an undirected
 * graph. The complete graph with N nodes has N(N-1)/2 edges. Hence, this class
 * can encode all graphs with N = 1..11. NAUTY, the most efficient graph
 * enumeration program as of today (May, 2009) can enumerate graphs with
 * non-isomorphs for N = 11 in 12 minutes (in an INTEL CORE2 DUO T9400 2.53GHz
 * with 4GB RAM). For N=12, it takes approximately 33h, so it basically falls
 * out of the online realm. Further, for N=12, the maximum number of edges is
 * 66, hence, we cannot use this class to encode the graphs anyway. This is a
 * rather convenient coincidence.
 * 
 * @author Demian Lessa
 */
public class BitmapOfLong extends AbstractBitwiseBitmap {

  private byte bitSize = 0;
  private long bitmap = 0;

  // ***********************
  // * PUBLIC CONSTRUCTORS
  // ***********************
  public BitmapOfLong(final Bitmap other) {

    this(other.bitSize(), false);
    copyFrom(other);
  }

  public BitmapOfLong(final int capacity) {

    this(capacity, false);
  }

  public BitmapOfLong(final int capacity, final boolean isSet) {

    setBitSize(capacity);
    allocateBitmap();
    setBits(isSet);
  }

  public BitmapOfLong(final String strBitmap) {

    this((byte) strBitmap.length(), false);
    copyFrom(strBitmap);
  }

  // ***********************
  // * PROTECTED CONSTRUCTOR
  // ***********************
  protected BitmapOfLong() {

    setBitSize((byte) 0);
  }

  // ***********************
  // * PUBLIC METHODS
  // ***********************
  /**
   * @TIME = O(1) for a compatible instance of other.
   */
  @Override
  public void and(final Bitmap other) {

    if (other instanceof BitmapOfLong) {
      BitmapOfLong bm = (BitmapOfLong) other;
      bitmap &= bm.bitmap;
    }
    else {
      super.and(other);
    }
  }

  public int bitSize() {

    return bitSize;
  }

  /**
   * The highest order bit is at bit 0 of bitmap, and the lowest order bit is at
   * bit bitSize() - 1 of bitmap.
   * 
   * @TIME = O(1).
   */
  public void clearBit(final int bitIndex) {

    bitmap &= ~maskSingleBit(bitIndex);
  }

  /**
   * @TIME = O(1) for an instance of BitmapOfLong.
   */
  @Override
  public boolean equals(Object other) {

    if (other instanceof BitmapOfLong) {
      BitmapOfLong bm = (BitmapOfLong) other;
      return ((bitmap & maskSingleLong()) == (bm.bitmap & bm.maskSingleLong()));
    }
    else {
      return super.equals(other);
    }
  }

  /**
   * @TIME = O(1).
   */
  @Override
  public int hashCode() {

    return (new Long(bitmap)).hashCode();
  }

  /**
   * @TIME = O(1) for a compatible instance.
   */
  @Override
  public void nand(final Bitmap other) {

    if (other instanceof BitmapOfLong) {
      BitmapOfLong bm = (BitmapOfLong) other;
      bitmap &= ~bm.bitmap;
    }
    else {
      super.nand(other);
    }
  }

  /**
   * @TIME = O(1).
   */
  @Override
  public void not() {

    bitmap = ~bitmap;
  }

  /**
   * @TIME = O(1) for a compatible instance.
   */
  @Override
  public void or(final Bitmap other) {

    if (other instanceof BitmapOfLong) {
      BitmapOfLong bm = (BitmapOfLong) other;
      bitmap |= bm.bitmap;
    }
    else {
      super.or(other);
    }
  }

  /**
   * The highest order bit is at bit 0 of bitmap, and the lowest order bit is at
   * bit bitSize() - 1 of bitmap.
   * 
   * @TIME = O(1).
   */
  public void setBit(final int bitIndex) {

    bitmap |= (maskSingleLong() & maskSingleBit(bitIndex));
  }

  /**
   * @TIME = O(1).
   */
  @Override
  public void setBits(final boolean value) {

    bitmap = value ? maskSingleLong() : Bitmap.LONG_ZERO;
  }

  /**
   * @TIME = O(1).
   */
  public boolean testBit(final int bitIndex) {

    long bm = maskSingleBit(bitIndex);
    return (bitmap & bm) == bm;
  }

  /**
   * @TIME = O(1) for a compatible instance.
   */
  @Override
  public void xor(final Bitmap other) {

    if (other instanceof BitmapOfLong) {
      BitmapOfLong bm = (BitmapOfLong) other;
      bitmap ^= bm.bitmap;
    }
    else {
      super.xor(other);
    }
  }

  // ***********************
  // * PROTECTED METHODS
  // ***********************
  @Override
  protected void allocateBitmap() {

    bitmap = 0;
  }

  /**
   * @TIME = O(1).
   */
  @Override
  protected void copyFrom(final Bitmap other) {

    if (other instanceof BitmapOfLong) {
      BitmapOfLong copy = (BitmapOfLong) other;
      bitmap = copy.bitmap;
    }
  }

  @Override
  protected Bitmap createInstance(final Bitmap other) {

    return new BitmapOfLong(other);
  }

  @Override
  protected Bitmap createInstance(final int capacity) {

    return new BitmapOfLong(capacity);
  }

  @Override
  protected Bitmap createInstance(final int capacity, final boolean isSet) {

    return new BitmapOfLong(capacity, isSet);
  }

  @Override
  protected Bitmap createInstance(final String strBitmap) {

    return new BitmapOfLong(strBitmap);
  }

  /**
   * Returns a long with a single bit set, corresponding to the bit offset of
   * the given bit index.
   * 
   * @TIME = O(1).
   */
  protected long maskSingleBit(final int bitIndex) {

    return Bitmap.LONG_ONE << (bitSize() - 1 - bitIndex);
  }

  /**
   * Returns a bitmask corresponding to the valid bits of the long with the
   * given index.
   * 
   * @TIME = O(1).
   */
  protected long maskSingleLong() {

    return Bitmap.LONG_MASK >>> (Bitmap.SZ_LONG - bitSize());
  }

  protected void setBitSize(final int capacity) {

    bitSize = (byte) capacity;
  }
}