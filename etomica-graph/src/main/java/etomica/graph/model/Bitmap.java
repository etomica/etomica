/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.model;

/**
 * Abstracts away the form in which bitmaps implement their operations and store their
 * data.
 *
 * @author Demian Lessa
 */
public interface Bitmap extends Comparable<Bitmap> {

  public final static char CH_ONE = '1';
  public final static char CH_ZERO = '0';
  public final static long LONG_MASK = 0xffffffffffffffffL;
  public final static long LONG_ONE = 0x0000000000000001L;
  public final static long LONG_ZERO = 0x0000000000000000L;
  public final static byte SZ_LONG = 64;

  /**
   * and: bitmap -> bitmap = (bitmap & other). Assumes the bit length of the other bitmap
   * is AT LEAST as large as the bit length of this bitmap.
   */
  public void and(final Bitmap other);

  /**
   * Returns the number of bits set in the bitmap.
   */
  public int bitCount();

  /**
   * Returns the length in bits of the bitmap.
   */
  public int bitSize();

  public void clearBit(final int bitIndex);

  /**
   * Returns a bitmap instance with the same number of bits as this bitmap. If other has
   * fewer bits than this bitmap, the output bitmap is 0-padded and has the same integer
   * value as other. If other has more bits than this bitmap, the output bitmap is the
   * result of truncating the highest bits of other.
   *
   * @return copy of the this bitmap
   */
  public Bitmap comparableInstance(Bitmap other);

  /**
   * Creates an instance-compatible copy of this bitmap. This bitmap and the returned
   * bitmap are both instances of the same class.
   *
   * @return copy of the this bitmap
   */
  public Bitmap copy();

  /**
   * Creates an instance-compatible copy of highest numBits bits of this bitmap. This
   * bitmap and the returned bitmap are both instances of the same class.
   *
   * @return copy of the this bitmap
   */
  public Bitmap copyHighest(int numBits);

  /**
   * Creates an instance-compatible copy of lowest numBits bits of this bitmap. This
   * bitmap and the returned bitmap are both instances of the same class.
   *
   * @return copy of the this bitmap
   */
  public Bitmap copyLowest(int numBits);

  /**
   * Decrements the integer value of this bitmap by one: bitmap--.
   *
   */
  public void dec();

  public void defBit(int bitIndex, boolean value);

  public void flipBit(int bitIndex);

  public int hsb();

  public int hub();

  /**
   * Increments the integer value of this bitmap by one: bitmap++.
   *
   */
  public void inc();

  public int lsb();

  public int lub();

  /**
   * Updates this bitmap with the result of the bitwise AND with the other bitmap: bitmap
   * = ~(bitmap & other). Assumes the bit length of the other bitmap is AT LEAST as large
   * as the bit length of this bitmap.
   *
   */
  public void nand(final Bitmap other);

  /**
   * Inverts the bits of this bitmap: bitmap = ~bitmap.
   *
   */
  public void not();

  /**
   * Updates this bitmap with the result of the bitwise AND with the other bitmap: bitmap
   * = (bitmap | other). Assumes the bit length of the other bitmap is AT LEAST as large
   * as the bit length of this bitmap.
   *
   */
  public void or(final Bitmap other);

  public void setBit(final int bitIndex);

  public void setBits(final boolean value);

  public boolean testBit(final int bitIndex);

  /**
   * Updates this bitmap with the result of the bitwise AND with the other bitmap: bitmap
   * = (bitmap ^ other). Assumes the bit length of the other bitmap is AT LEAST as large
   * as the bit length of this bitmap.
   *
   */
  public void xor(final Bitmap other);

  /**
   * Returns a string for the integer represented by this bitmap
   */
  public String toNumberString();
}