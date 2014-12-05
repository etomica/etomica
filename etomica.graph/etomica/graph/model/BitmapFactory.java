/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.model;

import etomica.graph.model.impl.AbstractBitmap;
import etomica.graph.model.impl.BitmapOfLong;
import etomica.graph.model.impl.BitmapOfLongVector;

public class BitmapFactory {

  public static final Bitmap ONE = new BitmapOfLong(1, true);
  public static final Bitmap ZERO = new BitmapOfLong(1, false);
  public static final Bitmap EMPTY = new AbstractBitmap() {

    @Override
    public boolean equals(final Object other) {

      return other == this;
    }

    @Override
    public int compareTo(final Bitmap other) {

      return (this == other) ? 0 : -1;
    }

    @Override
    protected void allocateBitmap() {

      // no-op
    }

    public int bitSize() {

      return 0;
    }

    public void clearBit(int bitIndex) {

    }

    public void setBit(int bitIndex) {

    }

    public boolean testBit(int bitIndex) {

      return false;
    }
  };

  public static final Bitmap createBitmap(final String bitmap) {

    if (bitmap.length() <= Bitmap.SZ_LONG) {
      return new BitmapOfLong(bitmap);
    }
    else {
      return new BitmapOfLongVector(bitmap);
    }
  }

  public static final Bitmap createBitmap(final byte[] bitmap) {

    if (bitmap.length <= Bitmap.SZ_LONG) {
      return new BitmapOfLong(bitmap);
    }
    else {
      return new BitmapOfLongVector(bitmap);
    }
  }

  public static final Bitmap createBitmap(final byte nodeCount, boolean isSet) {

    byte bitSize = (byte) (nodeCount * (nodeCount - 1) / 2);
    if (bitSize <= Bitmap.SZ_LONG) {
      return new BitmapOfLong(bitSize, isSet);
    }
    else {
      return new BitmapOfLongVector(bitSize, isSet);
    }
  }
}