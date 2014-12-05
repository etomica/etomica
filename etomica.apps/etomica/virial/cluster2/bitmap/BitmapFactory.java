/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.bitmap;

import etomica.virial.cluster2.bitmap.impl.BitmapOfBigInteger;
import etomica.virial.cluster2.bitmap.impl.BitmapOfLong;
import etomica.virial.cluster2.bitmap.impl.BitmapOfLongVector;

public class BitmapFactory {

  public static boolean useBigInteger = false;
  public static final Bitmap ONE = new BitmapOfLong(1, true);
  public static final Bitmap ZERO = new BitmapOfLong(1, false);
  public static final Bitmap EMPTY = new BitmapOfLong() {

    @Override
    public boolean equals(final Object other) {

      return other == this;
    }

    @Override
    public int compareTo(final Bitmap other) {

      return (this == other) ? 0 : -1;
    }
  };

  public static final Bitmap getBitmap(final String bitmap) {

    if (bitmap.length() <= Bitmap.SZ_LONG) {
      return new BitmapOfLong(bitmap);
    }
    else if (!useBigInteger) {
      return new BitmapOfLongVector(bitmap);
    }
    else {
      return new BitmapOfBigInteger(bitmap);
    }
  }

  public static final Bitmap getBitmap(final int capacity, boolean isSet) {

    // return new BitmapOfLongVector(capacity, isSet);
    if (capacity <= Bitmap.SZ_LONG) {
      return new BitmapOfLong(capacity, isSet);
    }
    else if (!useBigInteger) {
      return new BitmapOfLongVector(capacity, isSet);
    }
    else {
      return new BitmapOfBigInteger(capacity, isSet);
    }
  }
}