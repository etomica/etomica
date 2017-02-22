/* $RCSfile$
 * $Author$
 * $Date$
 * $Revision$
 *
 * Copyright (C) 2003-2005  The Jmol Development Team
 *
 * Contact: jmol-developers@lists.sf.net
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.jmol.util;


import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;


final public class ArrayUtil {

  /**
   * Very important that this not be used with Int32Array or Float32Array,
   * because it is not initialized to all zeros in MSIE 9.
   * 
   * @param array
   * @param minimumLength
   * @return array
   */
  public static Object ensureLength(Object array, int minimumLength) {
    if (array != null && getLength(array) >= minimumLength)
      return array;
    return arrayCopyObject(array, minimumLength);
  }

  public static String[] ensureLengthS(String[] array, int minimumLength) {
    if (array != null && array.length >= minimumLength)
      return array;
    return arrayCopyS(array, minimumLength);
  }

  public static float[] ensureLengthA(float[] array, int minimumLength) {
    if (array != null && array.length >= minimumLength)
      return array;
    return arrayCopyF(array, minimumLength);
  }

  public static int[] ensureLengthI(int[] array, int minimumLength) {
    if (array != null && array.length >= minimumLength)
      return array;
    return arrayCopyI(array, minimumLength);
  }

  public static short[] ensureLengthShort(short[] array, int minimumLength) {
    if (array != null && array.length >= minimumLength)
      return array;
    return arrayCopyShort(array, minimumLength);
  }

  public static byte[] ensureLengthByte(byte[] array, int minimumLength) {
    if (array != null && array.length >= minimumLength)
      return array;
    return arrayCopyByte(array, minimumLength);
  }

  /**
   * Very important that this not be used with Int32Array or Float32Array,
   * because it is not initialized to all zeros in MSIE 9.
   * 
   * @param array
   * @return array
   */
  public static Object doubleLength(Object array) {
    return arrayCopyObject(array, (array == null ? 16 : 2 * getLength(array)));
  }

  public static String[] doubleLengthS(String[] array) {
    return arrayCopyS(array, (array == null ? 16 : 2 * array.length));
  }

  public static float[] doubleLengthF(float[] array) {
    return arrayCopyF(array, (array == null ? 16 : 2 * array.length));
  }

  public static int[] doubleLengthI(int[] array) {
    return arrayCopyI(array, (array == null ? 16 : 2 * array.length));
  }

  public static short[] doubleLengthShort(short[] array) {
    return arrayCopyShort(array, (array == null ? 16 : 2 * array.length));
  }

  public static byte[] doubleLengthByte(byte[] array) {
    return arrayCopyByte(array, (array == null ? 16 : 2 * array.length));
  }

  public static boolean[] doubleLengthBool(boolean[] array) {
    return arrayCopyBool(array, (array == null ? 16 : 2 * array.length));
  }

  public static Object deleteElements(Object array, int firstElement,
                                     int nElements) {
    if (nElements == 0 || array == null)
      return array;
    int oldLength = getLength(array);
    if (firstElement >= oldLength)
      return array;
    int n = oldLength - (firstElement + nElements);
    if (n < 0)
      n = 0;
    Object t = newInstanceO(array, firstElement + n);
    if (firstElement > 0)
      System.arraycopy(array, 0, t, 0, firstElement);
    if (n > 0)
      System.arraycopy(array, firstElement + nElements, t, firstElement, n);
    return t;
  }

  /**
   * note -- cannot copy if array is null!
   * 
   * @param array
   * @param newLength
   * @return array
   */
  public static Object arrayCopyObject(Object array, int newLength) {
    //System.out.println("ArrayUtil.copy " + newLength + " " + array + "  ");
    if (array == null) {
      return null; // We can't allocate since we don't know the type of array
    }
    int oldLength = getLength(array);
    if (newLength == oldLength)
      return array;
    Object t = newInstanceO(array, newLength);
    System.arraycopy(array, 0, t, 0, oldLength < newLength ? oldLength
        : newLength);
    return t;

  }

  /**
   * Very important that this not be used with Int32Array or Float32Array,
   * because those need to be initialized to all zeros in MSIE 9, and
   * MSIE 9 cannot distinguish Int32Array or Float32Array from Array.
   * 
   * @param array
   * @param n
   * @return array
   */
  private static Object newInstanceO(Object array, int n) {
    /**
     * @j2sNative
     * 
     * return new Array(n);
     * 
     */
    {
      return Array.newInstance(array.getClass().getComponentType(), n);
    }
  }

  private static int getLength(Object array) {
    /**
     * @j2sNative
     * 
     *  return array.length
     *   
     */
    {
      return Array.getLength(array);
    }
  }

  public static String[] arrayCopyS(String[] array, int newLength) {
    if (newLength < 0)
      newLength = array.length;
    String[] t = new String[newLength];
    if (array != null) {
      int oldLength = array.length;
      System.arraycopy(array, 0, t, 0, oldLength < newLength ? oldLength
          : newLength);
    }
    return t;
  }

  public static int[][] arrayCopyII(int[][] array, int newLength) {
    int[][] t = newInt2(newLength);
    if (array != null) {
      int oldLength = array.length;
      System.arraycopy(array, 0, t, 0, oldLength < newLength ? oldLength
          : newLength);
    }
    return t;
  }

  public static Point3f[] arrayCopyPt(Point3f[] array, int newLength) {
    if (newLength < 0)
      newLength = array.length;
    Point3f[] t = new Point3f[newLength];
    if (array != null) {
      int oldLength = array.length;
      System.arraycopy(array, 0, t, 0, oldLength < newLength ? oldLength
          : newLength);
    }
    return t;
  }

  public static float[] arrayCopyF(float[] array, int newLength) {
    if (newLength < 0)
      newLength = array.length;
    float[] t = new float[newLength];
    if (array != null) {
      int oldLength = array.length;
      System.arraycopy(array, 0, t, 0, oldLength < newLength ? oldLength
          : newLength);
    }
    return t;
  }

  public static int[] arrayCopyI(int[] array, int newLength) {
    if (newLength < 0)
      newLength = array.length;
    int[] t = new int[newLength];
    if (array != null) {
      int oldLength = array.length;
      System.arraycopy(array, 0, t, 0, oldLength < newLength ? oldLength
          : newLength);
    }
    return t;
  }

  /**
   * a specialized method that allows copying from a starting point either
   * to the end or to the middle (color schemes, especially)
   * @param array
   * @param i0
   * @param n
   * @return array or null
   */
  public static int[] arrayCopyRangeI(int[] array, int i0, int n) {
    if (array == null)
      return null;
    int oldLength = array.length;
    if (n == -1) n = oldLength;
    if (n == -2) n = oldLength / 2;
    n = n - i0;
    int[] t = new int[n];
    System.arraycopy(array, i0, t, 0, n);
    return t;
  }

  public static int[] arrayCopyRangeRevI(int[] array, int i0, int n) {
    if (array == null)
      return null;
    int[] t = arrayCopyRangeI(array, i0, n);
    if (n < 0)
      n = array.length;
    for (int i = n / 2; --i >= 0;)
      swapInt(t, i, n - 1 - i);
    return t;
  }

  public static short[] arrayCopyShort(short[] array, int newLength) {
    if (newLength < 0)
      newLength = array.length;
    short[] t = new short[newLength];
    if (array != null) {
      int oldLength = array.length;
      System.arraycopy(array, 0, t, 0, oldLength < newLength ? oldLength
          : newLength);
    }
    return t;
  }

  public static byte[] arrayCopyByte(byte[] array, int newLength) {
    if (newLength < 0)
      newLength = array.length;
    byte[] t = new byte[newLength];
    if (array != null) {
      int oldLength = array.length;
      System.arraycopy(array, 0, t, 0, oldLength < newLength ? oldLength
          : newLength);
    }
    return t;
  }

  public static boolean[] arrayCopyBool(boolean[] array, int newLength) {
    if (newLength < 0)
      newLength = array.length;
    boolean[] t = new boolean[newLength];
    if (array != null) {
      int oldLength = array.length;
      System.arraycopy(array, 0, t, 0, oldLength < newLength ? oldLength
          : newLength);
    }
    return t;
  }

  public static void swapInt(int[] array, int indexA, int indexB) {
    int t = array[indexA];
    array[indexA] = array[indexB];
    array[indexB] = t;
  }

  /*
  public static void swap(short[] array, int indexA, int indexB) {
    short t = array[indexA];
    array[indexA] = array[indexB];
    array[indexB] = t;
  }

  public static void swap(float[] array, int indexA, int indexB) {
    float t = array[indexA];
    array[indexA] = array[indexB];
    array[indexB] = t;
  }
  */
  
  public static String dumpArray(String msg, float[][] A, int x1, int x2, int y1, int y2) {
    String s = "dumpArray: " + msg + "\n";
    for (int x = x1; x <= x2; x++)
      s += "\t*" + x + "*";
    for (int y = y2; y >= y1; y--) {
      s += "\n*" + y + "*";
      for (int x = x1; x <= x2; x++)
        s += "\t" + (x < A.length && y < A[x].length ? A[x][y] : Float.NaN);
    }
    return s;
  }

  public static String dumpIntArray(int[] A, int n) {
    String str = "";
    for (int i = 0; i < n; i++)
      str += " " + A[i];
    return str;
  }

  public static String sortedItem(List<String> v, int n) {
    if (v.size() == 0)
      return null;
    if (v.size() == 1)
      return v.get(0);
    String[] keys = v.toArray(new String[v.size()]);
    Arrays.sort(keys);
    return keys[n % keys.length];
  }

  /**
   * Helper method for creating a List<T>[] without warnings.
   * 
   * @param <T> Type of objects in the list.
   * @param size Array size.
   * @return Array of List<T>
   */
  @SuppressWarnings("unchecked")
  public static <T> List<T>[] createArrayOfArrayList(int size) {
    return new ArrayList[size];
  }

  /**
   * Helper method for creating a Map<K, V>[] without warnings.
   * 
   * @param <K> Type of object for the keys in the map.
   * @param <V> Type of object for the values in the map.
   * @param size Array size.
   * @return Array of Map<K, V>
   */
  @SuppressWarnings("unchecked")
  public static <K, V> Map<K, V>[] createArrayOfHashtable(int size) {
    return new Hashtable[size];
  }

  public static void swap(Object[] o, int i, int j) {
    Object oi = o[i];
    o[i] = o[j];
    o[j] = oi;
  }

  public static float[][] newFloat2(int n) {
    /**
     * @j2sNative
     * 
     * return Clazz.newArray(n, null);
     * 
     */
    {
    return new float[n][];
    }
  }

  public static int[][] newInt2(int n) {
    /**
     * @j2sNative
     * 
     * return Clazz.newArray(n, null);
     * 
     */
    {
    return new int[n][];
    }
  }

  public static float[][][] newFloat3(int nx, int ny) {
    /**
     * @j2sNative
     * 
     * return Clazz.newArray(nx, null);
     * 
     */
    {
      return (ny < 0 ? new float[nx][][] : new float[nx][ny][]);
    }
  }

  public static int[][][][] newInt4(int n) {
    /**
     * @j2sNative
     * 
     * return Clazz.newArray(n, null);
     * 
     */
    {
    return new int[n][][][];
    }
  }

  public static short[][] newShort2(int n) {
    /**
     * @j2sNative
     * 
     * return Clazz.newArray(n, null);
     * 
     */
    {
    return new short[n][];
    }
  }

  public static byte[][] newByte2(int n) {
    /**
     * @j2sNative
     * 
     * return Clazz.newArray(n, null);
     * 
     */
    {
    return new byte[n][];
    }
  }

}
