/*
   Copyright (C) 1997,1998,1999
   Kenji Hiranabe, Eiwa System Management, Inc.

   This program is free software.
   Implemented by Kenji Hiranabe(hiranabe@esm.co.jp),
   conforming to the Java(TM) 3D API specification by Sun Microsystems.

   Permission to use, copy, modify, distribute and sell this software
   and its documentation for any purpose is hereby granted without fee,
   provided that the above copyright notice appear in all copies and
   that both that copyright notice and this permission notice appear
   in supporting documentation. Kenji Hiranabe and Eiwa System Management,Inc.
   makes no representations about the suitability of this software for any
   purpose.  It is provided "AS IS" with NO WARRANTY.
*/
package org.jmol.util;

import java.io.Serializable;

/**
 * A single precision floating point 3 by 3 matrix.
 * 
 * @version specification 1.1, implementation $Revision$, $Date:
 *          2006/07/28 17:01:33 $
 * @author Kenji hiranabe
 * 
 *         additions by Bob Hanson hansonr@stolaf.edu 9/30/2012 for unique
 *         constructor and method names for the optimization of compiled
 *         JavaScript using Java2Script
 */
public class Matrix3f implements Serializable {

  /**
   * The first element of the first row.
   */
  public float m00;

  /**
   * The second element of the first row.
   */
  public float m01;

  /**
   * third element of the first row.
   */
  public float m02;

  /**
   * The first element of the second row.
   */
  public float m10;

  /**
   * The second element of the second row.
   */
  public float m11;

  /**
   * The third element of the second row.
   */
  public float m12;

  /**
   * The first element of the third row.
   */
  public float m20;

  /**
   * The second element of the third row.
   */
  public float m21;

  /**
   * The third element of the third row.
   */
  public float m22;

  /**
   * Constructs and initializes a Matrix3f to all zeros.
   */
  public Matrix3f() {
  }

  /**
   * Constructs and initializes a Matrix3f from the specified 9 element array.
   * this.m00 =v[0], this.m01=v[1], etc.
   * 
   * @param v
   *        the array of length 9 containing in order
   * @return m
   */
  public static Matrix3f newA(float[] v) {
    Matrix3f m = new Matrix3f();
    m.setA(v);
    return m;
  }

  /**
   * Constructs a new matrix with the same values as the Matrix3f parameter.
   * 
   * @param m1
   *        The source matrix.
   * @return m
   */
  public static Matrix3f newM(Matrix3f m1) {
    Matrix3f m = new Matrix3f();
    m.m00 = m1.m00;
    m.m01 = m1.m01;
    m.m02 = m1.m02;

    m.m10 = m1.m10;
    m.m11 = m1.m11;
    m.m12 = m1.m12;

    m.m20 = m1.m20;
    m.m21 = m1.m21;
    m.m22 = m1.m22;
    return m;
  }

  /**
   * Returns a string that contains the values of this Matrix3f.
   * 
   * @return the String representation
   */
  @Override
  public String toString() {
    return "[\n  [" + m00 + "\t" + m01 + "\t" + m02 + "]" + "\n  [" + m10
        + "\t" + m11 + "\t" + m12 + "]" + "\n  [" + m20 + "\t" + m21 + "\t"
        + m22 + "] ]";
  }

  /**
   * Sets this Matrix3f to identity.
   */
  public final void setIdentity() {
    m00 = 1.0f;
    m01 = 0.0f;
    m02 = 0.0f;
    m10 = 0.0f;
    m11 = 1.0f;
    m12 = 0.0f;
    m20 = 0.0f;
    m21 = 0.0f;
    m22 = 1.0f;
  }

  /**
   * Sets the specified element of this matrix3d to the value provided.
   * 
   * @param row
   *        the row number to be modified (zero indexed)
   * @param column
   *        the column number to be modified (zero indexed)
   * @param value
   *        the new value
   */
  public final void setElement(int row, int column, float value) {
    if (row == 0)
      if (column == 0)
        m00 = value;
      else if (column == 1)
        m01 = value;
      else if (column == 2)
        m02 = value;
      else
        throw new ArrayIndexOutOfBoundsException(
            "column must be 0 to 2 and is " + column);
    else if (row == 1)
      if (column == 0)
        m10 = value;
      else if (column == 1)
        m11 = value;
      else if (column == 2)
        m12 = value;
      else
        throw new ArrayIndexOutOfBoundsException(
            "column must be 0 to 2 and is " + column);
    else if (row == 2)
      if (column == 0)
        m20 = value;
      else if (column == 1)
        m21 = value;
      else if (column == 2)
        m22 = value;
      else
        throw new ArrayIndexOutOfBoundsException(
            "column must be 0 to 2 and is " + column);
    else
      throw new ArrayIndexOutOfBoundsException("row must be 0 to 2 and is "
          + row);
  }

  /**
   * Retrieves the value at the specified row and column of this matrix.
   * 
   * @param row
   *        the row number to be retrieved (zero indexed)
   * @param column
   *        the column number to be retrieved (zero indexed)
   * @return the value at the indexed element
   */
  public final float getElement(int row, int column) {
    if (row == 0)
      if (column == 0)
        return m00;
      else if (column == 1)
        return m01;
      else if (column == 2)
        return m02;
      else
        throw new ArrayIndexOutOfBoundsException(
            "column must be 0 to 2 and is " + column);
    else if (row == 1)
      if (column == 0)
        return m10;
      else if (column == 1)
        return m11;
      else if (column == 2)
        return m12;
      else
        throw new ArrayIndexOutOfBoundsException(
            "column must be 0 to 2 and is " + column);

    else if (row == 2)
      if (column == 0)
        return m20;
      else if (column == 1)
        return m21;
      else if (column == 2)
        return m22;
      else
        throw new ArrayIndexOutOfBoundsException(
            "column must be 0 to 2 and is " + column);
    else
      throw new ArrayIndexOutOfBoundsException("row must be 0 to 2 and is "
          + row);
  }

  /**
   * Sets the specified row of this matrix3d to the three values provided.
   * 
   * @param row
   *        the row number to be modified (zero indexed)
   * @param x
   *        the first column element
   * @param y
   *        the second column element
   * @param z
   *        the third column element
   */
  public final void setRow(int row, float x, float y, float z) {
    if (row == 0) {
      m00 = x;
      m01 = y;
      m02 = z;
    } else if (row == 1) {
      m10 = x;
      m11 = y;
      m12 = z;
    } else if (row == 2) {
      m20 = x;
      m21 = y;
      m22 = z;
    } else {
      throw new ArrayIndexOutOfBoundsException("row must be 0 to 2 and is "
          + row);
    }
  }

  /**
   * Sets the specified row of this matrix3d to the Vector provided.
   * 
   * @param row
   *        the row number to be modified (zero indexed)
   * @param v
   *        the replacement row
   */
  public final void setRowV(int row, Vector3f v) {
    if (row == 0) {
      m00 = v.x;
      m01 = v.y;
      m02 = v.z;
    } else if (row == 1) {
      m10 = v.x;
      m11 = v.y;
      m12 = v.z;
    } else if (row == 2) {
      m20 = v.x;
      m21 = v.y;
      m22 = v.z;
    } else {
      throw new ArrayIndexOutOfBoundsException("row must be 0 to 2 and is "
          + row);
    }
  }

  /**
   * Copies the matrix values in the specified row into the array parameter.
   * 
   * @param row
   *        the matrix row
   * @param v
   *        The array into which the matrix row values will be copied
   */
  public final void getRow(int row, float v[]) {
    if (row == 0) {
      v[0] = m00;
      v[1] = m01;
      v[2] = m02;
    } else if (row == 1) {
      v[0] = m10;
      v[1] = m11;
      v[2] = m12;
    } else if (row == 2) {
      v[0] = m20;
      v[1] = m21;
      v[2] = m22;
    } else {
      throw new ArrayIndexOutOfBoundsException("row must be 0 to 2 and is "
          + row);
    }
  }

  /**
   * Sets the specified row of this matrix3d to the four values provided.
   * 
   * @param row
   *        the row number to be modified (zero indexed)
   * @param v
   *        the replacement row
   */
  public final void setRowA(int row, float v[]) {
    if (row == 0) {
      m00 = v[0];
      m01 = v[1];
      m02 = v[2];
    } else if (row == 1) {
      m10 = v[0];
      m11 = v[1];
      m12 = v[2];
    } else if (row == 2) {
      m20 = v[0];
      m21 = v[1];
      m22 = v[2];
    } else {
      throw new ArrayIndexOutOfBoundsException("row must be 0 to 2 and is "
          + row);
    }
  }

  /**
   * Sets the specified column of this matrix3d to the three values provided.
   * 
   * @param column
   *        the column number to be modified (zero indexed)
   * @param x
   *        the first row element
   * @param y
   *        the second row element
   * @param z
   *        the third row element
   */
  public final void setColumn(int column, float x, float y, float z) {
    if (column == 0) {
      m00 = x;
      m10 = y;
      m20 = z;
    } else if (column == 1) {
      m01 = x;
      m11 = y;
      m21 = z;
    } else if (column == 2) {
      m02 = x;
      m12 = y;
      m22 = z;
    } else {
      throw new ArrayIndexOutOfBoundsException("column must be 0 to 2 and is "
          + column);
    }
  }

  /**
   * Sets the specified column of this matrix3d to the vector provided.
   * 
   * @param column
   *        the column number to be modified (zero indexed)
   * @param v
   *        the replacement column
   */
  public final void setColumnV(int column, Vector3f v) {
    if (column == 0) {
      m00 = v.x;
      m10 = v.y;
      m20 = v.z;
    } else if (column == 1) {
      m01 = v.x;
      m11 = v.y;
      m21 = v.z;
    } else if (column == 2) {
      m02 = v.x;
      m12 = v.y;
      m22 = v.z;
    } else {
      throw new ArrayIndexOutOfBoundsException("column must be 0 to 2 and is "
          + column);
    }
  }

  /**
   * Sets the specified column of this matrix3d to the four values provided.
   * 
   * @param column
   *        the column number to be modified (zero indexed)
   * @param v
   *        the replacement column
   */
  public final void setColumnA(int column, float v[]) {
    if (column == 0) {
      m00 = v[0];
      m10 = v[1];
      m20 = v[2];
    } else if (column == 1) {
      m01 = v[0];
      m11 = v[1];
      m21 = v[2];
    } else if (column == 2) {
      m02 = v[0];
      m12 = v[1];
      m22 = v[2];
    } else {
      throw new ArrayIndexOutOfBoundsException("column must be 0 to 2 and is "
          + column);
    }
  }

  /**
   * Copies the matrix values in the specified column into the vector parameter.
   * 
   * @param column
   *        the matrix column
   * @param v
   *        The vector into which the matrix row values will be copied
   */
  public final void getColumnV(int column, Vector3f v) {
    if (column == 0) {
      v.x = m00;
      v.y = m10;
      v.z = m20;
    } else if (column == 1) {
      v.x = m01;
      v.y = m11;
      v.z = m21;
    } else if (column == 2) {
      v.x = m02;
      v.y = m12;
      v.z = m22;
    } else {
      throw new ArrayIndexOutOfBoundsException("column must be 0 to 2 and is "
          + column);
    }
  }

  /**
   * Copies the matrix values in the specified column into the array parameter.
   * 
   * @param column
   *        the matrix column
   * @param v
   *        The array into which the matrix row values will be copied
   */
  public final void getColumn(int column, float v[]) {
    if (column == 0) {
      v[0] = m00;
      v[1] = m10;
      v[2] = m20;
    } else if (column == 1) {
      v[0] = m01;
      v[1] = m11;
      v[2] = m21;
    } else if (column == 2) {
      v[0] = m02;
      v[1] = m12;
      v[2] = m22;
    } else {
      throw new ArrayIndexOutOfBoundsException("column must be 0 to 2 and is "
          + column);
    }
  }

  /**
   * Sets the value of this matrix to sum of itself and matrix m1.
   * 
   * @param m1
   *        the other matrix
   */
  public final void add(Matrix3f m1) {
    m00 += m1.m00;
    m01 += m1.m01;
    m02 += m1.m02;
    m10 += m1.m10;
    m11 += m1.m11;
    m12 += m1.m12;
    m20 += m1.m20;
    m21 += m1.m21;
    m22 += m1.m22;
  }

  /**
   * Sets the value of this matrix to the matrix difference of itself and matrix
   * m1 (this = this - m1).
   * 
   * @param m1
   *        the other matrix
   */
  public final void sub(Matrix3f m1) {
    m00 -= m1.m00;
    m01 -= m1.m01;
    m02 -= m1.m02;
    m10 -= m1.m10;
    m11 -= m1.m11;
    m12 -= m1.m12;
    m20 -= m1.m20;
    m21 -= m1.m21;
    m22 -= m1.m22;
  }

  /**
   * Sets the value of this matrix to its transpose.
   */
  public final void transpose() {
    float tmp = m01;
    m01 = m10;
    m10 = tmp;

    tmp = m02;
    m02 = m20;
    m20 = tmp;

    tmp = m12;
    m12 = m21;
    m21 = tmp;

  }

  /**
   * Sets the value of this matrix to the transpose of the argument matrix
   * 
   * @param m1
   *        the matrix to be transposed
   */
  public final void transposeM(Matrix3f m1) {
    // alias-safe
    setM(m1);
    transpose();
  }

  /**
   * Sets the value of this matrix to the matrix conversion of the single
   * precision axis and angle argument.
   * 
   * @param a1
   *        the axis and angle to be converted
   */
  public final void setAA(AxisAngle4f a1) {
    setFromAxisAngle(a1.x, a1.y, a1.z, a1.angle);
  }

  private void setFromAxisAngle(double x, double y, double z, double angle) {
    // Taken from Rick's which is taken from Wertz. pg. 412
    // Bug Fixed and changed into right-handed by hiranabe
    double n = Math.sqrt(x * x + y * y + z * z);
    // zero-div may occur
    n = 1 / n;
    x *= n;
    y *= n;
    z *= n;
    double c = Math.cos(angle);
    double s = Math.sin(angle);
    double omc = 1.0 - c;
    m00 = (float) (c + x * x * omc);
    m11 = (float) (c + y * y * omc);
    m22 = (float) (c + z * z * omc);

    double tmp1 = x * y * omc;
    double tmp2 = z * s;
    m01 = (float) (tmp1 - tmp2);
    m10 = (float) (tmp1 + tmp2);

    tmp1 = x * z * omc;
    tmp2 = y * s;
    m02 = (float) (tmp1 + tmp2);
    m20 = (float) (tmp1 - tmp2);

    tmp1 = y * z * omc;
    tmp2 = x * s;
    m12 = (float) (tmp1 - tmp2);
    m21 = (float) (tmp1 + tmp2);
  }

  /**
   * Sets the value of this matrix to the double value of the Matrix3f argument.
   * 
   * @param m1
   *        the matrix3f
   */
  public final void setM(Matrix3f m1) {
    m00 = m1.m00;
    m01 = m1.m01;
    m02 = m1.m02;
    m10 = m1.m10;
    m11 = m1.m11;
    m12 = m1.m12;
    m20 = m1.m20;
    m21 = m1.m21;
    m22 = m1.m22;
  }

  /**
   * Sets the values in this Matrix3f equal to the row-major array parameter
   * (ie, the first four elements of the array will be copied into the first row
   * of this matrix, etc.).
   * 
   * @param m
   */
  public final void setA(float m[]) {
    m00 = m[0];
    m01 = m[1];
    m02 = m[2];
    m10 = m[3];
    m11 = m[4];
    m12 = m[5];
    m20 = m[6];
    m21 = m[7];
    m22 = m[8];
  }

  /**
   * Sets the value of this matrix to the matrix inverse of the passed matrix
   * m1.
   * 
   * @param m1
   *        the matrix to be inverted
   */
  public final void invertM(Matrix3f m1) {
    setM(m1);
    invert();
  }

  /**
   * Sets the value of this matrix to its inverse.
   */
  public final void invert() {
    double s = determinant();
    if (s == 0.0)
      return;
    s = 1 / s;
    // alias-safe way.
    set(m11 * m22 - m12 * m21, m02 * m21 - m01 * m22, m01 * m12 - m02 * m11,
        m12 * m20 - m10 * m22, m00 * m22 - m02 * m20, m02 * m10 - m00 * m12,
        m10 * m21 - m11 * m20, m01 * m20 - m00 * m21, m00 * m11 - m01 * m10);
    mulf((float) s);
  }

  /**
   * Computes the determinant of this matrix.
   * 
   * @return the determinant of the matrix
   */
  public final float determinant() {
    // less *,+,- calculation than expanded expression.
    return m00 * (m11 * m22 - m21 * m12) - m01 * (m10 * m22 - m20 * m12) + m02
        * (m10 * m21 - m20 * m11);
  }

  /**
   * Sets the value of this matrix to a scale matrix with the passed scale
   * amount.
   * 
   * @param scale
   *        the scale factor for the matrix
   */
  public final void setScale(float scale) {
    m00 = scale;
    m01 = 0.0f;
    m02 = 0.0f;
    m10 = 0.0f;
    m11 = scale;
    m12 = 0.0f;
    m20 = 0.0f;
    m21 = 0.0f;
    m22 = scale;
  }

  /**
   * Sets the value of this matrix to a rotation matrix about the x axis by the
   * passed angle.
   * 
   * @param angle
   *        the angle to rotate about the X axis in radians
   */
  public final void rotX(float angle) {
    double c = Math.cos(angle);
    double s = Math.sin(angle);
    m00 = 1.0f;
    m01 = 0.0f;
    m02 = 0.0f;
    m10 = 0.0f;
    m11 = (float) c;
    m12 = (float) -s;
    m20 = 0.0f;
    m21 = (float) s;
    m22 = (float) c;
  }

  /**
   * Sets the value of this matrix to a rotation matrix about the y axis by the
   * passed angle.
   * 
   * @param angle
   *        the angle to rotate about the Y axis in radians
   */
  public final void rotY(float angle) {
    double c = Math.cos(angle);
    double s = Math.sin(angle);
    m00 = (float) c;
    m01 = 0.0f;
    m02 = (float) s;
    m10 = 0.0f;
    m11 = 1.0f;
    m12 = 0.0f;
    m20 = (float) -s;
    m21 = 0.0f;
    m22 = (float) c;
  }

  /**
   * Sets the value of this matrix to a rotation matrix about the z axis by the
   * passed angle.
   * 
   * @param angle
   *        the angle to rotate about the Z axis in radians
   */
  public final void rotZ(float angle) {
    double c = Math.cos(angle);
    double s = Math.sin(angle);
    m00 = (float) c;
    m01 = (float) -s;
    m02 = 0.0f;
    m10 = (float) s;
    m11 = (float) c;
    m12 = 0.0f;
    m20 = 0.0f;
    m21 = 0.0f;
    m22 = 1.0f;
  }

  /**
   * Multiplies each element of this matrix by a scalar.
   * 
   * @param scalar
   *        The scalar multiplier.
   */
  public final void mulf(float scalar) {
    m00 *= scalar;
    m01 *= scalar;
    m02 *= scalar;
    m10 *= scalar;
    m11 *= scalar;
    m12 *= scalar;
    m20 *= scalar;
    m21 *= scalar;
    m22 *= scalar;
  }

  /**
   * Sets the value of this matrix to the result of multiplying itself with
   * matrix m1.
   * 
   * @param m1
   *        the other matrix
   */
  public final void mul(Matrix3f m1) {
    mul2(this, m1);
  }

  /**
   * Sets the value of this matrix to the result of multiplying the two argument
   * matrices together.
   * 
   * @param m1
   *        the first matrix
   * @param m2
   *        the second matrix
   */
  public final void mul2(Matrix3f m1, Matrix3f m2) {
    // alias-safe way.
    set(m1.m00 * m2.m00 + m1.m01 * m2.m10 + m1.m02 * m2.m20, m1.m00 * m2.m01
        + m1.m01 * m2.m11 + m1.m02 * m2.m21, m1.m00 * m2.m02 + m1.m01 * m2.m12
        + m1.m02 * m2.m22,

    m1.m10 * m2.m00 + m1.m11 * m2.m10 + m1.m12 * m2.m20, m1.m10 * m2.m01
        + m1.m11 * m2.m11 + m1.m12 * m2.m21, m1.m10 * m2.m02 + m1.m11 * m2.m12
        + m1.m12 * m2.m22,

    m1.m20 * m2.m00 + m1.m21 * m2.m10 + m1.m22 * m2.m20, m1.m20 * m2.m01
        + m1.m21 * m2.m11 + m1.m22 * m2.m21, m1.m20 * m2.m02 + m1.m21 * m2.m12
        + m1.m22 * m2.m22);
  }

  /**
   * Returns true if the Object o is of type Matrix3f and all of the data
   * members of t1 are equal to the corresponding data members in this Matrix3f.
   * 
   * @param o
   *        the object with which the comparison is made.
   */
  @Override
  public boolean equals(Object o) {
    if (!(o instanceof Matrix3f))
      return false;
    Matrix3f m = (Matrix3f) o;
    return m00 == m.m00 && m01 == m.m01 && m02 == m.m02 && m10 == m.m10
        && m11 == m.m11 && m12 == m.m12 && m20 == m.m20 && m21 == m.m21
        && m22 == m.m22;
  }

  /**
   * Returns a hash number based on the data values in this object. Two
   * different Matrix3f objects with identical data values (ie, returns true for
   * equals(Matrix3f) ) will return the same hash number. Two objects with
   * different data members may return the same hash value, although this is not
   * likely.
   * 
   * @return the integer hash value
   */
  @Override
  public int hashCode() {
    return Tuple3f.floatToIntBits0(m00) ^ Tuple3f.floatToIntBits0(m01)
        ^ Tuple3f.floatToIntBits0(m02) ^ Tuple3f.floatToIntBits0(m10)
        ^ Tuple3f.floatToIntBits0(m11) ^ Tuple3f.floatToIntBits0(m12)
        ^ Tuple3f.floatToIntBits0(m20) ^ Tuple3f.floatToIntBits0(m21)
        ^ Tuple3f.floatToIntBits0(m22);
  }

  /**
   * Sets this matrix to all zeros.
   */
  public final void setZero() {
    m00 = 0.0f;
    m01 = 0.0f;
    m02 = 0.0f;
    m10 = 0.0f;
    m11 = 0.0f;
    m12 = 0.0f;
    m20 = 0.0f;
    m21 = 0.0f;
    m22 = 0.0f;
  }

  public final void transform(Tuple3f t) {
    // alias-safe
    transform2(t, t);
  }

  /**
   * Transform the vector vec using this Matrix3f and place the result into
   * vecOut.
   * 
   * @param t
   *        the single precision vector to be transformed
   * @param result
   *        the vector into which the transformed values are placed
   */
  public final void transform2(Tuple3f t, Tuple3f result) {
    // alias-safe
    result.set(m00 * t.x + m01 * t.y + m02 * t.z, m10 * t.x + m11 * t.y + m12
        * t.z, m20 * t.x + m21 * t.y + m22 * t.z);
  }

  /**
   * Sets 9 values
   * 
   * @param m00
   * @param m01
   * @param m02
   * @param m10
   * @param m11
   * @param m12
   * @param m20
   * @param m21
   * @param m22
   */
  private void set(float m00, float m01, float m02, float m10, float m11,
                   float m12, float m20, float m21, float m22) {
    this.m00 = m00;
    this.m01 = m01;
    this.m02 = m02;
    this.m10 = m10;
    this.m11 = m11;
    this.m12 = m12;
    this.m20 = m20;
    this.m21 = m21;
    this.m22 = m22;
  }

}
