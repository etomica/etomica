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
 * A single precision floating point 4 by 4 matrix.
 * 
 * @version specification 1.1, implementation $Revision$, $Date:
 *          2006/08/01 16:08:49 $
 * @author Kenji hiranabe
 * 
 * additions by Bob Hanson hansonr@stolaf.edu 9/30/2012
 * for unique constructor and method names
 * for the optimization of compiled JavaScript using Java2Script
 */
public class Matrix4f implements Serializable {

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
   * The fourth element of the first row.
   */
  public float m03;

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
   * The fourth element of the second row.
   */
  public float m13;

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
   * The fourth element of the third row.
   */
  public float m23;

  /**
   * The first element of the fourth row.
   */
  public float m30;

  /**
   * The second element of the fourth row.
   */
  public float m31;

  /**
   * The third element of the fourth row.
   */
  public float m32;

  /**
   * The fourth element of the fourth row.
   */
  public float m33;

  public Matrix4f() {
  }

  /**
   * Constructs and initializes a Matrix4f from the specified 16 element array.
   * this.m00 =v[0], this.m01=v[1], etc.
   * 
   * @param v
   *        the array of length 16 containing in order
   * @return m
   */
  public static Matrix4f newA(float[] v) {
    Matrix4f m = new Matrix4f();
    m.m00 = v[0];
    m.m01 = v[1];
    m.m02 = v[2];
    m.m03 = v[3];

    m.m10 = v[4];
    m.m11 = v[5];
    m.m12 = v[6];
    m.m13 = v[7];

    m.m20 = v[8];
    m.m21 = v[9];
    m.m22 = v[10];
    m.m23 = v[11];

    m.m30 = v[12];
    m.m31 = v[13];
    m.m32 = v[14];
    m.m33 = v[15];

    return m;
  }

  /**
   * Constructs a new matrix with the same values as the Matrix4f parameter.
   * 
   * @param m1
   *        the source matrix
   * @return m
   */
  public static Matrix4f newM(Matrix4f m1) {
    Matrix4f m = new Matrix4f();
    m.m00 = m1.m00;
    m.m01 = m1.m01;
    m.m02 = m1.m02;
    m.m03 = m1.m03;

    m.m10 = m1.m10;
    m.m11 = m1.m11;
    m.m12 = m1.m12;
    m.m13 = m1.m13;

    m.m20 = m1.m20;
    m.m21 = m1.m21;
    m.m22 = m1.m22;
    m.m23 = m1.m23;

    m.m30 = m1.m30;
    m.m31 = m1.m31;
    m.m32 = m1.m32;
    m.m33 = m1.m33;
    return m;
  }

  /**
   * Constructs and initializes a Matrix4f from the rotation matrix
   * and translation.
   * @param m1  The rotation matrix representing the rotational components
   * @param t  The translational components of the matrix
   * @return m
   */
 public static Matrix4f newMV(Matrix3f m1, Vector3f t) {
   Matrix4f m = new Matrix4f();
   m.setMV(m1, t);
   return m;
 }

 /**
  * Initializes a Matrix4f from the rotation matrix
  * and translation.
  * @param m1  The rotation matrix representing the rotational components
  * @param t  The translational components of the matrix
  */
  public void setMV(Matrix3f m1, Vector3f t) {
    setM3(m1);
    setTranslation(t);
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
    double n = Math.sqrt(x*x + y*y + z*z);
    // zero-div may occur
    n = 1/n;
    x *= n;
    y *= n;
    z *= n;
    double c = Math.cos(angle);
    double s = Math.sin(angle);
    double omc = 1.0 - c;

    m00 = (float)(c + x*x*omc);
    m11 = (float)(c + y*y*omc);
    m22 = (float)(c + z*z*omc);

    double tmp1 = x*y*omc;
    double tmp2 = z*s;
    m01 = (float)(tmp1 - tmp2);
    m10 = (float)(tmp1 + tmp2);

    tmp1 = x*z*omc;
    tmp2 = y*s;
    m02 = (float)(tmp1 + tmp2);
    m20 = (float)(tmp1 - tmp2);

    tmp1 = y*z*omc;
    tmp2 = x*s;
    m12 = (float)(tmp1 - tmp2);
    m21 = (float)(tmp1 + tmp2);
      }

  /**
   * Sets the value of this matrix to a copy of the passed matrix m1.
   * 
   * @param m1
   *        the matrix to be copied
   */
  public final void setM(Matrix4f m1) {
    m00 = m1.m00;
    m01 = m1.m01;
    m02 = m1.m02;
    m03 = m1.m03;
    m10 = m1.m10;
    m11 = m1.m11;
    m12 = m1.m12;
    m13 = m1.m13;
    m20 = m1.m20;
    m21 = m1.m21;
    m22 = m1.m22;
    m23 = m1.m23;
    m30 = m1.m30;
    m31 = m1.m31;
    m32 = m1.m32;
    m33 = m1.m33;
  }

  /**
   * Returns a string that contains the values of this Matrix4f.
   * 
   * @return the String representation
   */
  @Override
  public String toString() {
    return "[\n  [" + m00 + "\t" + m01 + "\t" + m02 + "\t" + m03 + "]" +
    		"\n  [" + m10 + "\t" + m11 + "\t" + m12 + "\t" + m13 + "]" +
    		"\n  [" + m20 + "\t" + m21 + "\t" + m22 + "\t" + m23 + "]" + 
    		"\n  [" + m30 + "\t" + m31 + "\t" + m32 + "\t" + m33 + "] ]";
  }

  /**
   * Sets this Matrix4f to identity.
   */
  public final void setIdentity() {
    m00 = 1.0f;
    m01 = 0.0f;
    m02 = 0.0f;
    m03 = 0.0f;
    m10 = 0.0f;
    m11 = 1.0f;
    m12 = 0.0f;
    m13 = 0.0f;
    m20 = 0.0f;
    m21 = 0.0f;
    m22 = 1.0f;
    m23 = 0.0f;
    m30 = 0.0f;
    m31 = 0.0f;
    m32 = 0.0f;
    m33 = 1.0f;
  }

  /**
   * Sets the specified element of this matrix4f to the value provided.
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
      else if (column == 3)
        m03 = value;
      else
        throw new ArrayIndexOutOfBoundsException(
            "column must be 0 to 3 and is " + column);
    else if (row == 1)
      if (column == 0)
        m10 = value;
      else if (column == 1)
        m11 = value;
      else if (column == 2)
        m12 = value;
      else if (column == 3)
        m13 = value;
      else
        throw new ArrayIndexOutOfBoundsException(
            "column must be 0 to 3 and is " + column);
    else if (row == 2)
      if (column == 0)
        m20 = value;
      else if (column == 1)
        m21 = value;
      else if (column == 2)
        m22 = value;
      else if (column == 3)
        m23 = value;
      else
        throw new ArrayIndexOutOfBoundsException(
            "column must be 0 to 3 and is " + column);
    else if (row == 3)
      if (column == 0)
        m30 = value;
      else if (column == 1)
        m31 = value;
      else if (column == 2)
        m32 = value;
      else if (column == 3)
        m33 = value;
      else
        throw new ArrayIndexOutOfBoundsException(
            "column must be 0 to 3 and is " + column);
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
      else if (column == 3)
        return m03;
      else
        throw new ArrayIndexOutOfBoundsException(
            "column must be 0 to 3 and is " + column);
    else if (row == 1)
      if (column == 0)
        return m10;
      else if (column == 1)
        return m11;
      else if (column == 2)
        return m12;
      else if (column == 3)
        return m13;
      else
        throw new ArrayIndexOutOfBoundsException(
            "column must be 0 to 3 and is " + column);
    else if (row == 2)
      if (column == 0)
        return m20;
      else if (column == 1)
        return m21;
      else if (column == 2)
        return m22;
      else if (column == 3)
        return m23;
      else
        throw new ArrayIndexOutOfBoundsException(
            "column must be 0 to 3 and is " + column);
    else if (row == 3)
      if (column == 0)
        return m30;
      else if (column == 1)
        return m31;
      else if (column == 2)
        return m32;
      else if (column == 3)
        return m33;
      else
        throw new ArrayIndexOutOfBoundsException(
            "column must be 0 to 3 and is " + column);
    else
      throw new ArrayIndexOutOfBoundsException("row must be 0 to 3 and is "
          + row);
  }

  /**
   * Retrieves the translational components of this matrix.
   * 
   * @param trans
   *        the vector that will receive the translational component
   */
  public final void get(Vector3f trans) {
    trans.x = m03;
    trans.y = m13;
    trans.z = m23;
  }

  /**
   * Gets the upper 3x3 values of this matrix and places them into the matrix
   * m1.
   * 
   * @param m1
   *        The matrix that will hold the values
   */
  public final void getRotationScale(Matrix3f m1) {
    m1.m00 = m00;
    m1.m01 = m01;
    m1.m02 = m02;
    m1.m10 = m10;
    m1.m11 = m11;
    m1.m12 = m12;
    m1.m20 = m20;
    m1.m21 = m21;
    m1.m22 = m22;
  }

  /**
   * Replaces the upper 3x3 matrix values of this matrix with the values in the
   * matrix m1.
   * 
   * @param m1
   *        The matrix that will be the new upper 3x3
   */
  public final void setRotationScale(Matrix3f m1) {
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
   * Sets the specified row of this matrix4f to the four values provided.
   * 
   * @param row
   *        the row number to be modified (zero indexed)
   * @param v
   *        the replacement row
   */
  public final void setRow(int row, float v[]) {
    if (row == 0) {
      m00 = v[0];
      m01 = v[1];
      m02 = v[2];
      m03 = v[3];
    } else if (row == 1) {
      m10 = v[0];
      m11 = v[1];
      m12 = v[2];
      m13 = v[3];
    } else if (row == 2) {
      m20 = v[0];
      m21 = v[1];
      m22 = v[2];
      m23 = v[3];
    } else if (row == 3) {
      m30 = v[0];
      m31 = v[1];
      m32 = v[2];
      m33 = v[3];
    } else {
      throw new ArrayIndexOutOfBoundsException("row must be 0 to 3 and is "
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
      v[3] = m03;
    } else if (row == 1) {
      v[0] = m10;
      v[1] = m11;
      v[2] = m12;
      v[3] = m13;
    } else if (row == 2) {
      v[0] = m20;
      v[1] = m21;
      v[2] = m22;
      v[3] = m23;
    } else if (row == 3) {
      v[0] = m30;
      v[1] = m31;
      v[2] = m32;
      v[3] = m33;
    } else {
      throw new ArrayIndexOutOfBoundsException("row must be 0 to 3 and is "
          + row);
    }
  }

  /**
   * Sets the specified column of this matrix4f to the four values provided.
   * 
   * @param column
   *        the column number to be modified (zero indexed)
   * @param x
   *        the first row element
   * @param y
   *        the second row element
   * @param z
   *        the third row element
   * @param w
   *        the fourth row element
   */
  public final void setColumn4(int column, float x, float y, float z, float w) {
    if (column == 0) {
      m00 = x;
      m10 = y;
      m20 = z;
      m30 = w;
    } else if (column == 1) {
      m01 = x;
      m11 = y;
      m21 = z;
      m31 = w;
    } else if (column == 2) {
      m02 = x;
      m12 = y;
      m22 = z;
      m32 = w;
    } else if (column == 3) {
      m03 = x;
      m13 = y;
      m23 = z;
      m33 = w;
    } else {
      throw new ArrayIndexOutOfBoundsException("column must be 0 to 3 and is "
          + column);
    }
  }

  /**
   * Sets the specified column of this matrix4f to the four values provided.
   * 
   * @param column
   *        the column number to be modified (zero indexed)
   * @param v
   *        the replacement column
   */
  public final void setColumn(int column, float v[]) {
    if (column == 0) {
      m00 = v[0];
      m10 = v[1];
      m20 = v[2];
      m30 = v[3];
    } else if (column == 1) {
      m01 = v[0];
      m11 = v[1];
      m21 = v[2];
      m31 = v[3];
    } else if (column == 2) {
      m02 = v[0];
      m12 = v[1];
      m22 = v[2];
      m32 = v[3];
    } else if (column == 3) {
      m03 = v[0];
      m13 = v[1];
      m23 = v[2];
      m33 = v[3];
    } else {
      throw new ArrayIndexOutOfBoundsException("column must be 0 to 3 and is "
          + column);
    }
  }

  /**
   * Copies the matrix values in the specified column into the array parameter.
   * 
   * @param column
   *        the matrix column
   * @param v
   *        The array into which the matrix column values will be copied
   */
  public final void getColumn(int column, float v[]) {
    if (column == 0) {
      v[0] = m00;
      v[1] = m10;
      v[2] = m20;
      v[3] = m30;
    } else if (column == 1) {
      v[0] = m01;
      v[1] = m11;
      v[2] = m21;
      v[3] = m31;
    } else if (column == 2) {
      v[0] = m02;
      v[1] = m12;
      v[2] = m22;
      v[3] = m32;
    } else if (column == 3) {
      v[0] = m03;
      v[1] = m13;
      v[2] = m23;
      v[3] = m33;
    } else {
      throw new ArrayIndexOutOfBoundsException("column must be 0 to 3 and is "
          + column);
    }
  }

  /**
   * Sets the value of this matrix to the matrix difference of itself and matrix
   * m1 (this = this - m1).
   * 
   * @param m1
   *        the other matrix
   */
  public final void sub(Matrix4f m1) {
    m00 -= m1.m00;
    m01 -= m1.m01;
    m02 -= m1.m02;
    m03 -= m1.m03;
    m10 -= m1.m10;
    m11 -= m1.m11;
    m12 -= m1.m12;
    m13 -= m1.m13;
    m20 -= m1.m20;
    m21 -= m1.m21;
    m22 -= m1.m22;
    m23 -= m1.m23;
    m30 -= m1.m30;
    m31 -= m1.m31;
    m32 -= m1.m32;
    m33 -= m1.m33;
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

    tmp = m03;
    m03 = m30;
    m30 = tmp;

    tmp = m12;
    m12 = m21;
    m21 = tmp;

    tmp = m13;
    m13 = m31;
    m31 = tmp;

    tmp = m23;
    m23 = m32;
    m32 = tmp;
  }


  /**
   * Sets the value of this matrix to the matrix inverse of the passed matrix
   * m1.
   * 
   * @param m1
   *        the matrix to be inverted
   */
  public final void invertM(Matrix4f m1) {
    setM(m1);
    invert();
  }

  /**
   * Sets the value of this matrix to its inverse.
   */
  public final void invert() {
    float s = determinant();
    if (s == 0.0)
      return;
    s = 1 / s;
    // alias-safe way.
    // less *,+,- calculation than expanded expression.
    set(m11 * (m22 * m33 - m23 * m32) + m12 * (m23 * m31 - m21 * m33) + m13
        * (m21 * m32 - m22 * m31), m21 * (m02 * m33 - m03 * m32) + m22
        * (m03 * m31 - m01 * m33) + m23 * (m01 * m32 - m02 * m31), m31
        * (m02 * m13 - m03 * m12) + m32 * (m03 * m11 - m01 * m13) + m33
        * (m01 * m12 - m02 * m11), m01 * (m13 * m22 - m12 * m23) + m02
        * (m11 * m23 - m13 * m21) + m03 * (m12 * m21 - m11 * m22),

    m12 * (m20 * m33 - m23 * m30) + m13 * (m22 * m30 - m20 * m32) + m10
        * (m23 * m32 - m22 * m33), m22 * (m00 * m33 - m03 * m30) + m23
        * (m02 * m30 - m00 * m32) + m20 * (m03 * m32 - m02 * m33), m32
        * (m00 * m13 - m03 * m10) + m33 * (m02 * m10 - m00 * m12) + m30
        * (m03 * m12 - m02 * m13), m02 * (m13 * m20 - m10 * m23) + m03
        * (m10 * m22 - m12 * m20) + m00 * (m12 * m23 - m13 * m22),

    m13 * (m20 * m31 - m21 * m30) + m10 * (m21 * m33 - m23 * m31) + m11
        * (m23 * m30 - m20 * m33), m23 * (m00 * m31 - m01 * m30) + m20
        * (m01 * m33 - m03 * m31) + m21 * (m03 * m30 - m00 * m33), m33
        * (m00 * m11 - m01 * m10) + m30 * (m01 * m13 - m03 * m11) + m31
        * (m03 * m10 - m00 * m13), m03 * (m11 * m20 - m10 * m21) + m00
        * (m13 * m21 - m11 * m23) + m01 * (m10 * m23 - m13 * m20),

    m10 * (m22 * m31 - m21 * m32) + m11 * (m20 * m32 - m22 * m30) + m12
        * (m21 * m30 - m20 * m31), m20 * (m02 * m31 - m01 * m32) + m21
        * (m00 * m32 - m02 * m30) + m22 * (m01 * m30 - m00 * m31), m30
        * (m02 * m11 - m01 * m12) + m31 * (m00 * m12 - m02 * m10) + m32
        * (m01 * m10 - m00 * m11), m00 * (m11 * m22 - m12 * m21) + m01
        * (m12 * m20 - m10 * m22) + m02 * (m10 * m21 - m11 * m20));

    mul(s);
  }

  /**
   * Computes the determinant of this matrix.
   * 
   * @return the determinant of the matrix
   */
  public final float determinant() {
    // less *,+,- calculation than expanded expression.
    return (m00 * m11 - m01 * m10) * (m22 * m33 - m23 * m32)
        - (m00 * m12 - m02 * m10) * (m21 * m33 - m23 * m31)
        + (m00 * m13 - m03 * m10) * (m21 * m32 - m22 * m31)
        + (m01 * m12 - m02 * m11) * (m20 * m33 - m23 * m30)
        - (m01 * m13 - m03 * m11) * (m20 * m32 - m22 * m30)
        + (m02 * m13 - m03 * m12) * (m20 * m31 - m21 * m30);

  }

  /**
   * Sets the rotational component (upper 3x3) of this matrix to the matrix
   * values in the single precision Matrix3f argument; the other elements of
   * this matrix are initialized as if this were an identity matrix (ie, affine
   * matrix with no translational component).
   * 
   * @param m1
   *        the 3x3 matrix
   */
  public final void setM3(Matrix3f m1) {
    m00 = m1.m00;
    m01 = m1.m01;
    m02 = m1.m02;
    m03 = 0.0f;
    m10 = m1.m10;
    m11 = m1.m11;
    m12 = m1.m12;
    m13 = 0.0f;
    m20 = m1.m20;
    m21 = m1.m21;
    m22 = m1.m22;
    m23 = 0.0f;
    m30 = 0.0f;
    m31 = 0.0f;
    m32 = 0.0f;
    m33 = 1.0f;
  }

  /**
   * Sets the values in this Matrix4f equal to the row-major array parameter
   * (ie, the first four elements of the array will be copied into the first row
   * of this matrix, etc.).
   * 
   * @param m
   */
  public final void setA(float m[]) {
    m00 = m[0];
    m01 = m[1];
    m02 = m[2];
    m03 = m[3];
    m10 = m[4];
    m11 = m[5];
    m12 = m[6];
    m13 = m[7];
    m20 = m[8];
    m21 = m[9];
    m22 = m[10];
    m23 = m[11];
    m30 = m[12];
    m31 = m[13];
    m32 = m[14];
    m33 = m[15];
  }

  /**
   * Modifies the translational components of this matrix to the values of the
   * Vector3f argument; the other values of this matrix are not modified.
   * 
   * @param trans
   *        the translational component
   */
  public void setTranslation(Vector3f trans) {
    m03 = trans.x;
    m13 = trans.y;
    m23 = trans.z;
  }

  /**
   * Sets the value of this matrix to a rotation matrix about the x axis by the
   * passed angle.
   * 
   * @param angle
   *        the angle to rotate about the X axis in radians
   */
  public final void rotX(float angle) {
    float c = (float) Math.cos(angle);
    float s = (float) Math.sin(angle);
    m00 = 1.0f;
    m01 = 0.0f;
    m02 = 0.0f;
    m03 = 0.0f;
    m10 = 0.0f;
    m11 = c;
    m12 = -s;
    m13 = 0.0f;
    m20 = 0.0f;
    m21 = s;
    m22 = c;
    m23 = 0.0f;
    m30 = 0.0f;
    m31 = 0.0f;
    m32 = 0.0f;
    m33 = 1.0f;
  }

  /**
   * Sets the value of this matrix to a rotation matrix about the y axis by the
   * passed angle.
   * 
   * @param angle
   *        the angle to rotate about the Y axis in radians
   */
  public final void rotY(float angle) {
    float c = (float) Math.cos(angle);
    float s = (float) Math.sin(angle);
    m00 = c;
    m01 = 0.0f;
    m02 = s;
    m03 = 0.0f;
    m10 = 0.0f;
    m11 = 1.0f;
    m12 = 0.0f;
    m13 = 0.0f;
    m20 = -s;
    m21 = 0.0f;
    m22 = c;
    m23 = 0.0f;
    m30 = 0.0f;
    m31 = 0.0f;
    m32 = 0.0f;
    m33 = 1.0f;
  }

  /**
   * Sets the value of this matrix to a rotation matrix about the z axis by the
   * passed angle.
   * 
   * @param angle
   *        the angle to rotate about the Z axis in radians
   */
  public final void rotZ(float angle) {
    float c = (float) Math.cos(angle);
    float s = (float) Math.sin(angle);
    m00 = c;
    m01 = -s;
    m02 = 0.0f;
    m03 = 0.0f;
    m10 = s;
    m11 = c;
    m12 = 0.0f;
    m13 = 0.0f;
    m20 = 0.0f;
    m21 = 0.0f;
    m22 = 1.0f;
    m23 = 0.0f;
    m30 = 0.0f;
    m31 = 0.0f;
    m32 = 0.0f;
    m33 = 1.0f;
  }

  /**
   * Multiplies each element of this matrix by a scalar.
   * 
   * @param scalar
   *        The scalar multiplier.
   */
  public final void mul(float scalar) {
    m00 *= scalar;
    m01 *= scalar;
    m02 *= scalar;
    m03 *= scalar;
    m10 *= scalar;
    m11 *= scalar;
    m12 *= scalar;
    m13 *= scalar;
    m20 *= scalar;
    m21 *= scalar;
    m22 *= scalar;
    m23 *= scalar;
    m30 *= scalar;
    m31 *= scalar;
    m32 *= scalar;
    m33 *= scalar;
  }

  /**
   * Sets the value of this matrix to the result of multiplying itself with
   * matrix m1.
   * 
   * @param m1
   *        the other matrix
   */
  public final void mul(Matrix4f m1) {
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
  public final void mul2(Matrix4f m1, Matrix4f m2) {
    // alias-safe way.
    set(m1.m00 * m2.m00 + m1.m01 * m2.m10 + m1.m02 * m2.m20 + m1.m03 * m2.m30,
        m1.m00 * m2.m01 + m1.m01 * m2.m11 + m1.m02 * m2.m21 + m1.m03 * m2.m31,
        m1.m00 * m2.m02 + m1.m01 * m2.m12 + m1.m02 * m2.m22 + m1.m03 * m2.m32,
        m1.m00 * m2.m03 + m1.m01 * m2.m13 + m1.m02 * m2.m23 + m1.m03 * m2.m33,

        m1.m10 * m2.m00 + m1.m11 * m2.m10 + m1.m12 * m2.m20 + m1.m13 * m2.m30,
        m1.m10 * m2.m01 + m1.m11 * m2.m11 + m1.m12 * m2.m21 + m1.m13 * m2.m31,
        m1.m10 * m2.m02 + m1.m11 * m2.m12 + m1.m12 * m2.m22 + m1.m13 * m2.m32,
        m1.m10 * m2.m03 + m1.m11 * m2.m13 + m1.m12 * m2.m23 + m1.m13 * m2.m33,

        m1.m20 * m2.m00 + m1.m21 * m2.m10 + m1.m22 * m2.m20 + m1.m23 * m2.m30,
        m1.m20 * m2.m01 + m1.m21 * m2.m11 + m1.m22 * m2.m21 + m1.m23 * m2.m31,
        m1.m20 * m2.m02 + m1.m21 * m2.m12 + m1.m22 * m2.m22 + m1.m23 * m2.m32,
        m1.m20 * m2.m03 + m1.m21 * m2.m13 + m1.m22 * m2.m23 + m1.m23 * m2.m33,

        m1.m30 * m2.m00 + m1.m31 * m2.m10 + m1.m32 * m2.m20 + m1.m33 * m2.m30,
        m1.m30 * m2.m01 + m1.m31 * m2.m11 + m1.m32 * m2.m21 + m1.m33 * m2.m31,
        m1.m30 * m2.m02 + m1.m31 * m2.m12 + m1.m32 * m2.m22 + m1.m33 * m2.m32,
        m1.m30 * m2.m03 + m1.m31 * m2.m13 + m1.m32 * m2.m23 + m1.m33 * m2.m33);
  }

  /**
   * Returns true if the Object o is of type Matrix4f and all of the data
   * members of t1 are equal to the corresponding data members in this Matrix4f.
   * 
   * @param o
   *        the object with which the comparison is made.
   */
  @Override
  public boolean equals(Object o) {
    if (!(o instanceof Matrix4f))
      return false;
    Matrix4f m = (Matrix4f) o;
    return (this.m00 == m.m00 && this.m01 == m.m01 && this.m02 == m.m02
        && this.m03 == m.m03 && this.m10 == m.m10 && this.m11 == m.m11
        && this.m12 == m.m12 && this.m13 == m.m13 && this.m20 == m.m20
        && this.m21 == m.m21 && this.m22 == m.m22 && this.m23 == m.m23
        && this.m30 == m.m30 && this.m31 == m.m31 && this.m32 == m.m32 && this.m33 == m.m33);
  }

  /**
   * Returns a hash number based on the data values in this object. Two
   * different Matrix4f objects with identical data values (ie, returns true for
   * equals(Matrix4f) ) will return the same hash number. Two objects with
   * different data members may return the same hash value, although this is not
   * likely.
   * 
   * @return the integer hash value
   */
  @Override
  public int hashCode() {
    return Tuple3f.floatToIntBits0(m00) ^ Tuple3f.floatToIntBits0(m01)
        ^ Tuple3f.floatToIntBits0(m02) ^ Tuple3f.floatToIntBits0(m03)
        ^ Tuple3f.floatToIntBits0(m10) ^ Tuple3f.floatToIntBits0(m11)
        ^ Tuple3f.floatToIntBits0(m12) ^ Tuple3f.floatToIntBits0(m13)
        ^ Tuple3f.floatToIntBits0(m20) ^ Tuple3f.floatToIntBits0(m21)
        ^ Tuple3f.floatToIntBits0(m22) ^ Tuple3f.floatToIntBits0(m23)
        ^ Tuple3f.floatToIntBits0(m30) ^ Tuple3f.floatToIntBits0(m31)
        ^ Tuple3f.floatToIntBits0(m32) ^ Tuple3f.floatToIntBits0(m33);
  }

  /**
   * Transform the vector vec using this Matrix4f and place the result into
   * vecOut.
   * 
   * @param vec
   *        the single precision vector to be transformed
   * @param vecOut
   *        the vector into which the transformed values are placed
   */
  public final void transformT2(Tuple4f vec, Tuple4f vecOut) {
    // alias-safe
    vecOut.set(m00 * vec.x + m01 * vec.y + m02 * vec.z + m03 * vec.w, m10
        * vec.x + m11 * vec.y + m12 * vec.z + m13 * vec.w, m20 * vec.x + m21
        * vec.y + m22 * vec.z + m23 * vec.w, m30 * vec.x + m31 * vec.y + m32
        * vec.z + m33 * vec.w);
  }

  /**
   * Transform the vector vec using this Matrix4f and place the result back into
   * vec.
   * 
   * @param vec
   *        the single precision vector to be transformed
   */
  public final void transform4(Tuple4f vec) {
    transformT2(vec, vec);
  }

  /**
   * Transforms the point parameter with this Matrix4f and places the result
   * into pointOut. The fourth element of the point input paramter is assumed to
   * be one.
   * 
   * @param point
   *        the input point to be transformed.
   * @param pointOut
   *        the transformed point
   */
  public final void transform2(Point3f point, Point3f pointOut) {
    try {
      pointOut.set(m00 * point.x + m01 * point.y + m02 * point.z + m03, m10
          * point.x + m11 * point.y + m12 * point.z + m13, m20 * point.x + m21
          * point.y + m22 * point.z + m23);
    } catch (NullPointerException e) {
    }
  }

  /**
   * Transforms the point parameter with this Matrix4f and places the result
   * back into point. The fourth element of the point input paramter is assumed
   * to be one.
   * 
   * @param point
   *        the input point to be transformed.
   */
  public final void transform(Point3f point) {
    transform2(point, point);
  }

  /**
   * Transforms the normal parameter by this Matrix4f and places the value into
   * normalOut. The fourth element of the normal is assumed to be zero.
   * 
   * @param normal
   *        the input normal to be transformed.
   * @param normalOut
   *        the transformed normal
   */
  public final void transformV2(Vector3f normal, Vector3f normalOut) {
    normalOut.set(m00 * normal.x + m01 * normal.y + m02 * normal.z, m10
        * normal.x + m11 * normal.y + m12 * normal.z, m20 * normal.x + m21
        * normal.y + m22 * normal.z);
  }

  /**
   * Transforms the normal parameter by this transform and places the value back
   * into normal. The fourth element of the normal is assumed to be zero.
   * 
   * @param normal
   *        the input normal to be transformed.
   */
  public final void transformV(Vector3f normal) {
    transformV2(normal, normal);
  }

  /**
   * Sets this matrix to all zeros.
   */
  public final void setZero() {
    m00 = 0.0f;
    m01 = 0.0f;
    m02 = 0.0f;
    m03 = 0.0f;
    m10 = 0.0f;
    m11 = 0.0f;
    m12 = 0.0f;
    m13 = 0.0f;
    m20 = 0.0f;
    m21 = 0.0f;
    m22 = 0.0f;
    m23 = 0.0f;
    m30 = 0.0f;
    m31 = 0.0f;
    m32 = 0.0f;
    m33 = 0.0f;
  }

  /**
   * Sets 16 values
   * @param m00 
   * @param m01 
   * @param m02 
   * @param m03 
   * @param m10 
   * @param m11 
   * @param m12 
   * @param m13 
   * @param m20 
   * @param m21 
   * @param m22 
   * @param m23 
   * @param m30 
   * @param m31 
   * @param m32 
   * @param m33 
   */
  private void set(float m00, float m01, float m02, float m03, float m10,
                   float m11, float m12, float m13, float m20, float m21,
                   float m22, float m23, float m30, float m31, float m32,
                   float m33) {
    this.m00 = m00;
    this.m01 = m01;
    this.m02 = m02;
    this.m03 = m03;
    this.m10 = m10;
    this.m11 = m11;
    this.m12 = m12;
    this.m13 = m13;
    this.m20 = m20;
    this.m21 = m21;
    this.m22 = m22;
    this.m23 = m23;
    this.m30 = m30;
    this.m31 = m31;
    this.m32 = m32;
    this.m33 = m33;
  }

}
