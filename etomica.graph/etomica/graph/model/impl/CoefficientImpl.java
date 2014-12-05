/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.model.impl;

import etomica.graph.model.Coefficient;

public class CoefficientImpl implements Coefficient {

  private int numerator;
  private int denominator;
  private double dValue = Double.NaN;
  private static boolean hitOverflow = false;

  public CoefficientImpl(int value) {
    this(value, 1);
  }

  public CoefficientImpl(int numerator, int denominator) {
    this(numerator, denominator, Double.NaN);
  }

  protected CoefficientImpl(int numerator, int denominator, double dValue) {
    this.numerator = numerator;
    this.denominator = denominator;
    this.dValue = dValue;
    if (!Double.isNaN(dValue)) return;
    if (denominator < 0) {
      this.numerator = -numerator;
      this.denominator = -denominator;
    }
    reduce(this.numerator, this.denominator);
  }

  public boolean hasOverflow() {
    return !Double.isNaN(dValue);
  }

  public void add(Coefficient value) {
    if (!Double.isNaN(dValue) || value.hasOverflow()) {
      dValue += value.getValue();
      return;
    }

    long newNumerator = numerator, newDenominator = denominator;
    if (getDenominator() == value.getDenominator()) {
      newNumerator += value.getNumerator();
      newDenominator = denominator;
    }
    else {
      newNumerator *= value.getDenominator();
      newNumerator += value.getNumerator() * newDenominator;
      newDenominator *= value.getDenominator();
    }
    reduce(newNumerator, newDenominator);
  }

  public Coefficient copy() {
    return new CoefficientImpl(numerator, denominator, dValue);
  }

  public int getDenominator() {
    if (!Double.isNaN(dValue)) throw new RuntimeException("already broken");
    return denominator;
  }

  public int getNumerator() {
    if (!Double.isNaN(dValue)) throw new RuntimeException("already broken");
    return numerator;
  }

  public double getValue() {
    if (!Double.isNaN(dValue)) {
      return dValue;
    }

    return ((double)numerator)/denominator;
  }

  public void multiply(Coefficient value) {
    if (!Double.isNaN(dValue) || value.hasOverflow()) {
      dValue *= value.getValue();
      return;
    }
    long newNumerator = numerator;
    long newDenominator = denominator;
    reduce(newNumerator * value.getNumerator(), newDenominator * value.getDenominator());
  }

  public void divide(Coefficient value) {
    if (!Double.isNaN(dValue) || value.hasOverflow()) {
      dValue /= value.getValue();
      return;
    }

    long newNumerator = numerator;
    newNumerator *= value.getDenominator();
    long newDenominator = denominator;
    newDenominator *= value.getNumerator();
    if (denominator < 0) {
      newNumerator = -newNumerator;
      newDenominator = -newDenominator;
    }
    reduce(newNumerator, newDenominator);
  }

  protected void reduce(long newNumerator, long newDenominator) {
      if (newDenominator < 0) {
        newDenominator = -newDenominator;
        newNumerator = -newNumerator;
      }
      long min = Math.abs(newNumerator);
      long max = newDenominator;
      if (min > max) {
        long t = min;
        min = max;
        max = t;
      }
      if (min < 2) {
        if (newNumerator < Integer.MIN_VALUE || newNumerator > Integer.MAX_VALUE || newDenominator > Integer.MAX_VALUE) {
          dValue = (double)newNumerator / (double)newDenominator;
          if (!hitOverflow) System.err.println("integer overflow, using floating point");
          hitOverflow = true;
          return;
        }
        numerator = (int)newNumerator;
        denominator = (int)newDenominator;
        return;
      }
      // find factors for min, see if each factor is a factor of max
      int sqrt = (int)Math.sqrt(min);
      for (int i = 1; i<sqrt+1; i++) {
        long fac = min/i;
        if (fac*i == min) {
          // i and fac are factors of min
          long fac2 = max/fac;
          if (fac*fac2 == max) {
            // fac is a common factor
            newNumerator /= fac;
            newDenominator /= fac;
            // now start over
            reduce(newNumerator, newDenominator);
            return;
          }
          if (i>1) {
            fac2 = max/i;
            if (fac2*i == max) {
              // i is a common factor
              newNumerator /= i;
              newDenominator /= i;
              // now start over
              reduce(newNumerator, newDenominator);
              return;
            }
          }
        }
      }
      if (newNumerator < Integer.MIN_VALUE || newNumerator > Integer.MAX_VALUE || newDenominator > Integer.MAX_VALUE) {
        dValue = (double)newNumerator / (double)newDenominator;
        if (!hitOverflow) System.err.println("integer overflow, using floating point");
        hitOverflow = true;
        return;
      }
      numerator = (int)newNumerator;
      denominator = (int)newDenominator;
  }

  public void setDenominator(int value) {
    if (!Double.isNaN(dValue)) throw new RuntimeException("already broken");
    denominator = value;
  }

  public void setNumerator(int value) {
    if (!Double.isNaN(dValue)) throw new RuntimeException("already broken");
    numerator = value;
  }

  public String toString() {
    if (!Double.isNaN(dValue)) return dValue+"";
    return numerator + (denominator == 1 ? "" : "/" + denominator);
  }
}
