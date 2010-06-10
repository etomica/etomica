package etomica.graph.model.impl;

import etomica.graph.model.Coefficient;

public class CoefficientImpl implements Coefficient {

  private int numerator;
  private int denominator;
  private int sign;

  public CoefficientImpl(int value) {

    this(value, 1, 1);
  }

  public CoefficientImpl(int numerator, int denominator) {

    this(numerator, denominator, 1);
  }

  public CoefficientImpl(int numerator, int denominator, int sgn) {

    this.numerator = numerator;
    this.denominator = denominator;
    sign = sgn;
    reduce();
  }

  public void add(Coefficient value) {

    if (getDenominator() == value.getDenominator()) {
      numerator = sign * numerator + value.getSign() * value.getNumerator();
    }
    else {
      numerator = sign * numerator * value.getDenominator() + value.getSign()
          * value.getNumerator() * denominator;
      denominator = getDenominator() * value.getDenominator();
    }
    sign = 1;
    if (numerator < 0) {
      numerator = -numerator;
      sign = -1;
    }
    reduce();
  }

  public Coefficient copy() {

    return new CoefficientImpl(numerator, denominator, sign);
  }

  public int getDenominator() {

    return denominator;
  }

  public int getNumerator() {

    return numerator;
  }

  public int getSign() {

    return sign;
  }

  public void inc() {

    numerator += denominator;
  }

  public void multiply(Coefficient value) {

    numerator = numerator * value.getNumerator();
    denominator = denominator * value.getDenominator();
    sign = sign * value.getSign();
    reduce();
  }

  protected void reduce() {
outer: do {
      int min = numerator;
      int max = denominator;
      if (min > max) {
        int t = min;
        min = max;
        max = t;
      }
      if (min < 2) {
        return;
      }
      // find factors for min, see if each factor is a factor of max
      int sqrt = (int)Math.sqrt(min);
      for (int i = 1; i<sqrt+1; i++) {
        int fac = min/i;
        if (fac*i == min) {
          // i and fac are factors of min
          int fac2 = max/fac;
          if (fac*fac2 == max) {
            // fac is a common factor
            numerator /= fac;
            denominator /= fac;
            // now start over
            reduce();
            continue outer;
          }
          if (i>1) {
            fac2 = max/i;
            if (fac2*i == max) {
              // i is a common factor
              numerator /= i;
              denominator /= i;
              // now start over
              reduce();
              continue outer;
            }
          }
        }
      }
      // we want one pass, with the option starting over from the middle
    } while (false);
  }

  public void setDenominator(int value) {

    denominator = value;
  }

  public void setNumerator(int value) {

    numerator = value;
  }

  public void switchSign() {

    sign = -sign;
  }

  @Override
  public String toString() {

    return (sign < 0 ? "-" : "") + numerator
        + (denominator == 1 ? "" : "/" + denominator);
  }
}
