package etomica.graph.model.impl;

import etomica.graph.model.Coefficient;

public class CoefficientImpl implements Coefficient {

  private int numerator;
  private int denominator;

  public CoefficientImpl(int value) {

    this(value, 1);
  }

  public CoefficientImpl(int numerator, int denominator) {
    this.numerator = numerator;
    this.denominator = denominator;
    if (denominator < 0) {
      numerator = -numerator;
      denominator = -denominator;
    }
    reduce();
  }

  public void add(Coefficient value) {

    if (getDenominator() == value.getDenominator()) {
      numerator = numerator + value.getNumerator();
    }
    else {
      numerator = numerator * value.getDenominator() + 
          value.getNumerator() * denominator;
      denominator = getDenominator() * value.getDenominator();
    }
    reduce();
  }

  public Coefficient copy() {
    return new CoefficientImpl(numerator, denominator);
  }

  public int getDenominator() {

    return denominator;
  }

  public int getNumerator() {

    return numerator;
  }

  public void multiply(Coefficient value) {

    numerator = numerator * value.getNumerator();
    denominator = denominator * value.getDenominator();
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

  public String toString() {

    return numerator + (denominator == 1 ? "" : "/" + denominator);
  }
}
