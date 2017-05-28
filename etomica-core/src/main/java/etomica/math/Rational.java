/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math;

/**
 * Class handling rational numbers (numerator / denominator), including
 * addition, multiplication and reducing (15/6 => 5/2).
 *
 * @author Andrew Schultz
 */
public class Rational {

    public Rational() {
        this(1,1);
    }

    /**
     * Creates Rational instance with the given numerator and denominator.  If
     * appropriate, the rational is simplified (factoring).
     */
    public Rational(int numerator, int denominator) {
        this.numerator = numerator;
        this.denominator = denominator;
        simplify();
    }

    /**
     * Returns the numerator.
     */
    public int numerator() {
        return numerator;
    }

    /**
     * Returns the denominator.
     */
    public int denominator() {
        return denominator;
    }

    /**
     * Returns a new Rational instance, which is the result of multiplying this
     * Rational by the given Rational.
     */
    public Rational times(Rational r) {
        return new Rational(numerator*r.numerator, denominator*r.denominator);
    }

    /**
     * Returns a new Rational instance, which is the result of adding the given
     * Rational to this Rational.
     */
    public Rational plus(Rational r) {
        return new Rational(numerator*r.denominator+denominator*r.numerator, denominator*r.denominator);
    }

    /**
     * Simplify the rational.  If the denominator is negative, the signs of the
     * numerator and denominator are both flipped.  Numerator and denominator
     * are also divided by any common factors.
     */
    protected void simplify() {
        if (denominator < 0) {
            numerator = -numerator;
            denominator = -denominator;
        }

        int min = Math.abs(numerator);
        int max = Math.abs(denominator);
        if (min > max) {
            int t = min;
            min = max;
            max = t;
        }
        if (min < 2) {
            return;
        }
        int sqrt = (int)Math.sqrt(min);
        for (int i = 1; i<sqrt+1; i++) {
            int fac = min/i;
            if (fac*i == min) {
                int fac2 = max/fac;
                if (fac*fac2 == max) {
                    // fac is a common factor
                    numerator /= fac;
                    denominator /= fac;
                    // now start over
                    simplify();
                    return;
                }
                if (i>1) {
                    fac2 = max/i;
                    if (fac2*i == max) {
                        // i is a common factor
                        numerator /= i;
                        denominator /= i;
                        // now start over
                        simplify();
                        return;
                    }
                }
            }
        }
    }

    public String toString() {
        return denominator == 1 ? numerator+"" : numerator+"/"+denominator;
    }

    protected int numerator, denominator;
    public static final Rational ZERO = new Rational(0,1);
}
