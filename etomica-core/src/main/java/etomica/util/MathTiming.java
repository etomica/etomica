/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.util;

import etomica.math.SpecialFunctions;
import org.apache.commons.math3.special.Erf;

/**
 * This class benchmarks various mathematical functions.
 * erfc is also benchmarked in the top-level benchmarks module
 */
public class MathTiming {

    public static void main(String[] args) {
        long n = 100;
        if (args.length > 0) n = Long.parseLong(args[0]);
        n *= 1000000L;
        long nmult = 10L * n, nsin = n, nexp = n, nsqrt = 10L * n, npow = n, nerfc1 = n / 10L, nerfc2 = n;

        double sbase = 0, ssin = 0, sexp = 0, ssqrt = 0, spow = 0, serfc1 = 0, serfc2 = 0;

        double x = 0, dx = 3.0 / nmult;
        long t1 = System.nanoTime();
        for (long i = 0; i < nmult; i++) {
            sbase += x;
            x += dx;
        }
        double tbase = (System.nanoTime() - t1) / (double) nmult;

        x = 0;
        dx = 10.0 / nsin;
        t1 = System.nanoTime();
        for (long i = 0; i < nsin; i++) {
            ssin += Math.sin(x);
            x += dx;
        }
        double tsin = (System.nanoTime() - t1) / (double) nsin;

        x = -5;
        dx = 10.0 / nexp;
        t1 = System.nanoTime();
        for (long i = 0; i < nexp; i++) {
            sexp += Math.exp(x);
            x += dx;
        }
        double texp = (System.nanoTime() - t1) / (double) nexp;

        x = 0;
        dx = 10.0 / nsqrt;
        t1 = System.nanoTime();
        for (long i = 0; i < nsqrt; i++) {
            ssqrt += Math.sqrt(x);
            x += dx;
        }
        double tsqrt = (System.nanoTime() - t1) / (double) nsqrt;

        x = 0;
        dx = 10.0 / npow;
        t1 = System.nanoTime();
        for (long i = 0; i < npow; i++) {
            spow += Math.pow(x, 0.234);
            x += dx;
        }
        double tpow = (System.nanoTime() - t1) / (double) npow;

        x = 0;
        dx = 2.0 / nerfc1;
        t1 = System.nanoTime();
        for (long i = 0; i < nerfc1; i++) {
            serfc1 += Erf.erfc(x);
            x += dx;
        }
        double terfc1 = (System.nanoTime() - t1) / (double) nerfc1;

        x = 0;
        dx = 2.0 / nerfc2;
        t1 = System.nanoTime();
        for (long i = 0; i < nerfc2; i++) {
            serfc2 += SpecialFunctions.erfc(x);
            x += dx;
        }
        double terfc2 = (System.nanoTime() - t1) / (double) nerfc2;

        System.out.println("base: " + tbase);
        System.out.println("sin: " + tsin);
        System.out.println("exp: " + texp);
        System.out.println("sqrt: " + tsqrt);
        System.out.println("pow: " + tpow);
        System.out.println("erfc(apache): " + terfc1);
        System.out.println("erfc(approx): " + terfc2);
        // printed so that the compiler can't cheat
        System.out.println("sum base: " + sbase);
        System.out.println("sum sin: " + ssin);
        System.out.println("sum exp: " + sexp);
        System.out.println("sum sqrt: " + ssqrt);
        System.out.println("sum pow: " + spow);
        System.out.println("sum erfc(apache): " + serfc1);
        System.out.println("sum erfc(approx): " + serfc2);

    }
}
