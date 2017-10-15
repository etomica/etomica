/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import static org.junit.Assert.assertEquals;

/**
 * Created by alex on 4/24/17.
 */
public class SpecialFunctionsTest {
    private static final double EPSILON = 1e-10;
    @Rule
    public ExpectedException thrown = ExpectedException.none();
    private double epislon = 1E-7;

    @Test
    public void testErfc() throws Exception {
        for (double x = -10; x < 10.0001; x += 0.1) {
            assertEquals(org.apache.commons.math3.special.Erf.erfc(x), SpecialFunctions.erfc(x), 2E-7);
        }

    }

    @Test
    public void testLnGamma() throws Exception {
        for (double x = 0.1; x < 10.01; x += 0.1) {
            if (Math.abs(SpecialFunctions.lnGamma(x)) < 1.0) {
                assertEquals(org.apache.commons.math3.special.Gamma.logGamma(x), SpecialFunctions.lnGamma(x), 1E-14);
            } else {
                assertEquals(1.0, org.apache.commons.math3.special.Gamma.logGamma(x) / SpecialFunctions.lnGamma(x), 1E-14);
            }
        }
    }

    @Test
    public void testGamma() throws Exception {
        for (double x = 0.1; x < 10.01; x += 0.1) {
            if (Math.abs(SpecialFunctions.gamma(x)) < 1.0) {
                assertEquals(org.apache.commons.math3.special.Gamma.gamma(x), SpecialFunctions.gamma(x), 1E-14);
            } else {
                assertEquals(1.0, org.apache.commons.math3.special.Gamma.gamma(x) / SpecialFunctions.gamma(x), 1E-14);
            }
        }
    }

    @Test
    public void testGammaQ() throws Exception {
    }

    @Test
    public void testConfluentHypergeometric1F1() throws Exception {
        assertEquals(1.0, 1.32382847084462 / SpecialFunctions.confluentHypergeometric1F1(0.25, 1.5, 1.33333333333333), 1E-14);

        assertEquals(1.0, 1.18850394032265 / SpecialFunctions.confluentHypergeometric1F1(0.25, 1.5, 0.888888888888889), 1E-14);

        assertEquals(1.0, 1.11547393799773 / SpecialFunctions.confluentHypergeometric1F1(0.25, 1.5, 0.592592592592593), 1E-14);

        assertEquals(1.0, 1.07293902048305 / SpecialFunctions.confluentHypergeometric1F1(0.25, 1.5, 0.395061728395062), 1E-14);

        assertEquals(1.0, 1.04695720520900 / SpecialFunctions.confluentHypergeometric1F1(0.25, 1.5, 0.263374485596708), 1E-14);

        assertEquals(1.0, 1.03059829570291 / SpecialFunctions.confluentHypergeometric1F1(0.25, 1.5, 0.175582990397805), 1E-14);

        assertEquals(1.0, 1.02009476435734 / SpecialFunctions.confluentHypergeometric1F1(0.25, 1.5, 0.117055326931870), 1E-14);

        assertEquals(1.0, 1.01326419084245 / SpecialFunctions.confluentHypergeometric1F1(0.25, 1.5, 0.0780368846212468), 1E-14);

        assertEquals(1.0, 1.00878480732077 / SpecialFunctions.confluentHypergeometric1F1(0.25, 1.5, 0.0520245897474978), 1E-14);

        assertEquals(1.0, 1.00583100626732 / SpecialFunctions.confluentHypergeometric1F1(0.25, 1.5, 0.0346830598316652), 1E-14);

        assertEquals(0.018919655869471, SpecialFunctions.confluentHypergeometric1F1(-0.25, 0.5, 1.33333333333333), 1E-14);

        assertEquals(1.0, 0.432122622618571 / SpecialFunctions.confluentHypergeometric1F1(-0.25, 0.5, 0.888888888888889), 1E-14);

        assertEquals(1.0, 0.652955706161598 / SpecialFunctions.confluentHypergeometric1F1(-0.25, 0.5, 0.592592592592593), 1E-14);

        assertEquals(1.0, 0.781012238118942 / SpecialFunctions.confluentHypergeometric1F1(-0.25, 0.5, 0.395061728395062), 1E-14);

        assertEquals(1.0, 0.859080314343934 / SpecialFunctions.confluentHypergeometric1F1(-0.25, 0.5, 0.263374485596708), 1E-14);

        assertEquals(1.0, 0.908191340178865 / SpecialFunctions.confluentHypergeometric1F1(-0.25, 0.5, 0.175582990397805), 1E-14);

        assertEquals(1.0, 0.939711715868328 / SpecialFunctions.confluentHypergeometric1F1(-0.25, 0.5, 0.117055326931870), 1E-14);

        assertEquals(1.0, 0.960206262262829 / SpecialFunctions.confluentHypergeometric1F1(-0.25, 0.5, 0.0780368846212468), 1E-14);

        assertEquals(1.0, 0.973645236161197 / SpecialFunctions.confluentHypergeometric1F1(-0.25, 0.5, 0.0520245897474978), 1E-14);

        assertEquals(1.0, 0.98250688056046 / SpecialFunctions.confluentHypergeometric1F1(-0.25, 0.5, 0.0346830598316652), 1E-14);
    }

    @Test
    public void testCalcLegendrePolynomial() throws Exception {

        assertEquals(1.00000000000000, SpecialFunctions.calcLegendrePolynomial(0, -2.00000000000000), 1E-13);

        assertEquals(-2.00000000000000, SpecialFunctions.calcLegendrePolynomial(1, -2.00000000000000), 1E-13);

        assertEquals(5.50000000000000, SpecialFunctions.calcLegendrePolynomial(2, -2.00000000000000), 1E-13);

        assertEquals(-17.0000000000000, SpecialFunctions.calcLegendrePolynomial(3, -2.00000000000000), 1E-13);

        assertEquals(55.3750000000000, SpecialFunctions.calcLegendrePolynomial(4, -2.00000000000000), 1E-13);

        assertEquals(-1.75000000000000, SpecialFunctions.calcLegendrePolynomial(1, -1.75000000000000), 1E-13);

        assertEquals(4.09375000000000, SpecialFunctions.calcLegendrePolynomial(2, -1.75000000000000), 1E-13);

        assertEquals(-10.7734375000000, SpecialFunctions.calcLegendrePolynomial(3, -1.75000000000000), 1E-13);

        assertEquals(29.9233398437500, SpecialFunctions.calcLegendrePolynomial(4, -1.75000000000000), 1E-13);

        assertEquals(-1.50000000000000, SpecialFunctions.calcLegendrePolynomial(1, -1.50000000000000), 1E-13);

        assertEquals(2.87500000000000, SpecialFunctions.calcLegendrePolynomial(2, -1.50000000000000), 1E-13);

        assertEquals(-6.18750000000000, SpecialFunctions.calcLegendrePolynomial(3, -1.50000000000000), 1E-13);

        assertEquals(14.0859375000000, SpecialFunctions.calcLegendrePolynomial(4, -1.50000000000000), 1E-13);

        assertEquals(-1.25000000000000, SpecialFunctions.calcLegendrePolynomial(1, -1.25000000000000), 1E-13);

        assertEquals(1.84375000000000, SpecialFunctions.calcLegendrePolynomial(2, -1.25000000000000), 1E-13);

        assertEquals(-3.00781250000000, SpecialFunctions.calcLegendrePolynomial(3, -1.25000000000000), 1E-13);

        assertEquals(5.19677734375000, SpecialFunctions.calcLegendrePolynomial(4, -1.25000000000000), 1E-13);

        assertEquals(-1.00000000000000, SpecialFunctions.calcLegendrePolynomial(1, -1.00000000000000), 1E-13);

        assertEquals(1.00000000000000, SpecialFunctions.calcLegendrePolynomial(2, -1.00000000000000), 1E-13);

        assertEquals(-1.00000000000000, SpecialFunctions.calcLegendrePolynomial(3, -1.00000000000000), 1E-13);

        assertEquals(1.00000000000000, SpecialFunctions.calcLegendrePolynomial(4, -1.00000000000000), 1E-13);

        assertEquals(-0.750000000000000, SpecialFunctions.calcLegendrePolynomial(1, -0.750000000000000), 1E-13);

        assertEquals(0.343750000000000, SpecialFunctions.calcLegendrePolynomial(2, -0.750000000000000), 1E-13);

        assertEquals(0.0703125000000000, SpecialFunctions.calcLegendrePolynomial(3, -0.750000000000000), 1E-13);

        assertEquals(-0.350097656250000, SpecialFunctions.calcLegendrePolynomial(4, -0.750000000000000), 1E-13);

        assertEquals(-0.500000000000000, SpecialFunctions.calcLegendrePolynomial(1, -0.500000000000000), 1E-13);

        assertEquals(-0.125000000000000, SpecialFunctions.calcLegendrePolynomial(2, -0.500000000000000), 1E-13);

        assertEquals(0.437500000000000, SpecialFunctions.calcLegendrePolynomial(3, -0.500000000000000), 1E-13);

        assertEquals(-0.289062500000000, SpecialFunctions.calcLegendrePolynomial(4, -0.500000000000000), 1E-13);

        assertEquals(-0.250000000000000, SpecialFunctions.calcLegendrePolynomial(1, -0.250000000000000), 1E-13);

        assertEquals(-0.406250000000000, SpecialFunctions.calcLegendrePolynomial(2, -0.250000000000000), 1E-13);

        assertEquals(0.335937500000000, SpecialFunctions.calcLegendrePolynomial(3, -0.250000000000000), 1E-13);

        assertEquals(0.157714843750000, SpecialFunctions.calcLegendrePolynomial(4, -0.250000000000000), 1E-13);

        assertEquals(0, SpecialFunctions.calcLegendrePolynomial(1, 0), 1E-13);

        assertEquals(-0.500000000000000, SpecialFunctions.calcLegendrePolynomial(2, 0), 1E-13);

        assertEquals(0, SpecialFunctions.calcLegendrePolynomial(3, 0), 1E-13);

        assertEquals(0.375000000000000, SpecialFunctions.calcLegendrePolynomial(4, 0), 1E-13);

        assertEquals(0.250000000000000, SpecialFunctions.calcLegendrePolynomial(1, 0.250000000000000), 1E-13);

        assertEquals(-0.406250000000000, SpecialFunctions.calcLegendrePolynomial(2, 0.250000000000000), 1E-13);

        assertEquals(-0.335937500000000, SpecialFunctions.calcLegendrePolynomial(3, 0.250000000000000), 1E-13);

        assertEquals(0.157714843750000, SpecialFunctions.calcLegendrePolynomial(4, 0.250000000000000), 1E-13);

        assertEquals(0.500000000000000, SpecialFunctions.calcLegendrePolynomial(1, 0.500000000000000), 1E-13);

        assertEquals(-0.125000000000000, SpecialFunctions.calcLegendrePolynomial(2, 0.500000000000000), 1E-13);

        assertEquals(-0.437500000000000, SpecialFunctions.calcLegendrePolynomial(3, 0.500000000000000), 1E-13);

        assertEquals(-0.289062500000000, SpecialFunctions.calcLegendrePolynomial(4, 0.500000000000000), 1E-13);

        assertEquals(0.750000000000000, SpecialFunctions.calcLegendrePolynomial(1, 0.750000000000000), 1E-13);

        assertEquals(0.343750000000000, SpecialFunctions.calcLegendrePolynomial(2, 0.750000000000000), 1E-13);

        assertEquals(-0.0703125000000000, SpecialFunctions.calcLegendrePolynomial(3, 0.750000000000000), 1E-13);

        assertEquals(-0.350097656250000, SpecialFunctions.calcLegendrePolynomial(4, 0.750000000000000), 1E-13);

        assertEquals(1.00000000000000, SpecialFunctions.calcLegendrePolynomial(1, 1.00000000000000), 1E-13);

        assertEquals(1.00000000000000, SpecialFunctions.calcLegendrePolynomial(2, 1.00000000000000), 1E-13);

        assertEquals(1.00000000000000, SpecialFunctions.calcLegendrePolynomial(3, 1.00000000000000), 1E-13);

        assertEquals(1.00000000000000, SpecialFunctions.calcLegendrePolynomial(4, 1.00000000000000), 1E-13);

        assertEquals(1.25000000000000, SpecialFunctions.calcLegendrePolynomial(1, 1.25000000000000), 1E-13);

        assertEquals(1.84375000000000, SpecialFunctions.calcLegendrePolynomial(2, 1.25000000000000), 1E-13);

        assertEquals(3.00781250000000, SpecialFunctions.calcLegendrePolynomial(3, 1.25000000000000), 1E-13);

        assertEquals(5.19677734375000, SpecialFunctions.calcLegendrePolynomial(4, 1.25000000000000), 1E-13);

        assertEquals(1.50000000000000, SpecialFunctions.calcLegendrePolynomial(1, 1.50000000000000), 1E-13);

        assertEquals(2.87500000000000, SpecialFunctions.calcLegendrePolynomial(2, 1.50000000000000), 1E-13);

        assertEquals(6.18750000000000, SpecialFunctions.calcLegendrePolynomial(3, 1.50000000000000), 1E-13);

        assertEquals(14.0859375000000, SpecialFunctions.calcLegendrePolynomial(4, 1.50000000000000), 1E-13);

        assertEquals(1.75000000000000, SpecialFunctions.calcLegendrePolynomial(1, 1.75000000000000), 1E-13);

        assertEquals(4.09375000000000, SpecialFunctions.calcLegendrePolynomial(2, 1.75000000000000), 1E-13);

        assertEquals(10.7734375000000, SpecialFunctions.calcLegendrePolynomial(3, 1.75000000000000), 1E-13);

        assertEquals(29.9233398437500, SpecialFunctions.calcLegendrePolynomial(4, 1.75000000000000), 1E-13);

        assertEquals(2.00000000000000, SpecialFunctions.calcLegendrePolynomial(1, 2.00000000000000), 1E-13);

        assertEquals(5.50000000000000, SpecialFunctions.calcLegendrePolynomial(2, 2.00000000000000), 1E-13);

        assertEquals(17.0000000000000, SpecialFunctions.calcLegendrePolynomial(3, 2.00000000000000), 1E-13);

        assertEquals(55.3750000000000, SpecialFunctions.calcLegendrePolynomial(4, 2.00000000000000), 1E-13);

    }

    @Test
    public void testBessel() throws Exception {
        assertEquals(1.0, 0.1005008340281251 / SpecialFunctions.besselI( 1, 0.20000000000000000000), 1E-14);

        assertEquals(1.0, 0.005016687513894678 / SpecialFunctions.besselI( 2, 0.20000000000000000000), 1E-14);

        assertEquals(1.0, 0.0001670837502315642 / SpecialFunctions.besselI( 3, 0.20000000000000000000), 1E-14);

        assertEquals(1.0, 4.175006947752356E-6 / SpecialFunctions.besselI( 4, 0.20000000000000000000), 1E-14);

        assertEquals(1.0, 8.347232146991889E-8 / SpecialFunctions.besselI( 5, 0.20000000000000000000), 1E-14);

        assertEquals(1.0, 0.2040267557335706 / SpecialFunctions.besselI( 1, 0.40000000000000000000), 1E-14);

        assertEquals(1.0, 0.02026800356148826 / SpecialFunctions.besselI( 2, 0.40000000000000000000), 1E-14);

        assertEquals(1.0, 0.001346720118688000 / SpecialFunctions.besselI( 3, 0.40000000000000000000), 1E-14);

        assertEquals(1.0, 0.00006720178116825773 / SpecialFunctions.besselI( 4, 0.40000000000000000000), 1E-14);

        assertEquals(1.0, 2.684495322845460E-6 / SpecialFunctions.besselI( 5, 0.40000000000000000000), 1E-14);

        assertEquals(1.0, 0.3137040256049221 / SpecialFunctions.besselI( 1, 0.60000000000000000000), 1E-14);

        assertEquals(1.0, 0.04636527896759911 / SpecialFunctions.besselI( 2, 0.60000000000000000000), 1E-14);

        assertEquals(1.0, 0.004602165820928096 / SpecialFunctions.besselI( 3, 0.60000000000000000000), 1E-14);

        assertEquals(1.0, 0.0003436207583181480 / SpecialFunctions.besselI( 4, 0.60000000000000000000), 1E-14);

        assertEquals(1.0, 0.00002055571001945543 / SpecialFunctions.besselI( 5, 0.60000000000000000000), 1E-14);

        assertEquals(1.0, 0.4328648026206398 / SpecialFunctions.besselI( 1, 0.80000000000000000000), 1E-14);

        assertEquals(1.0, 0.08435291631820318 / SpecialFunctions.besselI( 2, 0.80000000000000000000), 1E-14);

        assertEquals(1.0, 0.01110022102962393 / SpecialFunctions.besselI( 3, 0.80000000000000000000), 1E-14);

        assertEquals(1.0, 0.001101258596023714 / SpecialFunctions.besselI( 4, 0.80000000000000000000), 1E-14);

        assertEquals(1.0, 0.00008763506938678689 / SpecialFunctions.besselI( 5, 0.80000000000000000000), 1E-14);

        assertEquals(1.0, 0.5651591039924850 / SpecialFunctions.besselI( 1, 1.0000000000000000000), 1E-14);

        assertEquals(1.0, 0.1357476697670383 / SpecialFunctions.besselI( 2, 1.0000000000000000000), 1E-14);

        assertEquals(1.0, 0.02216842492433190 / SpecialFunctions.besselI( 3, 1.0000000000000000000), 1E-14);

        assertEquals(1.0, 0.002737120221046866 / SpecialFunctions.besselI( 4, 1.0000000000000000000), 1E-14);

        assertEquals(1.0, 0.0002714631559569719 / SpecialFunctions.besselI( 5, 1.0000000000000000000), 1E-14);

        assertEquals(1.0, 0.7146779415526431 / SpecialFunctions.besselI( 1, 1.2000000000000000000), 1E-14);

        assertEquals(1.0, 0.2025956815463259 / SpecialFunctions.besselI( 2, 1.2000000000000000000), 1E-14);

        assertEquals(1.0, 0.03935900306489002 / SpecialFunctions.besselI( 3, 1.2000000000000000000), 1E-14);

        assertEquals(1.0, 0.005800666221875796 / SpecialFunctions.besselI( 4, 1.2000000000000000000), 1E-14);

        assertEquals(1.0, 0.0006878949190513823 / SpecialFunctions.besselI( 5, 1.2000000000000000000), 1E-14);

        assertEquals(1.0, 0.8860919814143274 / SpecialFunctions.besselI( 1, 1.4000000000000000000), 1E-14);

        assertEquals(1.0, 0.2875494119964631 / SpecialFunctions.besselI( 2, 1.4000000000000000000), 1E-14);

        assertEquals(1.0, 0.06452223285300407 / SpecialFunctions.besselI( 3, 1.4000000000000000000), 1E-14);

        assertEquals(1.0, 0.01102555691215997 / SpecialFunctions.besselI( 4, 1.4000000000000000000), 1E-14);

        assertEquals(1.0, 0.001519050497804235 / SpecialFunctions.besselI( 5, 1.4000000000000000000), 1E-14);

        assertEquals(1.0, 1.084810635129880 / SpecialFunctions.besselI( 1, 1.6000000000000000000), 1E-14);

        assertEquals(1.0, 0.3939673458265599 / SpecialFunctions.besselI( 2, 1.6000000000000000000), 1E-14);

        assertEquals(1.0, 0.09989227056347994 / SpecialFunctions.besselI( 3, 1.6000000000000000000), 1E-14);

        assertEquals(1.0, 0.01937133121351008 / SpecialFunctions.besselI( 4, 1.6000000000000000000), 1E-14);

        assertEquals(1.0, 0.003035614495929543 / SpecialFunctions.besselI( 5, 1.6000000000000000000), 1E-14);

        assertEquals(1.0, 1.317167230391899 / SpecialFunctions.besselI( 1, 1.8000000000000000000), 1E-14);

        assertEquals(1.0, 0.5260402117381632 / SpecialFunctions.besselI( 2, 1.8000000000000000000), 1E-14);

        assertEquals(1.0, 0.1481889820848698 / SpecialFunctions.besselI( 3, 1.8000000000000000000), 1E-14);

        assertEquals(1.0, 0.03207693812193060 / SpecialFunctions.besselI( 4, 1.8000000000000000000), 1E-14);

        assertEquals(1.0, 0.005624812654067089 / SpecialFunctions.besselI( 5, 1.8000000000000000000), 1E-14);

        assertEquals(1.0, 1.590636854637329 / SpecialFunctions.besselI( 1, 2.0000000000000000000), 1E-14);

        assertEquals(1.0, 0.6889484476987382 / SpecialFunctions.besselI( 2, 2.0000000000000000000), 1E-14);

        assertEquals(1.0, 0.2127399592398527 / SpecialFunctions.besselI( 3, 2.0000000000000000000), 1E-14);

        assertEquals(1.0, 0.05072856997918024 / SpecialFunctions.besselI( 4, 2.0000000000000000000), 1E-14);

        assertEquals(1.0, 0.009825679323131702 / SpecialFunctions.besselI( 5, 2.0000000000000000000), 1E-14);

        assertEquals(1.0, 1.914094650586386 / SpecialFunctions.besselI( 1, 2.2000000000000000000), 1E-14);

        assertEquals(1.0, 0.8890568175796904 / SpecialFunctions.besselI( 2, 2.2000000000000000000), 1E-14);

        assertEquals(1.0, 0.2976277095324036 / SpecialFunctions.besselI( 3, 2.2000000000000000000), 1E-14);

        assertEquals(1.0, 0.07734488249131686 / SpecialFunctions.besselI( 4, 2.2000000000000000000), 1E-14);

        assertEquals(1.0, 0.01637359138216051 / SpecialFunctions.besselI( 5, 2.2000000000000000000), 1E-14);

        assertEquals(1.0, 2.298123812543222 / SpecialFunctions.besselI( 1, 2.4000000000000000000), 1E-14);

        assertEquals(1.0, 1.134153480870062 / SpecialFunctions.besselI( 2, 2.4000000000000000000), 1E-14);

        assertEquals(1.0, 0.4078680110931191 / SpecialFunctions.besselI( 3, 2.4000000000000000000), 1E-14);

        assertEquals(1.0, 0.1144834531372640 / SpecialFunctions.besselI( 4, 2.4000000000000000000), 1E-14);

        assertEquals(1.0, 0.02625650063557234 / SpecialFunctions.besselI( 5, 2.4000000000000000000), 1E-14);

        assertEquals(1.0, 2.755384340504706 / SpecialFunctions.besselI( 1, 2.6000000000000000000), 1E-14);

        assertEquals(1.0, 1.433742488470821 / SpecialFunctions.besselI( 2, 2.6000000000000000000), 1E-14);

        assertEquals(1.0, 0.5496266659342133 / SpecialFunctions.besselI( 3, 2.6000000000000000000), 1E-14);

        assertEquals(1.0, 0.1653732593918667 / SpecialFunctions.besselI( 4, 2.6000000000000000000), 1E-14);

        assertEquals(1.0, 0.04078586780539262 / SpecialFunctions.besselI( 5, 2.6000000000000000000), 1E-14);

        assertEquals(1.0, 3.301055822635088 / SpecialFunctions.besselI( 1, 2.8000000000000000000), 1E-14);

        assertEquals(1.0, 1.799400687332901 / SpecialFunctions.besselI( 2, 2.8000000000000000000), 1E-14);

        assertEquals(1.0, 0.7304834121595154 / SpecialFunctions.besselI( 3, 2.8000000000000000000), 1E-14);

        assertEquals(1.0, 0.2340790898482246 / SpecialFunctions.besselI( 4, 2.8000000000000000000), 1E-14);

        assertEquals(1.0, 0.06168601259315954 / SpecialFunctions.besselI( 5, 2.8000000000000000000), 1E-14);

        assertEquals(1.0, 3.953370217402609 / SpecialFunctions.besselI( 1, 3.0000000000000000000), 1E-14);

        assertEquals(1.0, 2.245212440929951 / SpecialFunctions.besselI( 2, 3.0000000000000000000), 1E-14);

        assertEquals(1.0, 0.9597536294960079 / SpecialFunctions.besselI( 3, 3.0000000000000000000), 1E-14);

        assertEquals(1.0, 0.3257051819379354 / SpecialFunctions.besselI( 4, 3.0000000000000000000), 1E-14);

        assertEquals(1.0, 0.09120647766151335 / SpecialFunctions.besselI( 5, 3.0000000000000000000), 1E-14);

        assertEquals(1.0, 4.734253894709620 / SpecialFunctions.besselI( 1, 3.2000000000000000000), 1E-14);

        assertEquals(1.0, 2.788298502987037 / SpecialFunctions.besselI( 2, 3.2000000000000000000), 1E-14);

        assertEquals(1.0, 1.248880765975824 / SpecialFunctions.besselI( 3, 3.2000000000000000000), 1E-14);

        assertEquals(1.0, 0.4466470667823664 / SpecialFunctions.besselI( 4, 3.2000000000000000000), 1E-14);

        assertEquals(1.0, 0.1322630990199083 / SpecialFunctions.besselI( 5, 3.2000000000000000000), 1E-14);

        assertEquals(1.0, 5.670102192635220 / SpecialFunctions.besselI( 1, 3.4000000000000000000), 1E-14);

        assertEquals(1.0, 3.449458929469693 / SpecialFunctions.besselI( 2, 3.4000000000000000000), 1E-14);

        assertEquals(1.0, 1.611915216788522 / SpecialFunctions.besselI( 3, 3.4000000000000000000), 1E-14);

        assertEquals(1.0, 0.6049026645487712 / SpecialFunctions.besselI( 4, 3.4000000000000000000), 1E-14);

        assertEquals(1.0, 0.1886148296149430 / SpecialFunctions.besselI( 5, 3.4000000000000000000), 1E-14);

        assertEquals(1.0, 6.792714601361299 / SpecialFunctions.besselI( 1, 3.6000000000000000000), 1E-14);

        assertEquals(1.0, 4.253954212964399 / SpecialFunctions.besselI( 2, 3.6000000000000000000), 1E-14);

        assertEquals(1.0, 2.066098809178633 / SpecialFunctions.besselI( 3, 3.6000000000000000000), 1E-14);

        assertEquals(1.0, 0.8104561976666769 / SpecialFunctions.besselI( 4, 3.6000000000000000000), 1E-14);

        assertEquals(1.0, 0.2650850365860180 / SpecialFunctions.besselI( 5, 3.6000000000000000000), 1E-14);

        assertEquals(1.0, 8.140424578907956 / SpecialFunctions.besselI( 1, 3.8000000000000000000), 1E-14);

        assertEquals(1.0, 5.232454037200033 / SpecialFunctions.besselI( 2, 3.8000000000000000000), 1E-14);

        assertEquals(1.0, 2.632578223960553 / SpecialFunctions.besselI( 3, 3.8000000000000000000), 1E-14);

        assertEquals(1.0, 1.075751578314950 / SpecialFunctions.besselI( 4, 3.8000000000000000000), 1E-14);

        assertEquals(1.0, 0.3678380590869744 / SpecialFunctions.besselI( 5, 3.8000000000000000000), 1E-14);

        assertEquals(1.0, 9.759465153704450 / SpecialFunctions.besselI( 1, 4.0000000000000000000), 1E-14);

        assertEquals(1.0, 6.422189375284106 / SpecialFunctions.besselI( 2, 4.0000000000000000000), 1E-14);

        assertEquals(1.0, 3.337275778420344 / SpecialFunctions.besselI( 3, 4.0000000000000000000), 1E-14);

        assertEquals(1.0, 1.416275707653589 / SpecialFunctions.besselI( 4, 4.0000000000000000000), 1E-14);

        assertEquals(1.0, 0.5047243631131664 / SpecialFunctions.besselI( 5, 4.0000000000000000000), 1E-14);

        assertEquals(1.0, 11.70562014305162 / SpecialFunctions.besselI( 1, 4.2000000000000000000), 1E-14);

        assertEquals(1.0, 7.868351333273067 / SpecialFunctions.besselI( 2, 4.2000000000000000000), 1E-14);

        assertEquals(1.0, 4.211952206601076 / SpecialFunctions.besselI( 3, 4.2000000000000000000), 1E-14);

        assertEquals(1.0, 1.851276752414387 / SpecialFunctions.besselI( 4, 4.2000000000000000000), 1E-14);

        assertEquals(1.0, 0.6857107734308141 / SpecialFunctions.besselI( 5, 4.2000000000000000000), 1E-14);

        assertEquals(1.0, 14.04622133753311 / SpecialFunctions.besselI( 1, 4.4000000000000000000), 1E-14);

        assertEquals(1.0, 9.625789462431949 / SpecialFunctions.besselI( 2, 4.4000000000000000000), 1E-14);

        assertEquals(1.0, 5.295503644413152 / SpecialFunctions.besselI( 3, 4.4000000000000000000), 1E-14);

        assertEquals(1.0, 2.404648129141286 / SpecialFunctions.besselI( 4, 4.4000000000000000000), 1E-14);

        assertEquals(1.0, 0.9234161368835410 / SpecialFunctions.besselI( 5, 4.4000000000000000000), 1E-14);

        assertEquals(1.0, 16.86256476197666 / SpecialFunctions.besselI( 1, 4.6000000000000000000), 1E-14);

        assertEquals(1.0, 11.76107358300787 / SpecialFunctions.besselI( 2, 4.6000000000000000000), 1E-14);

        assertEquals(1.0, 6.635544255013292 / SpecialFunctions.besselI( 3, 4.6000000000000000000), 1E-14);

        assertEquals(1.0, 3.106015859077489 / SpecialFunctions.besselI( 4, 4.6000000000000000000), 1E-14);

        assertEquals(1.0, 1.233777543574181 / SpecialFunctions.besselI( 5, 4.6000000000000000000), 1E-14);

        assertEquals(1.0, 20.25283460023856 / SpecialFunctions.besselI( 1, 4.8000000000000000000), 1E-14);

        assertEquals(1.0, 14.35499690967306 / SpecialFunctions.besselI( 2, 4.8000000000000000000), 1E-14);

        assertEquals(1.0, 8.290337175511006 / SpecialFunctions.besselI( 3, 4.8000000000000000000), 1E-14);

        assertEquals(1.0, 3.992075440284307 / SpecialFunctions.besselI( 4, 4.8000000000000000000), 1E-14);

        assertEquals(1.0, 1.636878108370495 / SpecialFunctions.besselI( 5, 4.8000000000000000000), 1E-14);

        assertEquals(1.0, 24.33564214245053 / SpecialFunctions.besselI( 1, 5.0000000000000000000), 1E-14);

        assertEquals(1.0, 17.50561496662424 / SpecialFunctions.besselI( 2, 5.0000000000000000000), 1E-14);

        assertEquals(1.0, 10.33115016915114 / SpecialFunctions.besselI( 3, 5.0000000000000000000), 1E-14);

        assertEquals(1.0, 5.108234763642870 / SpecialFunctions.besselI( 4, 5.0000000000000000000), 1E-14);

        assertEquals(1.0, 2.157974547322546 / SpecialFunctions.besselI( 5, 5.0000000000000000000), 1E-14);

        assertEquals(1.0, 29.25430988179835 / SpecialFunctions.besselI( 1, 5.2000000000000000000), 1E-14);

        assertEquals(1.0, 21.33193506376818 / SpecialFunctions.besselI( 2, 5.2000000000000000000), 1E-14);

        assertEquals(1.0, 12.84512906351513 / SpecialFunctions.besselI( 3, 5.2000000000000000000), 1E-14);

        assertEquals(1.0, 6.510632298173797 / SpecialFunctions.besselI( 4, 5.2000000000000000000), 1E-14);

        assertEquals(1.0, 2.828771681709292 / SpecialFunctions.besselI( 5, 5.2000000000000000000), 1E-14);

        assertEquals(1.0, 35.18205850608358 / SpecialFunctions.besselI( 1, 5.4000000000000000000), 1E-14);

        assertEquals(1.0, 25.97839574633562 / SpecialFunctions.besselI( 2, 5.4000000000000000000), 1E-14);

        assertEquals(1.0, 15.93880239768683 / SpecialFunctions.besselI( 3, 5.4000000000000000000), 1E-14);

        assertEquals(1.0, 8.268615304461366 / SpecialFunctions.besselI( 4, 5.4000000000000000000), 1E-14);

        assertEquals(1.0, 3.689001946632952 / SpecialFunctions.besselI( 5, 5.4000000000000000000), 1E-14);

        assertEquals(1.0, 42.32828803246685 / SpecialFunctions.besselI( 1, 5.6000000000000000000), 1E-14);

        assertEquals(1.0, 31.62030556675627 / SpecialFunctions.besselI( 2, 5.6000000000000000000), 1E-14);

        assertEquals(1.0, 19.74235548478380 / SpecialFunctions.besselI( 3, 5.6000000000000000000), 1E-14);

        assertEquals(1.0, 10.46778183305934 / SpecialFunctions.besselI( 4, 5.6000000000000000000), 1E-14);

        assertEquals(1.0, 4.788381437556167 / SpecialFunctions.besselI( 5, 5.6000000000000000000), 1E-14);

        assertEquals(1.0, 50.94618497877481 / SpecialFunctions.besselI( 1, 5.8000000000000000000), 1E-14);

        assertEquals(1.0, 38.47044689994190 / SpecialFunctions.besselI( 2, 5.8000000000000000000), 1E-14);

        assertEquals(1.0, 24.41484228915970 / SpecialFunctions.besselI( 3, 5.8000000000000000000), 1E-14);

        assertEquals(1.0, 13.21371349736290 / SpecialFunctions.besselI( 4, 5.8000000000000000000), 1E-14);

        assertEquals(1.0, 6.189030568659158 / SpecialFunctions.besselI( 5, 5.8000000000000000000), 1E-14);

        assertEquals(1.0, 61.34193677764024 / SpecialFunctions.besselI( 1, 6.0000000000000000000), 1E-14);

        assertEquals(1.0, 46.78709471726456 / SpecialFunctions.besselI( 2, 6.0000000000000000000), 1E-14);

        assertEquals(1.0, 30.15054029946386 / SpecialFunctions.besselI( 3, 6.0000000000000000000), 1E-14);

        assertEquals(1.0, 16.63655441780070 / SpecialFunctions.besselI( 4, 6.0000000000000000000), 1E-14);

        assertEquals(1.0, 7.968467742396263 / SpecialFunctions.besselI( 5, 6.0000000000000000000), 1E-14);


        assertEquals(0.09950083263923600, SpecialFunctions.besselJ( 1, 0.20000000000000000000), 1E-13);

        assertEquals(0.004983354152783563, SpecialFunctions.besselJ( 2, 0.20000000000000000000), 1E-13);

        assertEquals(0.0001662504164352678, SpecialFunctions.besselJ( 3, 0.20000000000000000000), 1E-13);

        assertEquals(4.158340274471933E-6, SpecialFunctions.besselJ( 4, 0.20000000000000000000), 1E-13);

        assertEquals(8.319454360946915E-8, SpecialFunctions.besselJ( 5, 0.20000000000000000000), 1E-13);

        assertEquals(0.1960265779553187, SpecialFunctions.besselJ( 1, 0.40000000000000000000), 1E-13);

        assertEquals(0.01973466311703027, SpecialFunctions.besselJ( 2, 0.40000000000000000000), 1E-13);

        assertEquals(0.001320053214983958, SpecialFunctions.besselJ( 3, 0.40000000000000000000), 1E-13);

        assertEquals(0.00006613510772909677, SpecialFunctions.besselJ( 4, 0.40000000000000000000), 1E-13);

        assertEquals(2.648939597977585E-6, SpecialFunctions.besselJ( 5, 0.40000000000000000000), 1E-13);

        assertEquals(0.2867009880639157, SpecialFunctions.besselJ( 1, 0.60000000000000000000), 1E-13);

        assertEquals(0.04366509671584169, SpecialFunctions.besselJ( 2, 0.60000000000000000000), 1E-13);

        assertEquals(0.004399656708362193, SpecialFunctions.besselJ( 3, 0.60000000000000000000), 1E-13);

        assertEquals(0.0003314703677802370, SpecialFunctions.besselJ( 4, 0.60000000000000000000), 1E-13);

        assertEquals(0.00001994819537430024, SpecialFunctions.besselJ( 5, 0.60000000000000000000), 1E-13);

        assertEquals(0.3688420460941700, SpecialFunctions.besselJ( 1, 0.80000000000000000000), 1E-13);

        assertEquals(0.07581776248494472, SpecialFunctions.besselJ( 2, 0.80000000000000000000), 1E-13);

        assertEquals(0.01024676633055360, SpecialFunctions.besselJ( 3, 0.80000000000000000000), 1E-13);

        assertEquals(0.001032984994207302, SpecialFunctions.besselJ( 4, 0.80000000000000000000), 1E-13);

        assertEquals(0.00008308361151942143, SpecialFunctions.besselJ( 5, 0.80000000000000000000), 1E-13);

        assertEquals(0.4400505857449335, SpecialFunctions.besselJ( 1, 1.0000000000000000000), 1E-13);

        assertEquals(0.1149034849319005, SpecialFunctions.besselJ( 2, 1.0000000000000000000), 1E-13);

        assertEquals(0.01956335398266841, SpecialFunctions.besselJ( 3, 1.0000000000000000000), 1E-13);

        assertEquals(0.002476638964109955, SpecialFunctions.besselJ( 4, 1.0000000000000000000), 1E-13);

        assertEquals(0.0002497577302112344, SpecialFunctions.besselJ( 5, 1.0000000000000000000), 1E-13);

        assertEquals(0.4982890575672155, SpecialFunctions.besselJ( 1, 1.2000000000000000000), 1E-13);

        assertEquals(0.1593490183476631, SpecialFunctions.besselJ( 2, 1.2000000000000000000), 1E-13);

        assertEquals(0.03287433692499494, SpecialFunctions.besselJ( 3, 1.2000000000000000000), 1E-13);

        assertEquals(0.005022666277311587, SpecialFunctions.besselJ( 4, 1.2000000000000000000), 1E-13);

        assertEquals(0.0006101049237489684, SpecialFunctions.besselJ( 5, 1.2000000000000000000), 1E-13);

        assertEquals(0.5419477139308545, SpecialFunctions.besselJ( 1, 1.4000000000000000000), 1E-13);

        assertEquals(0.2073558995269320, SpecialFunctions.besselJ( 2, 1.4000000000000000000), 1E-13);

        assertEquals(0.05049771328895130, SpecialFunctions.besselJ( 3, 1.4000000000000000000), 1E-13);

        assertEquals(0.009062871711430658, SpecialFunctions.besselJ( 4, 1.4000000000000000000), 1E-13);

        assertEquals(0.001290125062081034, SpecialFunctions.besselJ( 5, 1.4000000000000000000), 1E-13);

        assertEquals(0.5698959352616804, SpecialFunctions.besselJ( 1, 1.6000000000000000000), 1E-13);

        assertEquals(0.2569677514377197, SpecialFunctions.besselJ( 2, 1.6000000000000000000), 1E-13);

        assertEquals(0.07252344333261900, SpecialFunctions.besselJ( 3, 1.6000000000000000000), 1E-13);

        assertEquals(0.01499516105960151, SpecialFunctions.besselJ( 4, 1.6000000000000000000), 1E-13);

        assertEquals(0.002452361965388557, SpecialFunctions.besselJ( 5, 1.6000000000000000000), 1E-13);

        assertEquals(0.5815169517311652, SpecialFunctions.besselJ( 1, 1.8000000000000000000), 1E-13);

        assertEquals(0.3061435353254030, SpecialFunctions.besselJ( 2, 1.8000000000000000000), 1E-13);

        assertEquals(0.09880201565861918, SpecialFunctions.besselJ( 3, 1.8000000000000000000), 1E-13);

        assertEquals(0.02319651686999431, SpecialFunctions.besselJ( 4, 1.8000000000000000000), 1E-13);

        assertEquals(0.004293614874688868, SpecialFunctions.besselJ( 5, 1.8000000000000000000), 1E-13);

        assertEquals(0.5767248077568734, SpecialFunctions.besselJ( 1, 2.0000000000000000000), 1E-13);

        assertEquals(0.3528340286156377, SpecialFunctions.besselJ( 2, 2.0000000000000000000), 1E-13);

        assertEquals(0.1289432494744021, SpecialFunctions.besselJ( 3, 2.0000000000000000000), 1E-13);

        assertEquals(0.03399571980756843, SpecialFunctions.besselJ( 4, 2.0000000000000000000), 1E-13);

        assertEquals(0.007039629755871685, SpecialFunctions.besselJ( 5, 2.0000000000000000000), 1E-13);

        assertEquals(0.5559630498190639, SpecialFunctions.besselJ( 1, 2.2000000000000000000), 1E-13);

        assertEquals(0.3950586874587933, SpecialFunctions.besselJ( 2, 2.2000000000000000000), 1E-13);

        assertEquals(0.1623254728332875, SpecialFunctions.besselJ( 3, 2.2000000000000000000), 1E-13);

        assertEquals(0.04764714754108161, SpecialFunctions.besselJ( 4, 2.2000000000000000000), 1E-13);

        assertEquals(0.01093688186155476, SpecialFunctions.besselJ( 5, 2.2000000000000000000), 1E-13);

        assertEquals(0.5201852681819310, SpecialFunctions.besselJ( 1, 2.4000000000000000000), 1E-13);

        assertEquals(0.4309800401876987, SpecialFunctions.besselJ( 2, 2.4000000000000000000), 1E-13);

        assertEquals(0.1981147987975668, SpecialFunctions.besselJ( 3, 2.4000000000000000000), 1E-13);

        assertEquals(0.06430695680621835, SpecialFunctions.besselJ( 4, 2.4000000000000000000), 1E-13);

        assertEquals(0.01624172388982766, SpecialFunctions.besselJ( 5, 2.4000000000000000000), 1E-13);

        assertEquals(0.4708182665175787, SpecialFunctions.besselJ( 1, 2.6000000000000000000), 1E-13);

        assertEquals(0.4589728517182526, SpecialFunctions.besselJ( 2, 2.6000000000000000000), 1E-13);

        assertEquals(0.2352938130489638, SpecialFunctions.besselJ( 3, 2.6000000000000000000), 1E-13);

        assertEquals(0.08401287070243310, SpecialFunctions.besselJ( 4, 2.6000000000000000000), 1E-13);

        assertEquals(0.02320732757390727, SpecialFunctions.besselJ( 5, 2.6000000000000000000), 1E-13);

        assertEquals(0.4097092468522887, SpecialFunctions.besselJ( 1, 2.8000000000000000000), 1E-13);

        assertEquals(0.4776854954017364, SpecialFunctions.besselJ( 2, 2.8000000000000000000), 1E-13);

        assertEquals(0.2726986037216204, SpecialFunctions.besselJ( 3, 2.8000000000000000000), 1E-13);

        assertEquals(0.1066686554303074, SpecialFunctions.besselJ( 4, 2.8000000000000000000), 1E-13);

        assertEquals(0.03206898322211490, SpecialFunctions.besselJ( 5, 2.8000000000000000000), 1E-13);

        assertEquals(0.3390589585259365, SpecialFunctions.besselJ( 1, 3.0000000000000000000), 1E-13);

        assertEquals(0.4860912605858911, SpecialFunctions.besselJ( 2, 3.0000000000000000000), 1E-13);

        assertEquals(0.3090627222552516, SpecialFunctions.besselJ( 3, 3.0000000000000000000), 1E-13);

        assertEquals(0.1320341839246122, SpecialFunctions.besselJ( 4, 3.0000000000000000000), 1E-13);

        assertEquals(0.04302843487704758, SpecialFunctions.besselJ( 5, 3.0000000000000000000), 1E-13);

        assertEquals(0.2613432487805048, SpecialFunctions.besselJ( 1, 3.2000000000000000000), 1E-13);

        assertEquals(0.4835277001449384, SpecialFunctions.besselJ( 2, 3.2000000000000000000), 1E-13);

        assertEquals(0.3430663764006682, SpecialFunctions.besselJ( 3, 3.2000000000000000000), 1E-13);

        assertEquals(0.1597217556063144, SpecialFunctions.besselJ( 4, 3.2000000000000000000), 1E-13);

        assertEquals(0.05623801261511791, SpecialFunctions.besselJ( 5, 3.2000000000000000000), 1E-13);

        assertEquals(0.1792258516815071, SpecialFunctions.besselJ( 1, 3.4000000000000000000), 1E-13);

        assertEquals(0.4697225683393576, SpecialFunctions.besselJ( 2, 3.4000000000000000000), 1E-13);

        assertEquals(0.3733889346000901, SpecialFunctions.besselJ( 3, 3.4000000000000000000), 1E-13);

        assertEquals(0.1891990809549190, SpecialFunctions.besselJ( 4, 3.4000000000000000000), 1E-13);

        assertEquals(0.07178537352913107, SpecialFunctions.besselJ( 5, 3.4000000000000000000), 1E-13);

        assertEquals(0.09546554717787640, SpecialFunctions.besselJ( 1, 3.6000000000000000000), 1E-13);

        assertEquals(0.4448053987996180, SpecialFunctions.besselJ( 2, 3.6000000000000000000), 1E-13);

        assertEquals(0.3987626737105880, SpecialFunctions.besselJ( 3, 3.6000000000000000000), 1E-13);

        assertEquals(0.2197990573846954, SpecialFunctions.besselJ( 4, 3.6000000000000000000), 1E-13);

        assertEquals(0.08967967603317951, SpecialFunctions.besselJ( 5, 3.6000000000000000000), 1E-13);

        assertEquals(0.01282100292673163, SpecialFunctions.besselJ( 1, 3.8000000000000000000), 1E-13);

        assertEquals(0.4093043064557913, SpecialFunctions.besselJ( 2, 3.8000000000000000000), 1E-13);

        assertEquals(0.4180256354477856, SpecialFunctions.besselJ( 3, 3.8000000000000000000), 1E-13);

        assertEquals(0.2507361705670280, SpecialFunctions.besselJ( 4, 3.8000000000000000000), 1E-13);

        assertEquals(0.1098399867985891, SpecialFunctions.besselJ( 5, 3.8000000000000000000), 1E-13);

        assertEquals(-0.06604332802354914, SpecialFunctions.besselJ( 1, 4.0000000000000000000), 1E-13);

        assertEquals(0.3641281458520728, SpecialFunctions.besselJ( 2, 4.0000000000000000000), 1E-13);

        assertEquals(0.4301714738756219, SpecialFunctions.besselJ( 3, 4.0000000000000000000), 1E-13);

        assertEquals(0.2811290649613601, SpecialFunctions.besselJ( 4, 4.0000000000000000000), 1E-13);

        assertEquals(0.1320866560470983, SpecialFunctions.besselJ( 5, 4.0000000000000000000), 1E-13);

        assertEquals(-0.1386469421260462, SpecialFunctions.besselJ( 1, 4.2000000000000000000), 1E-13);

        assertEquals(0.3105347009742123, SpecialFunctions.besselJ( 2, 4.2000000000000000000), 1E-13);

        assertEquals(0.4343942763872008, SpecialFunctions.besselJ( 3, 4.2000000000000000000), 1E-13);

        assertEquals(0.3100285510075031, SpecialFunctions.besselJ( 4, 4.2000000000000000000), 1E-13);

        assertEquals(0.1561362969604241, SpecialFunctions.besselJ( 5, 4.2000000000000000000), 1E-13);

        assertEquals(-0.2027755219230866, SpecialFunctions.besselJ( 1, 4.4000000000000000000), 1E-13);

        assertEquals(0.2500860982206644, SpecialFunctions.besselJ( 2, 4.4000000000000000000), 1E-13);

        assertEquals(0.4301265203055088, SpecialFunctions.besselJ( 3, 4.4000000000000000000), 1E-13);

        assertEquals(0.3364500658323021, SpecialFunctions.besselJ( 4, 4.4000000000000000000), 1E-13);

        assertEquals(0.1816008721168587, SpecialFunctions.besselJ( 5, 4.4000000000000000000), 1E-13);

        assertEquals(-0.2565528360974446, SpecialFunctions.besselJ( 1, 4.6000000000000000000), 1E-13);

        assertEquals(0.1845931052274261, SpecialFunctions.besselJ( 2, 4.6000000000000000000), 1E-13);

        assertEquals(0.4170685797734673, SpecialFunctions.besselJ( 3, 4.6000000000000000000), 1E-13);

        assertEquals(0.3594093901292703, SpecialFunctions.besselJ( 4, 4.6000000000000000000), 1E-13);

        assertEquals(0.2079912291470029, SpecialFunctions.besselJ( 5, 4.6000000000000000000), 1E-13);

        assertEquals(-0.2984998580995579, SpecialFunctions.besselJ( 1, 4.8000000000000000000), 1E-13);

        assertEquals(0.1160503864163677, SpecialFunctions.besselJ( 2, 4.8000000000000000000), 1E-13);

        assertEquals(0.3952085134465309, SpecialFunctions.besselJ( 3, 4.8000000000000000000), 1E-13);

        assertEquals(0.3779602553917960, SpecialFunctions.besselJ( 4, 4.8000000000000000000), 1E-13);

        assertEquals(0.2347252455397957, SpecialFunctions.besselJ( 5, 4.8000000000000000000), 1E-13);

        assertEquals(-0.3275791375914652, SpecialFunctions.besselJ( 1, 5.0000000000000000000), 1E-13);

        assertEquals(0.04656511627775222, SpecialFunctions.besselJ( 2, 5.0000000000000000000), 1E-13);

        assertEquals(0.3648312306136670, SpecialFunctions.besselJ( 3, 5.0000000000000000000), 1E-13);

        assertEquals(0.3912323604586482, SpecialFunctions.besselJ( 4, 5.0000000000000000000), 1E-13);

        assertEquals(0.2611405461201701, SpecialFunctions.besselJ( 5, 5.0000000000000000000), 1E-13);

        assertEquals(-0.3432230058719219, SpecialFunctions.besselJ( 1, 5.2000000000000000000), 1E-13);

        assertEquals(-0.02171840862129112, SpecialFunctions.besselJ( 2, 5.2000000000000000000), 1E-13);

        assertEquals(0.3265165377016980, SpecialFunctions.besselJ( 3, 5.2000000000000000000), 1E-13);

        assertEquals(0.3984682598155580, SpecialFunctions.besselJ( 4, 5.2000000000000000000), 1E-13);

        assertEquals(0.2865115543222374, SpecialFunctions.besselJ( 5, 5.2000000000000000000), 1E-13);

        assertEquals(-0.3453447907795863, SpecialFunctions.besselJ( 1, 5.4000000000000000000), 1E-13);

        assertEquals(-0.08669537682152215, SpecialFunctions.besselJ( 2, 5.4000000000000000000), 1E-13);

        assertEquals(0.2811259931340144, SpecialFunctions.besselJ( 3, 5.4000000000000000000), 1E-13);

        assertEquals(0.3990575914148714, SpecialFunctions.besselJ( 4, 5.4000000000000000000), 1E-13);

        assertEquals(0.3100704385917211, SpecialFunctions.besselJ( 5, 5.4000000000000000000), 1E-13);

        assertEquals(-0.3343328362910075, SpecialFunctions.besselJ( 1, 5.6000000000000000000), 1E-13);

        assertEquals(-0.1463754690747600, SpecialFunctions.besselJ( 2, 5.6000000000000000000), 1E-13);

        assertEquals(0.2297789298090361, SpecialFunctions.besselJ( 3, 5.6000000000000000000), 1E-13);

        assertEquals(0.3925671795844415, SpecialFunctions.besselJ( 4, 5.6000000000000000000), 1E-13);

        assertEquals(0.3310313267401661, SpecialFunctions.besselJ( 5, 5.6000000000000000000), 1E-13);

        assertEquals(-0.3110277443039424, SpecialFunctions.besselJ( 1, 5.8000000000000000000), 1E-13);

        assertEquals(-0.1989535138865205, SpecialFunctions.besselJ( 2, 5.8000000000000000000), 1E-13);

        assertEquals(0.1738184243822042, SpecialFunctions.besselJ( 3, 5.8000000000000000000), 1E-13);

        assertEquals(0.3787656770405248, SpecialFunctions.besselJ( 4, 5.8000000000000000000), 1E-13);

        assertEquals(0.3486169922254162, SpecialFunctions.besselJ( 5, 5.8000000000000000000), 1E-13);

        assertEquals(-0.2766838581275656, SpecialFunctions.besselJ( 1, 6.0000000000000000000), 1E-13);

        assertEquals(-0.2428732099601855, SpecialFunctions.besselJ( 2, 6.0000000000000000000), 1E-13);

        assertEquals(0.1147683848207753, SpecialFunctions.besselJ( 3, 6.0000000000000000000), 1E-13);

        assertEquals(0.3576415947809608, SpecialFunctions.besselJ( 4, 6.0000000000000000000), 1E-13);

        assertEquals(0.3620870748871724, SpecialFunctions.besselJ( 5, 6.0000000000000000000), 1E-13);

    }


    @Test
    public void testFactorial() {
        assertEquals(1, SpecialFunctions.factorial(0));
        assertEquals(1, SpecialFunctions.factorial(1));
        assertEquals(120, SpecialFunctions.factorial(5));
        assertEquals(2432902008176640000L, SpecialFunctions.factorial(20));

    }

    @Test
    public void testFactorialException() {
        // this asserts that an exception will be thrown for -1!
        thrown.expect(IllegalArgumentException.class);
        SpecialFunctions.factorial(-1);
    }

    @Test
    public void testLnFactorial() {
        assertEquals(0, SpecialFunctions.lnFactorial(1), EPSILON);
        assertEquals(10.60460290274525, SpecialFunctions.lnFactorial(8), EPSILON);
        assertEquals(42.335616460753485, SpecialFunctions.lnFactorial(20), EPSILON);
        assertEquals(363.73937555556349, SpecialFunctions.lnFactorial(100), EPSILON);
        assertEquals(82108.92783681436, SpecialFunctions.lnFactorial(10000), 1e-9);
    }


}
