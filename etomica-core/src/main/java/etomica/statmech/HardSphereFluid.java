/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.statmech;

/**
 * Computes EOS for HS fluid.  Uses exponential of a polynomial with low order
 * terms take from virial coefficients and high order terms fit to simulation
 * data.  The fitted terms take a different form that allows them to more
 * efficiently describe the high-density behavior.
 * 
 * @author Andrew Schultz
 */
public class HardSphereFluid {
    
    public final static double[] NexpVirial = {2.094395102393195,0.5483113556160761,-0.04333676720888957,0.05751876349502738,0.07502830097299396,
        -0.014474727265714415,-0.003183349327605778,0.017074836705042406,-0.0026047002163831975,-0.005498931252679095};
    public final static double[] NexpStretchFit = {6.960601433965988e-03,-2.263564675011192e-06,5.464112392968496e-10,-3.004093460638223e-14};
    public final static double a = 7123.0875;

    /**
     * Computes Z using the following equation
     * Z = exp([Ni * rho^(i-1)] + rho^11 [Mi x^i]) where x = a^rho
     *
     * @param rho
     * @return Z
     */

    public static double zFluid(double rho) {
        double expSum = 0;
        double rhoi = rho;
        for (int i=0; i<NexpVirial.length; i++) {
            expSum += NexpVirial[i]*rhoi;
            rhoi *= rho;
        }
        double rhon = rhoi;
        double x = Math.pow(a, rho);
        double xi = 1;
        for (int i=0; i<NexpStretchFit.length; i++) {
            expSum += rhon*NexpStretchFit[i]*xi;
            xi *= x;
        }
        return Math.exp(expSum);
    }
    
    /**
     * Returns the difference between absolute free energies between the given
     * densities for solid (FCC) hard spheres.  The result is computed using
     * numeric integration with n points.
     * @param rho Density
     * @return Absolute free energy difference
     */
    public static double Afluid(double rho){
        return Afluid(rho, 10000);
    }

    /**
     * Returns the difference between absolute free energies between the given
     * densities for solid (FCC) hard spheres.  The result is computed using
     * numeric integration with n points.
     * @param rho Density
     * @param n number of points
     * @return Absolute free energy difference
     */
    public static double Afluid(double rho, int n) {
        // use trapezoid rule to integrate the correction
        double h = rho/n;
        double sum = 0.0;
        for (int k=1; k<n; k++) {
            sum += (zFluid(k*h)-1)/(k*h);
        }
        sum += 0.5*(zFluid(rho)-1)/rho;
        return Math.log(rho)-1 + h*sum;
    }
}
