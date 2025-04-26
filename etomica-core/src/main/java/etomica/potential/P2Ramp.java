/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;
import etomica.space.Space;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;

/**
 * Lennard-Jones interatomic potential.
 * Spherically symmetric potential of the form u(r) = 4*epsilon*[(sigma/r)^12 - (sigma/r)^6]
 * where epsilon describes the strength of the pair interaction,
 * and sigma is the atom size parameter.
 *
 * @author David Kofke
 */
public class P2Ramp implements IPotential2 {

    public P2Ramp(double sigma, double epsilon) {
        super();
        setSigma(sigma);
        setEpsilon(epsilon);
    }

    /**
     * The energy u.
     */
    public double u(double r2) {
        if (r2 > sigma2) return 0;
        double r = Math.sqrt(r2);
        return epsilon*(sigma - r);
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        if (r2 > sigma2) return 0;
        double r = Math.sqrt(r2);
        return -epsilon*r;
    }

    /**
     * The second derivative of the pair energy, times the square of the
     * separation:  r^2 d^2u/dr^2.
     */
    public double d2u(double r2) {
        return 0;
    }

    public void u012add(double r2, double[] u012) {
        if (r2 > sigma2) return;
        double r = Math.sqrt(r2);
        u012[0] += epsilon*(sigma-r);
        u012[1] += -epsilon*r;
        u012[2] += 0;
    }

    /**
     * Accessor method for Lennard-Jones size parameter.
     */
    public double getSigma() {return sigma;}
    /**
     * Mutator method for Lennard-Jones size parameter.
     * Does not adjust potential cutoff for change in sigma.
     */
    public final void setSigma(double s) {
        sigma = s;
        sigma2 = s*s;
    }
    public Dimension getSigmaDimension() {return Length.DIMENSION;}
    
    /**
     * Accessor method for Lennard-Jones energy parameter
     */
    public double getEpsilon() {return epsilon;}
    /**
     * Mutator method for Lennard-Jones energy parameter
     */
    public final void setEpsilon(double eps) {
        epsilon = eps;
    }
    public Dimension getEpsilonDimension() {return Energy.DIMENSION;}

    @Override
    public double getRange() {
        return sigma;
    }

    private double sigma, sigma2;
    private double epsilon;
}
