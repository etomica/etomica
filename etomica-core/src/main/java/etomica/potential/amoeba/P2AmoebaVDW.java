/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential.amoeba;
import etomica.potential.IPotential2;
import etomica.potential.P2LennardJones;
import etomica.potential.TruncationFactory;
import etomica.space.Space;
import etomica.units.*;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;

/**
 * Ameoba potential for van der Waals interactions
 * Spherically symmetric potential of the form u(r) = epsilon*a^7*((1+gamma)/(s^7 + gamma) - 2)
 * with a = (1 + delta)/(s + delta)
 * and s = r/sigma
 *
 * Although u(sigma) = -epsilon, this is not the minimum.
 * Because of gamma and delta, minimum is not quite at sigma, nor is the well depth quite equal to epsilon
 *
 * @author David Kofke
 */
public class P2AmoebaVDW implements IPotential2 {

    public static IPotential2 makeTruncated(double sigma, double epsilon, double delta, double gamma, TruncationFactory tf) {
        return tf.make(new P2AmoebaVDW(sigma, epsilon, delta, gamma));
    }

    public P2AmoebaVDW() {
        this(1.0, 1.0, 0, 0);
    }

    public P2AmoebaVDW(double sigma, double epsilon, double delta, double gamma) {
        super();
        setSigma(sigma);
        setEpsilon(epsilon);
        setDelta(delta);
        setGamma(gamma);
    }

    /**
     * The energy u.
     */
    public double u(double r2) {
        double r = Math.sqrt(r2);
        double a = (1 + delta)/(r/sigma + delta);
        double a2 = a*a;
        double a7 = a2*a2*a2*a;
        double s = r/sigma;
        double s2 = s*s;
        double s6 = s2*s2*s2;
        Unit kcalpmole = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT), Mole.UNIT);
        double u = epsilon*a7*((1+gamma)/(s6*s + gamma) - 2);
        System.out.println(r+" "+" "+r2+" "+u+" "+kcalpmole.fromSim(u));
        return epsilon*a7*((1+gamma)/(s6*s + gamma) - 2);
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        double r = Math.sqrt(r2);
        double aa = (r/sigma + delta);
        double a = (1 + delta)/aa;
        double a2 = a*a;
        double a7 = a2*a2*a2*a;
        double s = r/sigma;
        double s2 = s*s;
        double s6 = s2*s2*s2;
        return -7*epsilon*(a7 / aa * s + a7 * (1 + gamma)/(s6*s + gamma) * s6*s);
    }

    public void u012add(double r2, double[] u012) {
        double r = Math.sqrt(r2);
        double s = r/sigma;
        double aa = (s + delta);
        double a = (1 + delta)/aa;
        double a2 = a*a;
        double a7 = a2*a2*a2*a;
        double s2 = s*s;
        double s6 = s2*s2*s2;
        double bb = s6*s + gamma;
        double c = a7*((1+gamma)/bb - 2);
        u012[0] += epsilon*c;
        u012[1] -= 7*epsilon*(c / aa * s + a7 * (1 + gamma)/(bb*bb) * s6*s);
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
        throw new RuntimeException("nope");
    }
            
    /**
     * Integral used for corrections to potential truncation.
     */
    public double integral(Space space, double rC) {
        throw new RuntimeException("nope");
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
    }

    public Dimension getSigmaDimension() {return Length.DIMENSION;}
    
    /**
     * Accessor method for Lennard-Jones energy parameter
     */
    public double getEpsilon() {return epsilon;}
    /**
     * Mutator method for Lennard-Jones energy parameter
     */
    public void setEpsilon(double eps) {
        epsilon = eps;
    }
    public Dimension getEpsilonDimension() {return Energy.DIMENSION;}

    public void setDelta(double d) {
        delta = d;
    }

    public void setGamma(double g) {
        gamma = g;
    }
   
    private double sigma, epsilon;
    private double delta, gamma;

    public static void main(String[] args) {
        double s = 1.5;
        P2AmoebaVDW p2 = new P2AmoebaVDW(s, 2, 0.07, 0.12);
        P2LennardJones p2LJ = new P2LennardJones(s/Math.pow(2, 1.0/6.0), 2);
        for (int i=10; i<200; i++) {
            double r = s*i*0.01;
            double[] u012 = new double[3];
            p2.u012add(r*r, u012);
            double u = u012[0];
            u012[0] = 0;
            p2LJ.u012add(r*r, u012);
            System.out.println(r+" "+u+" "+u012[0]);
        }
    }
}
