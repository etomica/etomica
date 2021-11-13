/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.potential;

import etomica.atom.IAtomList;
import etomica.space.Space;
import etomica.units.Kelvin;
import etomica.util.Constants;

/**
 * Simplified pair potential for Helium.  The potential has the form of an
 * exponential-6 model with an additional r^-10 term.
 *
 * @author Andrew Schultz
 */
public class P2HeSimplified extends Potential2SoftSpherical {

    public static Potential2Soft makeTruncated(Space space, TruncationFactory tf) {
        return tf.make(new P2HeSimplified(space));
    }

    public P2HeSimplified(Space space) {
        super(space);
    }

    /**
     * The energy u.
     */
    public double u(double r2) {

    	if (r2 < sigmaHC2) {
    		return Double.POSITIVE_INFINITY;
    	}

    	double r = Math.sqrt(r2);
        double r4 = r2*r2;
        double r6 = r4*r2;
        double r8 = r4*r4;
        double r10 = r6*r4;
        return A0*Math.exp(-A1*r)-A2/r6-A3/(useC10 ? r10 : r8);
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        double r = Math.sqrt(r2);
        double r4 = r2*r2;
        double r6 = r4*r2;
        double r8 = r4*r4;
        double r10 = r6*r4;
        return -A0*A1*Math.exp(-A1*r)*r + 6*A2/(r6) + A3*(useC10? 10/r10 : 8/r8);
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
        double r = Math.sqrt(r2);
        double r4 = r2*r2;
        double r6 = r4*r2;
        double r8 = r4*r4;
        double r10 = r6*r4;
        return A0*A1*A1*Math.exp(-A1*r)*r2 - 42*A2/(r6) - A3*(useC10? 110/r10 : 72/r8);
    }

    /**
     * Integral used for corrections to potential truncation.
     */
    public double integral(double rC) {
        double A = 4 * Math.PI;
        double sc = 1 / rC;
        double sc2 = sc * sc;
        double sc3 = sc * sc * sc;
        double sc6 = sc3 * sc3;
        double sc8 = sc6 * sc2;
        double sc10 = sc8 * sc2;
        double a = useC10 ? 10 : 8;
        // ignore exp term
        return -A * (A3 * (useC10 ? sc10 : sc8) / (a - 3) + A2 * sc6 / (6. - 3)) / sc3;
    }

    public P2HeQFH makeQFH(double temperature) {
        return this.new P2HeQFH(temperature);
    }

    /**
     * This inner class can calculates the Feynman-Hibbs semiclassical
     * approximation for the potential.  Results should be the same as using
     * P2EffectiveFeynmanHibbs, but should be faster, because u, du and d2u
     * are not called.  Much of the work is duplicated between those methods,
     * and they are all computed at the same time in this class.
     * 
     * This class is only ~10% slower than the classical potential.
     *
     * @author Andrew Schultz
     */
    public class P2HeQFH implements Potential2Spherical {

        protected final double temperature;
        protected final double mass = 4.002602;
        protected double fac;

        public P2HeQFH(double temperature) {
            this.temperature = temperature;
            double hbar = Constants.PLANCK_H/(2*Math.PI);
            fac = hbar*hbar/(24*mass/2)/temperature;
        }

        public double energy(IAtomList atoms) {
            dr.Ev1Mv2(atoms.get(1).getPosition(),atoms.get(0).getPosition());
            boundary.nearestImage(dr);
            return u(dr.squared());
        }

        public double getRange() {
            return P2HeSimplified.this.getRange();
        }

        public int nBody() {
            return 2;
        }

        public double u(double r2) {
            if (r2 < sigmaHC2) {
                return Double.POSITIVE_INFINITY;
            }

            double r = Math.sqrt(r2);

            double A1r = A1*r;
            double expA1r = Math.exp(-A1r);

            double s2 = 1.0/r2;
            double s4 = s2*s2;
            double s6 = s4*s2;
            double s8 = s4*s4;
            double s10 = s6*s4;

            double u = A0*expA1r*(1+fac*(-2*A1/r + A1*A1));
            if (!Double.isInfinite(u)) {
                // repulsion is not infinite.  add attractive dispersion
                u += A2*(-1 + fac*(2*6 - 42)*s2)*s6 + A3*(useC10 ? (-1 + fac*(2*10 - 110)*s2)*s10 : (-1 + fac*(2*8 - 72)*s2)*s8);
            }
            // if the classical potential is repulsive, the semiclassical potential
            // should be more repulsive.  In nearly all cases, it is, but this often
            // fails for very short separations.  Just enforce it here.
            // this is never needed for the potential here.  it would happen within our core
            //if (uc > 0 && u < uc) return uc;
            return u;
        }
    }
    
    public P2HeTI makeTI(double temperature) {
        return this.new P2HeTI(temperature);
    }
    public class P2HeTI implements Potential2Spherical {

        protected final double temperature;
        protected final double mass = 4.002602;
        protected double fac;

        public P2HeTI(double temperature) {
            this.temperature = temperature;
            double hbar = Constants.PLANCK_H/(2*Math.PI);
            fac = hbar*hbar/(24*mass/2)/(temperature*temperature);
        }

        public double energy(IAtomList atoms) {
            dr.Ev1Mv2(atoms.get(1).getPosition(),atoms.get(0).getPosition());
            boundary.nearestImage(dr);
            return u(dr.squared());
        }

        public double getRange() {
            return P2HeSimplified.this.getRange();
        }

        public int nBody() {
            return 2;
        }

        public double u(double r2) {
            if (r2 < sigmaHC2) {
                return Double.POSITIVE_INFINITY;
            }

            double r = Math.sqrt(r2);

            double A1r = A1*r;
            double expA1r = Math.exp(-A1r);

            double s2 = 1.0/r2;
            double s4 = s2*s2;
            double s6 = s4*s2;
            double s8 = s4*s4;
            double s10 = s6*s4;
          
            double uc = A0*expA1r - A2*s6  - A3*(useC10 ? s10 : s8);
            double duc = -A0*A1*r*expA1r + A2*6*s6 + A3*(useC10 ? 10*s10 : 8*s8);
            double usc = uc + fac/r2*(duc*duc);
            return usc;
        }
    }

    /**
     * Sets the ith parameter value
     * u = A0 * exp(-A1*r) - A2/r^6 - A3/r^X
     * 
     * where X = 10 if useC10 is true, X=8 otherwise
     */
    public void setA(int i, double a) {
        switch (i) {
            case 0:
                A0 = a;
                break;
            case 1:
                A1 = a;
                break;
            case 2:
                A2 = a;
                break;
            case 3: 
                A3 = a;
                break;
            default:
                throw new RuntimeException("oops");
        }
    }
    
    /**
     * Returns the ith parameter value
     * u = A0 * exp(-A1*r) - A2/r^6 - A3/r^X
     * 
     * where X = 10 if useC10 is true, X=8 otherwise
     */
    public double getA(int i) {
        switch (i) {
            case 0:
                return A0;
            case 1:
                return A1;
            case 2:
                return A2;
            case 3:
                return A3;
            default:
                throw new RuntimeException("oops");
        }
    }
    
    private static final long serialVersionUID = 1L;
    public static final boolean useC10 = false;
    protected double A0 = Kelvin.UNIT.toSim(7.09384366e6);
    protected double A1 = 4.47853246;
    protected double A2 = Kelvin.UNIT.toSim(8.70770913e3);
    protected double A3 = Kelvin.UNIT.toSim(6.14848434e4);
    protected final double sigmaHC2 = 1.6*1.6;
}
