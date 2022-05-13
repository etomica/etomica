/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.exception.MethodNotImplementedException;
import etomica.space.Vector;

/**
 * Methods for properties obtained for a soft, differentiable pair potential.
 *
 * @author David Kofke
 */
public interface IPotential3 {

    /**
     * Returns the range over which the potential applies.  IAtoms with a
     * greater separation do not interact.
     */
    default double getRange() {return Double.POSITIVE_INFINITY;}

    /**
     * The pair energy u(r^2) with no truncation applied.
     *  with r2xy the square of the distance between the particle x and y.
     */
    default double u(double r212, double r213, double r223) {
        throw new MethodNotImplementedException();
    }

    /**
     * Returns the energy between IAtoms atom1 and atom2 separated by vector
     * dr12 (which goes from atom1 to atom2; dr12 = r2 - r1).  PBC should not
     * be applied to dr12; PBC has already been accounted for and dr12 may even
     * exceed the Box dimensions when a lattice sum is used.  Likewise, the
     * IAtoms' actual positions (from getPosition()) ought to be ignored.
     */
    default double u(Vector dr12, Vector dr13, Vector dr23, IAtom atom1, IAtom atom2, IAtom atom3, double[] virial) {
        return u(dr12.squared(), dr13.squared(), dr23.squared());
    }

    default double udu(Vector dr12, Vector dr13, Vector dr23, IAtom atom1, IAtom atom2, IAtom atom3, double[] virial, Vector f1, Vector f2, Vector f3) {
        throw new MethodNotImplementedException();
    }

    /**
     * This method handles the vector math to compute the force contributions for each atom
     */
    default void forceHelper(Vector drAB, Vector drAC, Vector drBC, double rAB, double rAC, double rBC,
                             double costhetaA, double costhetaB, double costhetaC,
                             double rABdudrAB, double rACdudrAC, double rBCdudrBC,
                             double cosAdudcosA, double cosBdudcosB, double cosCdudcosC,
                             Vector fA, Vector fB, Vector fC) {
        double rAB2 = rAB*rAB;
        double rAC2 = rAC*rAC;
        double rBC2 = rBC*rBC;
                
        Vector tmp = Vector.d(drAB.getD());
        tmp.Ea1Tv1(-rABdudrAB/rAB2, drAB);
        tmp.PEa1Tv1(-rACdudrAC/rAC2, drAC);

        Vector tmp2 = Vector.d(drAB.getD());  // dcosthetaBdrA
        tmp2.Ea1Tv1(1/rAB/rBC, drBC);
        tmp2.PEa1Tv1(costhetaB/rAB2, drAB);
        Vector dcosthetaBdrA = Vector.d(drAB.getD());
        dcosthetaBdrA.E(tmp2);

        tmp.PEa1Tv1(cosBdudcosB/costhetaB, tmp2);

        // dcosthetaCdrA
        tmp2.Ea1Tv1(1/rAC/rBC, drBC);
        tmp2.PEa1Tv1(-costhetaC/rAC2, drAC);
        tmp.PEa1Tv1(-cosCdudcosC/costhetaC, tmp2);

        // dcothetaAdrA = -(dcosthetaAdrB + dcosthetaAdrC)
        tmp2.Ea1Tv1(1/rAB/rAC, drAC); // dcosthetaAdrB
        tmp2.PEa1Tv1(-costhetaA/rAB2, drAB);
        Vector dcosthetaAdrB = Vector.d(drAB.getD());
        dcosthetaAdrB.E(tmp2);
        tmp2.PEa1Tv1(1/rAC/rAB, drAB); // dcosthetaAdrC
        tmp2.PEa1Tv1(-costhetaA/rAC2, drAC);
        tmp.PEa1Tv1(-cosAdudcosA/costhetaA, tmp2);
        // our tmp is the gradient for A, so subtract to get force
        fA.ME(tmp);
        fC.PE(tmp);

        // fB
        tmp.Ea1Tv1(rABdudrAB/rAB2, drAB);
        tmp.PEa1Tv1(-rBCdudrBC/rBC2, drBC);

        // dcosthetaAdrB
        tmp.PEa1Tv1(cosAdudcosA/costhetaA, dcosthetaAdrB);

        // dcosthetaCdrB
        tmp2.Ea1Tv1(1/rBC/rAC, drAC);
        tmp2.PEa1Tv1(-costhetaC/rBC2, drBC);
        tmp.PEa1Tv1(-cosCdudcosC/costhetaC, tmp2);

        // dcothetaBdrB = -(dcosthetaBdrA + dcosthetaBdrC)
        tmp2.E(dcosthetaBdrA); // dcosthetaBdrA
        tmp2.PEa1Tv1(-1/rBC/rAB, drAB); // dcosthetaBdrC
        tmp2.PEa1Tv1(-costhetaB/rBC2, drBC);
        tmp.PEa1Tv1(-cosBdudcosB/costhetaB, tmp2);
        // our tmp is the gradient for B, so subtract to get force
        fB.ME(tmp);
        fC.PE(tmp);
    }
}