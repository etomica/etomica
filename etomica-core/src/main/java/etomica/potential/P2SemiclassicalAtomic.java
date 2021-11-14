/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.Constants;

import java.util.HashMap;
import java.util.Map;

/**
 * Effective semiclassical molecular potential using the approach of Takahashi and Imada
 * <p>
 * http://dx.doi.org/10.1143/JPSJ.53.3765
 * <p>
 * as described by Schenter
 * <p>
 * http://dx.doi.org/10.1063/1.1505441
 *
 * @author Andrew Schultz
 */
public class P2SemiclassicalAtomic implements Potential2Soft {

    protected final Potential2Soft p2Classy;
    protected final Map<AtomType, AtomInfo> agents;
    protected final Space space;
    protected double temperature, fac;

    public P2SemiclassicalAtomic(Space space, Potential2Soft p2Classy, double temperature) {
        this.space = space;
        this.p2Classy = p2Classy;
        agents = new HashMap<>();
        setTemperature(temperature);
    }

    public void setAtomInfo(AtomType species, AtomInfo moleculeInfo) {
        agents.put(species, moleculeInfo);
    }

    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
        double hbar = Constants.PLANCK_H / (2 * Math.PI);
        fac = hbar * hbar / (24 * temperature * temperature);
    }

    public double getRange() {
        return p2Classy.getRange();
    }

    public double u(Vector dr12, IAtom atom1, IAtom atom2) {
        double uC = p2Classy.u(dr12, atom1, atom2);
        if (uC / temperature > 100) return Double.POSITIVE_INFINITY;
        Vector f1 = space.makeVector();
        Vector f2 = space.makeVector();
        Vector t1 = space.makeVector();
        Vector t2 = space.makeVector();
        p2Classy.uduTorque(dr12, atom1, atom2, f1, f2, t1, t2);
        Vector[][] gradAndTorque = new Vector[][]{{f1,f2},{t1,t2}};
        double sum = 0;
        for (int i = 0; i < 2; i++) {
            IAtom iMol = i == 0 ? atom1 : atom2;
            double mi = iMol.getType().getMass();
            sum += gradAndTorque[0][i].squared() / mi;
            if (iMol instanceof IAtomOriented) {
                AtomInfo atomInfo = agents.get(iMol.getType());
                if (atomInfo == null) {
                    throw new RuntimeException("You must provide AtomInfo for oriented atoms");
                }
                Vector[] momentAndAxes = atomInfo.getMomentAndAxes((IAtomOriented) iMol);
                Vector moment = momentAndAxes[0];
                for (int j = 0; j < 3; j++) {
                    if (moment.getX(j) < 1e-10) continue;
                    Vector axis = momentAndAxes[j + 1];
                    double torque = gradAndTorque[1][i].dot(axis);
                    sum += torque * torque / moment.getX(j);
                }
            }
        }
        double uFull = uC + fac * sum;
//        System.out.println(uC+" "+uFull);
        return uFull;
    }

    public double energy(IAtomList atoms) {
        IAtom atom1 = atoms.get(0);
        IAtom atom2 = atoms.get(1);
        Vector dr12 = space.makeVector();
        dr12.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
        return u(dr12, atoms.get(0), atoms.get(1));
    }

    public interface AtomInfo {

        /**
         * Returns the moment of inertia as a vector containing Ix, Iy, Iz and also the principle axes.  The 0 element
         * is the moment of inertia vector while elements 1, 2 and 3 are the principle axes.
         */
        Vector[] getMomentAndAxes(IAtomOriented molecule);
    }
}
