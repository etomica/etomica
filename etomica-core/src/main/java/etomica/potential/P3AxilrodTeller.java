/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;

import java.util.Map;

/**
 * Axilrod-Teller potential.  The potential is atomic.  Ionization energy and
 * polarizability are required as input for each atom type. 
 * 
 * @author Andrew Schultz
 */
public class P3AxilrodTeller implements IPotentialAtomic, Potential3Soft {

    protected final Map<AtomType, MyAgent> paramsManager;
    protected final Space space;
    protected final Vector dr1, dr2;
    protected Boundary boundary;
    
    public P3AxilrodTeller(Space space, Map<AtomType, MyAgent> paramsManager) {
        this.space = space;
        this.paramsManager = paramsManager;
        dr1 = space.makeVector();
        dr2 = space.makeVector();
    }

    @Override
    public double u(Vector dr12, Vector dr13, Vector dr23, IAtom atom1, IAtom atom2, IAtom atom3) {

        double RAB2 = dr12.squared();
        double RAC2 = dr13.squared();
        double RBC2 = dr23.squared();
        if (RAB2*RAC2*RBC2 == 0) return 0;

        double RAB = Math.sqrt(RAB2);
        double RAC = Math.sqrt(RAC2);
        double RBC = Math.sqrt(RBC2);

        double costhetaA = (RAB2 + RAC2 - RBC2)/(2*RAC*RAB);
        double costhetaB = (RAB2 + RBC2 - RAC2)/(2*RAB*RBC);
        double costhetaC = (RAC2 + RBC2 - RAB2)/(2*RAC*RBC);

        IAtom[] atoms = new IAtom[]{atom1, atom2, atom3};
        double ep = 1;
        double es = 0;
        double[] eps = new double[3];
        double ap = 1;
        double cp = 3 * costhetaA * costhetaB * costhetaC;
        for (int i=0; i<3; i++) {
            MyAgent ag = paramsManager.get(atoms[i].getType());
            eps[(i+1)%3] += ag.E;
            eps[(i+2)%3] += ag.E;
            ep *= ag.E;
            es += ag.E;

            ap *= ag.alpha;
        }
        double e123 = ep*es;
        for (int i=0; i<3; i++) {
            e123 /= eps[i];
        }

        return 1.5*e123*ap*(cp+1)/(RAB*RAB2*RAC*RAC2*RBC*RBC2);
    }

    public double energy(IAtomList atoms) {
        double ep = 1;
        double es = 0;
        double[] eps = new double[3];
        double ap = 1;
        double cp = 3;
        double rp = 1;
        for (int i=0; i<3; i++) {
            MyAgent ag = paramsManager.get(atoms.get(i).getType());
            eps[(i+1)%3] += ag.E;
            eps[(i+2)%3] += ag.E;
            ep *= ag.E;
            es += ag.E;

            ap *= ag.alpha;
            dr1.Ev1Mv2(atoms.get(i).getPosition(),atoms.get((i+1)%3).getPosition());
            boundary.nearestImage(dr1);
            double r2 = dr1.squared();
            rp *= r2*Math.sqrt(r2);
            dr2.Ev1Mv2(atoms.get((i+2)%3).getPosition(),atoms.get((i+1)%3).getPosition());
            boundary.nearestImage(dr2);
            double cos = dr1.dot(dr2)/Math.sqrt(r2*dr2.squared());
            cp *= cos;
        }
        double e123 = ep*es;
        for (int i=0; i<3; i++) {
            e123 /= eps[i];
        }
        
        return 1.5*e123*ap*(cp+1)/rp;
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    // CO2: alpha=2.913, alpha(omega)=2.639, I=13.7eV
    // H2O: alpha=1.444, I=12.6eV
    public static class MyAgent {
        public final double alpha, E;
        /**
         * @params alpha polarizability
         * @params E ionization energy
         */
        public MyAgent(double alpha, double E) {
            this.alpha = alpha;
            this.E = E;
        }
    }
}
