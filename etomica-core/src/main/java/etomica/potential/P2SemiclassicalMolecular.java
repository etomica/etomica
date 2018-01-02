/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.Constants;

import java.util.HashMap;
import java.util.Map;

/**
 * Effective semiclassical molecular potential using the approach of
 * Takahashi and Imada
 * 
 * http://dx.doi.org/10.1143/JPSJ.53.3765
 * 
 * as described by Schenter
 * 
 * http://dx.doi.org/10.1063/1.1505441
 * 
 * @author Andrew Schultz
 */
public class P2SemiclassicalMolecular implements IPotentialMolecular {

    protected final IPotentialMolecularTorque p2Classy;
    protected double temperature, fac;
    private final Map<ISpecies, MoleculeInfo> agents;
    protected final Space space;
    
    public P2SemiclassicalMolecular(Space space, IPotentialMolecularTorque p2Classy) {
        this.space = space;
        this.p2Classy = p2Classy;
        if (p2Classy.nBody() != 2) throw new RuntimeException("I would really rather have a 2-body potential");
        this.agents = new HashMap<>();
    }
    
    public void setMoleculeInfo(ISpecies species, MoleculeInfo moleculeInfo) {
        agents.put(species, moleculeInfo);
    }
    
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
        double hbar = Constants.PLANCK_H/(2*Math.PI);
        fac = hbar*hbar/(24*temperature*temperature);
    }
    
    public double getRange() {
        return p2Classy.getRange();
    }

    public void setBox(Box box) {
        p2Classy.setBox(box);
    }

    public int nBody() {
        return 2;
    }

    public double energy(IMoleculeList molecules) {
        double uC = p2Classy.energy(molecules);
        Vector[][] gradAndTorque = p2Classy.gradientAndTorque(molecules);
        double sum = 0;
        for (int i=0; i<2; i++) {
            IMolecule iMol = molecules.getMolecule(i);
            MoleculeInfo molInfo = agents.get(iMol.getType());
            if (molInfo == null) {
                molInfo = new MoleculeInfoBrute(space);
                agents.put(iMol.getType(), molInfo);
            }
            double mi = molInfo.getMass(iMol);
            sum -= gradAndTorque[i][0].squared()/mi;
            Vector[] momentAndAxes = molInfo.getMomentAndAxes(iMol);
            Vector moment = momentAndAxes[0];
            for (int j=0; j<3; j++) {
                if (moment.getX(j) < 1e-10) continue;
                Vector axis = momentAndAxes[j+1];
                double torque = gradAndTorque[i][1].dot(axis);
                sum += torque*torque/moment.getX(j);
            }
        }
        double uFull = uC + fac*sum;
        return uFull;
    }

    public interface MoleculeInfo {
        /**
         * Returns the mass of the molecule.
         */
        public double getMass(IMolecule molecule);
        
        /**
         * Returns the moment of inertia as a vector containing Ix, Iy, Iz
         * and also the principle axes.  The 0 element is the moment of inertia
         * vector while elements 1, 2 and 3 are the principle axes.
         */
        public Vector[] getMomentAndAxes(IMolecule molecule);
    }
    
    public static class MoleculeInfoBrute implements MoleculeInfo {
        protected final Vector cm, rj;
        protected final Tensor id, moment, momentj, rjrj;
        protected final Vector[] rv;
        
        public MoleculeInfoBrute(Space space) {
            moment = space.makeTensor();
            momentj = space.makeTensor();
            rj = space.makeVector();
            rjrj = space.makeTensor();
            cm = space.makeVector();
            id = space.makeTensor();
            id.setComponent(0,0,1);
            id.setComponent(1,1,1);
            id.setComponent(2,2,1);
            rv = new Vector[4];
            for (int i=0; i<4; i++) {
                rv[i] = space.makeVector();
            }
        }

        public double getMass(IMolecule molecule) {
            double m = 0;
            moment.E(0);
            IAtomList atoms = molecule.getChildList();
            for (int j=0; j<atoms.getAtomCount(); j++) {
                IAtom a = atoms.getAtom(j);
                double mj = a.getType().getMass();
                m += mj;
            }
            return m;
        }
        
        public Vector[] getMomentAndAxes(IMolecule molecule) {
            double m = 0;
            moment.E(0);
            IAtomList atoms = molecule.getChildList();
            for (int j=0; j<atoms.getAtomCount(); j++) {
                IAtom a = atoms.getAtom(j);
                double mj = a.getType().getMass();
                cm.PEa1Tv1(mj, a.getPosition());
                m += mj;
            }
            cm.TE(1.0/m);
            for (int j=0; j<atoms.getAtomCount(); j++) {
                IAtom a = atoms.getAtom(j);
                double mj = a.getType().getMass();
                rj.Ev1Mv2(a.getPosition(), cm);
                momentj.E(id);
                momentj.TE(rj.squared());
                momentj.PEv1v2(rj, rj);
                momentj.TE(mj);
                moment.PE(momentj);
            }
            Matrix matrix = new Matrix(moment.toArray(), 3);
            EigenvalueDecomposition ed = new EigenvalueDecomposition(matrix);
            double[] evals = ed.getRealEigenvalues();
            double[][] evecs = ed.getV().getArray();
            rv[0].E(evals);
            for (int i=0; i<3; i++) {
                rv[i+1].E(evecs[i]);
            }
            return rv;
        }
    }
}
