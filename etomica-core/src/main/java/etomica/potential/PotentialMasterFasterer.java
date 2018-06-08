/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.species.ISpecies;

public class PotentialMasterFasterer {
    protected final Potential2Soft[][] pairPotentials;
    protected final Box box;
    protected final Vector dr;
    protected double[] uAtom, duAtom;
    protected int[] uAtomsChanged2;
    protected int nChanged;
    protected Vector zero;
    protected double virialTot;
    protected Vector[] forces;

    public PotentialMasterFasterer(Simulation sim, Box box) {
        ISpecies species = sim.getSpecies(sim.getSpeciesCount() - 1);
        int lastTypeIndex = species.getAtomType(species.getAtomTypeCount() - 1).getIndex();
        pairPotentials = new Potential2Soft[lastTypeIndex + 1][lastTypeIndex + 1];
        this.box = box;
        dr = sim.getSpace().makeVector();
        uAtom = new double[box.getLeafList().size()];
        uAtomsChanged2 = new int[box.getLeafList().size()];
        nChanged = 0;
        zero = box.getSpace().makeVector();
    }

    public double getLastVirial() {
        return virialTot;
    }

    public double getOldEnergy() {
        double uTot = 0;
        for (double iuAtom : uAtom) {
            uTot += iuAtom;
        }
        return uTot;
    }

    public void setPairPotential(AtomType atomType1, AtomType atomType2, Potential2SoftSpherical p12) {
        pairPotentials[atomType1.getIndex()][atomType2.getIndex()] = p12;
    }

    public void updateAtom(IAtom atom) {

    }

    protected double handleComputeAll(boolean doForces, int iAtom, int jAtom, Vector ri, Vector rj, Vector jbo, Potential2Soft pij) {
        dr.Ev1Mv2(rj, ri);
        dr.PE(jbo);
        double[] u = {0};
        double[] du = {0};
        double r2 = dr.squared();
        pij.udu(r2, u, du);
        double uij = u[0];
//        double uij = pij.u(dr.squared());
        if (uij == 0) return 0;
//        System.out.println(iAtom+" "+jAtom+" "+uij);
        uAtom[iAtom] += 0.5 * uij;
        uAtom[jAtom] += 0.5 * uij;
        double duij = du[0];
        virialTot += duij;
        if (doForces) {
            duij /= r2;
            dr.TE(duij / r2);
            forces[iAtom].PE(dr);
            forces[iAtom].ME(dr);
        }
        return uij;
    }

    public double computeAll(boolean doForces) {
        virialTot = 0;
        int numAtoms = box.getLeafList().size();
        if (numAtoms > uAtom.length) {
            uAtom = new double[numAtoms];
            duAtom = new double[numAtoms];
            uAtomsChanged2 = new int[numAtoms];
        } else {
            for (int i = 0; i < numAtoms; i++) uAtom[i] = 0;
        }
        IAtomList atoms = box.getLeafList();
        double u = 0;
        Boundary boundary = box.getBoundary();
        for (int i = 0; i < atoms.size(); i++) {
            IAtom iAtom = atoms.get(i);
            Vector ri = iAtom.getPosition();
            int iType = iAtom.getType().getIndex();
            Potential2Soft[] ip = pairPotentials[iType];
            for (int j = i + 1; j < atoms.size(); j++) {
                IAtom jAtom = atoms.get(j);
                int jType = jAtom.getType().getIndex();
                Potential2Soft pij = ip[jType];
                if (pij == null) continue;
                dr.Ev1Mv2(ri, jAtom.getPosition());
                boundary.nearestImage(dr);
                u += handleComputeAll(doForces, i, j, zero, dr, zero, pij);
            }
        }
        return u;
    }

    public double computeOneOld(IAtom iAtom) {
        return uAtom[iAtom.getLeafIndex()] * 2;
    }

    protected double handleComputeOne(Potential2Soft pij, Vector ri, Vector rj, Vector jbo, int iAtom, int jAtom) {
        dr.Ev1Mv2(rj, ri);
        dr.PE(jbo);
        double uij = pij.u(dr.squared());
        if (uij == 0) return 0;
        duAtom[0] += 0.5 * uij;
        duAtom[nChanged] += 0.5 * uij;
        uAtomsChanged2[nChanged] = jAtom;
        nChanged++;
        return uij;
    }

    public double computeOne(IAtom iAtom) {
        int i = iAtom.getLeafIndex();
        Vector ri = iAtom.getPosition();
        int iType = iAtom.getType().getIndex();
        uAtomsChanged2[0] = i;
        nChanged = 1;
        duAtom[0] = 0;
        Potential2Soft[] ip = pairPotentials[iType];
        IAtomList atoms = box.getLeafList();
        double u = 0;
        Boundary boundary = box.getBoundary();
        int s = atoms.size();
        for (int j = 0; j < s; j++) {
            if (i == j) continue;
            IAtom jAtom = atoms.get(j);
            int jType = jAtom.getType().getIndex();
            Potential2Soft pij = ip[jType];
            if (pij == null) continue;
            dr.Ev1Mv2(ri, jAtom.getPosition());
            boundary.nearestImage(dr);
            u += handleComputeOne(pij, zero, dr, zero, i, j);
        }
        return u;
    }

    public void processAtomU(double fac) {
        for (int j = 0; j < nChanged; j++) {
            int jj = uAtomsChanged2[j];
            uAtom[jj] += fac * duAtom[j];
            duAtom[j] = 0;
        }
        nChanged = 0;
    }

    public void resetAtomDU() {
        for (int j = 0; j < nChanged; j++) {
            duAtom[j] = 0;
        }
        nChanged = 0;
    }
}
