/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.nbr.list;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.potential.BondingInfo;
import etomica.potential.Potential2Soft;
import etomica.potential.PotentialMasterFasterer;
import etomica.potential.compute.PotentialCallback;
import etomica.space.Vector;
import etomica.species.SpeciesManager;

public class PotentialMasterListFasterer extends PotentialMasterFasterer implements IntegratorListener {

    private final NeighborListManagerFasterer nbrManager;

    public PotentialMasterListFasterer(SpeciesManager sm, Box box, int cellRange, double nbrRange, BondingInfo bondingInfo) {
        this(sm, box, cellRange, nbrRange, bondingInfo, new Potential2Soft[sm.getAtomTypeCount()][sm.getAtomTypeCount()]);
    }

    public PotentialMasterListFasterer(SpeciesManager sm, Box box, int cellRange, double nbrRange, BondingInfo bondingInfo, Potential2Soft[][] pairPotentials) {
        super(sm, box, bondingInfo, pairPotentials);
        this.nbrManager = new NeighborListManagerFasterer(sm, box, cellRange, nbrRange, bondingInfo);
        this.nbrManager.setPairPotentials(this.pairPotentials);
    }

    public NeighborListManagerFasterer getNeighborManager() {
        return nbrManager;
    }

    public void init() {
        this.nbrManager.init();
    }

    public void setNeighborRange(double newRange) {
        this.nbrManager.setNeighborRange(newRange);
    }

    // this determines how big the cells are.  PMCell just returns potential range
    public double getNeighborRange() {
        return this.nbrManager.getNeighborRange();
    }

    public void setDoDownNbrs(boolean doDown) {
        this.nbrManager.setDoDownNeighbors(doDown);
    }

    public double computeAll(boolean doForces, PotentialCallback pc) {
        zeroArrays(doForces);

        double uTot = 0;
        IAtomList atoms = box.getLeafList();
        int numAtoms = atoms.size();

        for (int i = 0; i < numAtoms; i++) {
            IAtom iAtom = atoms.get(i);
            Vector ri = iAtom.getPosition();
            int iType = iAtom.getType().getIndex();
            Potential2Soft[] iPotentials = pairPotentials[iType];
            int iNumNbrs = nbrManager.numAtomNbrsUp[i];
            int[] iNbrs = nbrManager.nbrs[i];
            Vector[] iNbrBoxOffsets = nbrManager.nbrBoxOffsets[i];
            for (int j = 0; j < iNumNbrs; j++) {
                int jj = iNbrs[j];
                IAtom jAtom = atoms.get(jj);
                int jType = jAtom.getType().getIndex();
                Potential2Soft pij = iPotentials[jType];
                Vector rj = jAtom.getPosition();
                Vector jbo = iNbrBoxOffsets[j];
                uTot += handleComputeAll(doForces, i, jj, ri, rj, jbo, pij, pc, false);
            }
        }

        double[] uCorrection = new double[1];
        double[] duCorrection = new double[1];
        this.computeAllTruncationCorrection(uCorrection, duCorrection);
        uTot += uCorrection[0];
        virialTot += duCorrection[0];
        energyTot = uTot;
        return uTot;
    }

    @Override
    public IntegratorListener makeIntegratorListener() {
        return this;
    }

    @Override
    public void integratorInitialized(IntegratorEvent e) {
        init();
    }

    @Override
    public void integratorStepStarted(IntegratorEvent e) {
        this.nbrManager.checkUpdateNbrs();
    }

    @Override
    public void integratorStepFinished(IntegratorEvent e) {
    }

    protected double computeOneInternal(IAtom atom) {
        int iType = atom.getType().getIndex();
        int i = atom.getLeafIndex();
        Potential2Soft[] ip = pairPotentials[iType];
        double u = 0;
        long t1 = System.nanoTime();
        Vector ri = atom.getPosition();

        IAtomList atoms = box.getLeafList();
        int iNumNbrs = nbrManager.numAtomNbrsUp[i];
        int[] iNbrs = nbrManager.nbrs[i];
        Vector[] iNbrBoxOffsets = nbrManager.nbrBoxOffsets[i];
        for (int j = 0; j < iNumNbrs; j++) {
            int jj = iNbrs[j];
            IAtom jAtom = atoms.get(jj);
            int jType = jAtom.getType().getIndex();
            Potential2Soft pij = ip[jType];
            Vector rj = jAtom.getPosition();
            Vector jbo = iNbrBoxOffsets[j];
            boolean skipIntra = bondingInfo.skipBondedPair(isPureAtoms, atom, jAtom);
            u += handleComputeOne(pij, ri, rj, jbo, i, jj, skipIntra);
        }

        iNumNbrs = nbrManager.numAtomNbrsDn[i];
        for (int j = 0; j < iNumNbrs; j++) {
            int jj = iNbrs[iNbrs.length - 1 - j];
            IAtom jAtom = atoms.get(jj);
            int jType = jAtom.getType().getIndex();
            Potential2Soft pij = ip[jType];
            Vector rj = jAtom.getPosition();
            Vector jbo = iNbrBoxOffsets[iNbrs.length - 1 - j];
            boolean skipIntra = bondingInfo.skipBondedPair(isPureAtoms, atom, jAtom);
            u += handleComputeOne(pij, rj, ri, jbo, i, jj, skipIntra);
        }

        tMC += System.nanoTime() - t1;
        return u;
    }

}
