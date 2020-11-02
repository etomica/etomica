/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.nbr.list;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.nbr.cell.PotentialMasterCellFasterer;
import etomica.potential.BondingInfo;
import etomica.potential.Potential2Soft;
import etomica.potential.PotentialMasterFasterer;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.util.Debug;

public class PotentialMasterListFasterer extends PotentialMasterFasterer implements IntegratorListener {

    private final NeighborListManagerFasterer nbrManager;

    public PotentialMasterListFasterer(Simulation sim, Box box, int cellRange, double nbrRange, BondingInfo bondingInfo) {
        super(sim, box, bondingInfo);
        this.nbrManager = new NeighborListManagerFasterer(sim, box, cellRange, nbrRange, bondingInfo);
        this.nbrManager.setPairPotentials(this.pairPotentials);
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

    public double computeAll(boolean doForces) {
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
                uTot += handleComputeAll(doForces, i, jj, ri, rj, jbo, pij);
            }
        }

        double[] uCorrection = new double[1];
        double[] duCorrection = new double[1];
        this.computeAllTruncationCorrection(uCorrection, duCorrection);
        uTot += uCorrection[0];
        virialTot += duCorrection[0];

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
}
