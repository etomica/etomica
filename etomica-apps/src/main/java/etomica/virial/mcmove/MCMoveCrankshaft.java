/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.virial.mcmove;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.potential.compute.NeighborManagerCell;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;

/**
 * Crankshaft Monte Carlo move for Mayer sampling of star polymers: a core atom with numArms linear
 * arms of armLength atoms each, with the atoms ordered core first, then arm 0 from the core
 * outward, then arm 1, and so on (as made by ConformationStarPolymerGraft with SpeciesBuilder).
 * <p>
 * Each trial picks a molecule, an arm and two atoms a and b on that arm, b being 2 to maxSpan bonds
 * farther from the core than a (a can be the core itself), and rotates the atoms strictly between
 * a and b about the axis through a and b by a random angle in [-stepSize, stepSize].  Since a and b
 * are on the axis, all bond lengths are kept, and nothing outside the segment moves.  Unlike a pivot
 * (MCMoveClusterAngle), which swings the whole outer part of the arm, this changes the conformation
 * locally, so it can take larger steps in crowded regions such as near the core.
 * <p>
 * The position of a along the arm is chosen uniformly in [aMin, aMax] (0 is the core, 1 the grafted
 * monomer, armLength the arm end), limited so that b stays on the arm, and the span uniformly in
 * [2, maxSpan], likewise limited.  The reverse move picks the same a and b with the opposite angle,
 * with the same probability, so the acceptance is (wNew/wOld) exp(-(uNew - uOld)/T), as for
 * MCMoveClusterAngle, where w is the sampling cluster's value (unless skipW) and u the energy from
 * the PotentialCompute (the intramolecular energy for the cluster drivers).
 * <p>
 * Not set up for flexible-correction alternate molecules (constraint maps): each trial moves atoms
 * of one molecule only.
 */
public class MCMoveCrankshaft extends MCMoveBoxStep {

    protected final PotentialCompute potential;
    protected final IRandom random;
    protected final Space space;
    protected final int numArms, armLength, aMin, aMax, maxSpan;
    protected final RotationTensor3D rotationTensor;
    protected final Vector axis, center;
    protected Vector[] oldPositions;
    protected int[] movedAtoms;
    protected int nMoved, iMolecule;
    protected double uOld, uNew, wOld, wNew;
    protected NeighborManagerCell cellManager;
    public boolean skipW = false;

    /**
     * @param aMin    smallest position along the arm for atom a (0 = core)
     * @param aMax    largest position along the arm for atom a
     * @param maxSpan largest number of bonds between a and b (at least 2)
     */
    public MCMoveCrankshaft(PotentialCompute potentialCompute, Space space, IRandom random, double stepSize,
                            int numArms, int armLength, int aMin, int aMax, int maxSpan) {
        super();
        if (aMin < 0 || aMax < aMin || aMin > armLength - 2 || maxSpan < 2) {
            throw new IllegalArgumentException("need 0 <= aMin <= aMax, aMin <= armLength-2 and maxSpan >= 2");
        }
        this.potential = potentialCompute;
        this.space = space;
        this.random = random;
        this.stepSize = stepSize;
        this.numArms = numArms;
        this.armLength = armLength;
        this.aMin = aMin;
        this.aMax = Math.min(aMax, armLength - 2);
        this.maxSpan = maxSpan;
        setStepSizeMax(Math.PI);
        rotationTensor = new RotationTensor3D();
        axis = space.makeVector();
        center = space.makeVector();
    }

    public void setBox(Box p) {
        super.setBox(p);
        oldPositions = space.makeVectorArray(armLength);
        movedAtoms = new int[armLength];
    }

    public void setCellManager(NeighborManagerCell cellManager) {
        this.cellManager = cellManager;
    }

    /** index within the molecule of the atom at position p of the given arm (p = 0 is the core) */
    protected int atomIndex(int arm, int p) {
        return p == 0 ? 0 : 1 + armLength * arm + p - 1;
    }

    public boolean doTrial() {
        uOld = potential.computeAll(false);
        wOld = wNew = 1;
        if (!skipW) wOld = ((BoxCluster) box).getSampleCluster().value((BoxCluster) box);

        iMolecule = random.nextInt(box.getMoleculeList().size());
        IAtomList atoms = box.getMoleculeList().get(iMolecule).getChildList();
        int arm = random.nextInt(numArms);
        int pa = aMin + random.nextInt(aMax - aMin + 1);
        int span = 2 + random.nextInt(Math.min(maxSpan, armLength - pa) - 1);
        int pb = pa + span;

        center.E(atoms.get(atomIndex(arm, pa)).getPosition());
        axis.Ev1Mv2(atoms.get(atomIndex(arm, pb)).getPosition(), center);
        axis.normalize();
        double dt = stepSize * (2 * random.nextDouble() - 1);
        rotationTensor.setRotationAxis(axis, dt);

        nMoved = 0;
        for (int p = pa + 1; p < pb; p++) {
            int i = atomIndex(arm, p);
            Vector r = atoms.get(i).getPosition();
            oldPositions[nMoved].E(r);
            movedAtoms[nMoved++] = i;
            r.ME(center);
            rotationTensor.transform(r);
            r.PE(center);
        }

        if (cellManager != null) cellManager.assignCellAll();
        ((BoxCluster) box).trialNotify();
        if (!skipW) wNew = ((BoxCluster) box).getSampleCluster().value((BoxCluster) box);
        uNew = potential.computeAll(false);
        return true;
    }

    public double getChi(double temperature) {
        return (wOld == 0 ? 1 : wNew / wOld) * Math.exp(-(uNew - uOld) / temperature);
    }

    public double energyChange() {
        return 0;
    }

    public void acceptNotify() {
        ((BoxCluster) box).acceptNotify();
    }

    public void rejectNotify() {
        IAtomList atoms = box.getMoleculeList().get(iMolecule).getChildList();
        for (int k = 0; k < nMoved; k++) {
            atoms.get(movedAtoms[k]).getPosition().E(oldPositions[k]);
        }
        if (cellManager != null) cellManager.assignCellAll();
        ((BoxCluster) box).rejectNotify();
    }
}
