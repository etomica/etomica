/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.starpolymer;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.potential.PotentialCalculationEnergySum;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;

/**
 * An MC Move for cluster simulations that "wiggles" a chain molecule.  If the
 * first or last atom in the chain is chosen, it is moved to a new position
 * with the same bond length as before, but perturbed by some angle from its
 * original position.  If an Atom in the middle of the chain is chosen, a
 * crankshaft move is performed that maintains its distances with its
 * neighbors.  If a middle Atom has a bond angle too close to 180 degrees
 * (such that rotation does nothing) the Atom is not moved at all.
 * In each doTrial, wiggle moves are attempted on all molecules in the box.
 *
 * @author Andrew Schultz
 */
public class MCMoveBondLength extends MCMoveMolecule {

    public MCMoveBondLength(Simulation sim, PotentialMaster potentialMaster, int armLength, Space _space) {
        this(potentialMaster, sim.getRandom(), 1.0, armLength, _space);
    }

    public MCMoveBondLength(PotentialMaster potentialMaster,
                            IRandom random, double stepSize, int armLength, Space _space) {
        super(potentialMaster, random, _space, stepSize, Double.POSITIVE_INFINITY);
        this.space = _space;
        this.armLength = armLength;
        translateVector = _space.makeVector();
        energyMeter = new MeterPotentialEnergy(potential);
    }

    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
    }

    public void setSpecies(ISpecies newSpecies) {
        species = newSpecies;
    }

    //note that total energy is calculated
    public boolean doTrial() {
        molecule = moleculeSource.getMolecule();
        energyMeter.setTarget(molecule);
        uOld = energyMeter.getDataAsScalar();
        if (uOld > 1e8) {
            PotentialCalculationEnergySum.debug = true;
            energyMeter.getDataAsScalar();
            throw new RuntimeException("molecule " + molecule + " in box " + box + " has an overlap");
        }
        IAtomList childList = molecule.getChildList();
        int numChildren = childList.size();

        startAtom = random.nextInt(numChildren - 1) + 1;

        IAtom fixedAtom = childList.get(startAtom - 1);
        if (startAtom % armLength == 1 || armLength == 1) {
//            System.out.println("star atom " + startAtom);
            fixedAtom = childList.get(0);
        }
        r0 = fixedAtom.getPosition();
//            System.out.println(selectedAtoms[i]+" "+j+" before "+selectedAtoms[i].coord.position());
        dr = stepSize * 2 * (random.nextDouble() - 0.5);
        rejectme = false;
        doTransform();
        if (rejectme) return false;

        return true;
    }


    protected void doTransform() {
        IAtomList childList = molecule.getChildList();
        int stopAtom = (((startAtom - 1) / armLength) + 1) * armLength;
        Vector startPos = childList.get(startAtom).getPosition();
        translateVector.Ev1Mv2(startPos, r0);
        bl = Math.sqrt(translateVector.squared());
        if (dr < -bl) {
            rejectme = true;
            return;
        }
        translateVector.TE(dr / bl);
        for (int iChild = startAtom; iChild <= stopAtom; iChild++) {
            IAtom a = childList.get(iChild);
            Vector r = a.getPosition();
            r.PE(translateVector);
        }
    }

    public void acceptNotify() {
//        if(r0.isZero()) System.out.println("accepted! ");
//        System.out.println("accepting bond length move");
    }

    public void rejectNotify() {
        dr = -dr;
//        System.out.println("rejecting bond length move");
        doTransform();
    }

    public double getChi(double temperature) {
        uNew = energyMeter.getDataAsScalar();
        double blNew = bl + dr;
        double ratio = blNew / bl;
        return Math.exp(-(uNew - uOld) / temperature) * ratio * ratio;
    }

    protected final MeterPotentialEnergy energyMeter;
    protected Vector r0;
    protected double dr;
    protected int startAtom;
    protected final Space space;
    protected ISpecies species;
    protected final int armLength;
    protected final Vector translateVector;
    protected double bl;
    protected boolean rejectme;
}
