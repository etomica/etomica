/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.action.MoleculeActionTranslateTo;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.atom.iterator.AtomIteratorNull;
import etomica.box.Box;
import etomica.box.RandomPositionSource;
import etomica.box.RandomPositionSourceRectangular;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;
import etomica.molecule.MoleculeSetSinglet;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.util.Debug;
import etomica.util.random.IRandom;

/**
 * Elementary Monte Carlo move in which a molecule of a specified species is
 * inserted into or removed from a box.
 *
 * @author David Kofke
 */
public class MCMoveInsertDelete extends MCMoveBox {

    protected final PotentialCompute potentialMaster;
    protected final AtomIteratorArrayListSimple affectedAtomIterator = new AtomIteratorArrayListSimple();
    protected final MoleculeArrayList reservoir;
    protected final MoleculeActionTranslateTo atomTranslator;
    //chemical potential
    protected double mu;
    protected ISpecies species;
    protected IMolecule testMolecule;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected boolean insert;
    protected IMoleculeList moleculeList;
    protected IRandom random;
    protected RandomPositionSource positionSource;

    /**
     * Constructs using default isPotentialHard = false
     *
     * @param potentialMaster used for calculation of energies
     * @param random          random number generator used by the simulation
     * @param space           governing space for the simulation
     */
    public MCMoveInsertDelete(PotentialCompute potentialMaster, IRandom random, Space space) {
        super();
        this.potentialMaster = potentialMaster;
        setMu(0.0);
        atomTranslator = new MoleculeActionTranslateTo(space);
        reservoir = new MoleculeArrayList();
        this.random = random;
        perParticleFrequency = true;
        positionSource = new RandomPositionSourceRectangular(space, random);
    }

    public ISpecies getSpecies() {
        return species;
    }

    public void setSpecies(ISpecies s) {
        species = s;
        if (box != null) {
            moleculeList = box.getMoleculeList(species);
        }
    }

    public void setBox(Box p) {
        super.setBox(p);
        if (species != null) {
            moleculeList = box.getMoleculeList(species);
        }
        positionSource.setBox(box);
    }

    /**
     * Returns the RandomPositionSource used by this move.
     */
    public RandomPositionSource getPositionSource() {
        return positionSource;
    }

    /**
     * Sets a new RandomPositionSource for this move to use.  By default, a
     * position source is used which assumes rectangular boundaries.
     */
    public void setPositionSource(RandomPositionSource newPositionSource) {
        positionSource = newPositionSource;
        if (box != null) {
            positionSource.setBox(box);
        }
    }

    /**
     * Chooses and performs with equal probability an elementary molecule insertion
     * or deletion.
     */
    public boolean doTrial() {
        insert = (random.nextInt(2) == 0);
        if (insert) {
            uOld = 0.0;

            if (!reservoir.isEmpty()) testMolecule = reservoir.remove(reservoir.size() - 1);
            else testMolecule = species.makeMolecule();

            atomTranslator.setDestination(positionSource.randomPosition());
            atomTranslator.actionPerformed(testMolecule);
            box.addMolecule(testMolecule);
        } else {//delete
            if (box.getNMolecules(species) == 0) {
                testMolecule = null;
                return false;
            }
            testMolecule = moleculeList.get(random.nextInt(moleculeList.size()));
            //delete molecule only upon accepting trial

            uOld = potentialMaster.computeOneOldMolecule(testMolecule);
        }
        uNew = Double.NaN;
        return true;
    }//end of doTrial

    public double getChi(double temperature) {//note that moleculeCount() gives the number of molecules after the trial is attempted
        int numMolecules = box.getNMolecules(species);
        double a = box.getBoundary().volume() / numMolecules;

        if (insert) {
            uNew = potentialMaster.computeOneMolecule(testMolecule);
        } else {
            uNew = 0;
        }
        double b = uOld - uNew;
        if (insert) b += mu;
        else b -= mu;

        return (insert ? a : 1.0 / a) * Math.exp(b / temperature);
    }

    public void acceptNotify() {
        if (insert && Debug.ON && Debug.DEBUG_NOW && Debug.allAtoms(new MoleculeSetSinglet(testMolecule)))
            System.out.println("accepted insertion of " + testMolecule);
        if (insert) {
            potentialMaster.processAtomU(1);
        } else {
            potentialMaster.computeOneMolecule(testMolecule);
            potentialMaster.processAtomU(-1);
            // accepted deletion - remove from box and add to reservoir
            box.removeMolecule(testMolecule);
            reservoir.add(testMolecule);
        }
    }

    public void rejectNotify() {
        if (insert) {
            // rejected insertion - remove from box and return to reservoir
            box.removeMolecule(testMolecule);
            reservoir.add(testMolecule);
            // test molecule is no longer in the simulation and should not be 
            // returned by affectedAtoms
            testMolecule = null;
        }
    }

    public double energyChange() {
        return uNew - uOld;
    }

    /**
     * Returns an iterator giving molecule that is being added or deleted
     * in the current or most recent trial.
     */
    public AtomIterator affectedAtoms() {
        if (testMolecule == null) return AtomIteratorNull.INSTANCE;
        affectedAtomIterator.setList(testMolecule.getChildList());
        return affectedAtomIterator;
    }

    public boolean lastMoveInsert() {
        return insert;
    }

    /**
     * Accessor method for the chemical potential of th insertion/deletion species.
     */
    public final double getMu() {
        return mu;
    }

    /**
     * Mutator method for the chemical potential of the insertion/deletion species.
     */
    public void setMu(double mu) {
        this.mu = mu;
    }

    /**
     * Indicates that chemical potential has dimensions of energy.
     */
    public final Dimension getMuDimension() {
        return Energy.DIMENSION;
    }

}//end of MCMoveInsertDelete
