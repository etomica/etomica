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
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorMC;
import etomica.molecule.*;
import etomica.space.Space;
import etomica.util.random.IRandom;

/**
 * Performs a trial that results in the exchange of a molecule from one box to another.
 * Primary use is as an elementary move in a Gibbs ensemble simulation
 *
 * @author David Kofke
 */

public class MCMoveMoleculeExchange extends MCMove {

    protected Box box1;
    protected Box box2;
    protected final IntegratorBox integrator1, integrator2;
    private final AtomIteratorArrayListSimple affectedAtomIterator = new AtomIteratorArrayListSimple();
    private final MoleculeActionTranslateTo moleculeTranslator;
    private final IRandom random;
    private MoleculeSource moleculeSource;
    protected RandomPositionSource positionSource;

    private transient IMolecule dMolecule, iMolecule;
    private transient Box iBox, dBox;
    private IntegratorBox iIntegrator;
    private transient double uOld;
    private transient double uNew = Double.NaN;


    public MCMoveMoleculeExchange(IRandom random,
                                  Space space,
                                  IntegratorBox integrator1,
                                  IntegratorBox integrator2) {
        super();
        this.random = random;
        moleculeTranslator = new MoleculeActionTranslateTo(space);
        setAtomPositionDefinition(new MoleculePositionCOM(space));
        this.integrator1 = integrator1;
        this.integrator2 = integrator2;
        box1 = integrator1.getBox();
        box2 = integrator2.getBox();
        moleculeSource = new MoleculeSourceRandomMolecule();
        ((MoleculeSourceRandomMolecule) moleculeSource).setRandomNumberGenerator(random);
        positionSource = new RandomPositionSourceRectangular(space, random);
    }

    public boolean doTrial() {
        IntegratorBox dIntegrator;
        if (random.nextInt(2) == 0) {
            iBox = box1;
            dBox = box2;
            iIntegrator = integrator1;
            dIntegrator = integrator2;
        } else {
            iBox = box2;
            dBox = box1;
            iIntegrator = integrator2;
            dIntegrator = integrator1;
        }
        if (dBox.getMoleculeList().size() == 0) { //no molecules to delete; trial is over
            uNew = uOld = 0.0;
            return false;
        }

        moleculeSource.setBox(dBox);
        dMolecule = moleculeSource.getMolecule();  //select random molecule to delete

        uOld = dIntegrator.getPotentialCompute().computeOneOldMolecule(dMolecule);
        // copy molecule and insert copy; remove from dBox in acceptNotify
//        dBox.removeMolecule(molecule);
        iMolecule = dMolecule.getType().makeMolecule();
        iMolecule.copyCoordinatesFrom(dMolecule);

        positionSource.setBox(iBox);
        moleculeTranslator.setDestination(positionSource.randomPosition());         //place at random in insertion box
        moleculeTranslator.actionPerformed(iMolecule);
        iBox.addMolecule(iMolecule);
        uNew = Double.NaN;
        return true;
    }//end of doTrial

    /**
     * Sets a new RandomPositionSource for this move to use.  By default, a
     * position source is used which assumes rectangular boundaries.
     */
    public void setPositionSource(RandomPositionSource newPositionSource) {
        positionSource = newPositionSource;
    }

    /**
     * Returns the RandomPositionSource used by this move.
     */
    public RandomPositionSource getPositionSource() {
        return positionSource;
    }

    /**
     * Sets the AtomSource this class uses to pick molecules to delete.
     */
    public void setMoleculeSource(MoleculeSource newMoleculeSource) {
        moleculeSource = newMoleculeSource;
    }

    /**
     * Returns the AtomSource this class uses to pick molecules to delete.
     */
    public MoleculeSource getMoleculeSource() {
        return moleculeSource;
    }

    public double getChi(double temperature) {
        uNew = iIntegrator.getPotentialCompute().computeOneMolecule(iMolecule);
        double B = -(uNew - uOld);
        // assume both integrators have the same temperature
        //note that dSpecies.nMolecules has been decremented
        //and iSpecies.nMolecules has been incremented
        return Math.exp(B / temperature) * (dBox.getNMolecules(dMolecule.getType()) + 1) / dBox.getBoundary().volume()
                * iBox.getBoundary().volume() / iBox.getNMolecules(iMolecule.getType());
    }

    public void acceptNotify() {
        IntegratorBox dIntegrator = integrator2;
        if (dIntegrator.getBox() == iBox) {
            dIntegrator = integrator1;
        }
        iIntegrator.getPotentialCompute().processAtomU(1);

        dIntegrator.getPotentialCompute().computeOneMolecule(dMolecule);
        dIntegrator.getPotentialCompute().processAtomU(-1);
        // accepted deletion - remove from box
        dBox.removeMolecule(dMolecule);

        if (iIntegrator instanceof IntegratorMC) {
            ((IntegratorMC) iIntegrator).notifyEnergyChange(uNew);
        } else {
            //XXX grossly inefficient
            iIntegrator.reset();
        }
        if (dIntegrator instanceof IntegratorMC) {
            ((IntegratorMC) dIntegrator).notifyEnergyChange(-uOld);
        } else {
            //XXX grossly inefficient
            dIntegrator.reset();
        }
    }

    public void rejectNotify() {
        iBox.removeMolecule(iMolecule);
    }

    public final AtomIterator affectedAtoms(Box box) {
        if (this.box1 != box && this.box2 != box) return AtomIteratorNull.INSTANCE;
        affectedAtomIterator.setList(iMolecule.getChildList());
        return affectedAtomIterator;
    }

    public double energyChange(Box box) {
        if (box == iBox) return uNew;
        else if (box == dBox) return -uOld;
        else return 0.0;
    }


    /**
     * @return Returns the atomPositionDefinition.
     */
    public IMoleculePositionDefinition getAtomPositionDefinition() {
        return moleculeTranslator.getAtomPositionDefinition();
    }

    /**
     * @param atomPositionDefinition The atomPositionDefinition to set.
     */
    public void setAtomPositionDefinition(
            IMoleculePositionDefinition atomPositionDefinition) {
        moleculeTranslator.setAtomPositionDefinition(atomPositionDefinition);
    }
}//end of MCMoveMoleculeExchange
