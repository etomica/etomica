/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeChildAtomAction;
import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.box.Box;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculePair;
import etomica.molecule.MoleculeSource;
import etomica.molecule.MoleculeSourceRandomMolecule;
import etomica.potential.IPotentialMolecular;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Standard Monte Carlo molecule-displacement trial move.  Two molecules are moved at a
 * time in such a way that the geometric center of the system is not changed.
 *
 * @author Nancy Cribbin
 * 	       Tai Boon Tan
 */
public class MCMoveMoleculeCoupledFasterer extends MCMoveBoxStep {

    protected final PotentialCompute potentialCompute;
    protected final MoleculeChildAtomAction moveMoleculeAction;
    protected final Vector groupTransVect;
    protected IMolecule molecule0, molecule1;
    protected MoleculeSource moleculeSource;
    protected double uOld, uNew;
    protected final IRandom random;
    protected final AtomIteratorArrayListSimple affectedMoleculeIterator;
    protected final AtomArrayList affectedMoleculeList;
    protected final AtomActionTranslateBy singleAction;
    protected final MoleculePair pair;
    protected IPotentialMolecular pairPotential;
    protected boolean doExcludeNonNeighbors, doIncludePair;
    protected IAtom[] atoms;

    public MCMoveMoleculeCoupledFasterer(PotentialCompute potentialCompute, IRandom nRandom,
                                         Space _space){
        super();
        this.potentialCompute = potentialCompute;
        this.random = nRandom;
        moleculeSource = new MoleculeSourceRandomMolecule();
        ((MoleculeSourceRandomMolecule)moleculeSource).setRandomNumberGenerator(random);

        affectedMoleculeList = new AtomArrayList();
        affectedMoleculeIterator = new AtomIteratorArrayListSimple(affectedMoleculeList);
        
        singleAction = new AtomActionTranslateBy(_space);
        groupTransVect = singleAction.getTranslationVector();
        
        moveMoleculeAction = new MoleculeChildAtomAction(singleAction);
        
        pair = new MoleculePair();
        
        perParticleFrequency = true;
        //energyMeter.setIncludeLrc(false);
    }

    public void setBox(Box newBox) {
        super.setBox(newBox);
        moleculeSource.setBox(newBox);
    }
    
    public void setPotential(IPotentialMolecular newPotential){
        pairPotential = newPotential;
    }
    
    public AtomIterator affectedAtoms() {
        affectedMoleculeList.clear();
        affectedMoleculeList.addAll(molecule0.getChildList());
        affectedMoleculeList.addAll(molecule1.getChildList());
        return affectedMoleculeIterator;
    }

    public double energyChange() {return uNew - uOld;}

    public boolean doTrial() {
//        System.out.println("doTrial MCMoveMoleculeCoupled called");
        
        molecule0 = moleculeSource.getMolecule();
        molecule1 = moleculeSource.getMolecule();
        if(molecule0==null || molecule1==null || molecule0==molecule1) return false;
        
        //make sure we don't double count the molecule0-molecule1 interaction
        pair.mol0 = molecule0;
        pair.mol1 = molecule1;
        atoms = new IAtom[molecule0.getChildList().size()+molecule1.getChildList().size()];
        int idx = 0;
        for (IAtom a : molecule0.getChildList()) atoms[idx++] = a;
        for (IAtom a : molecule1.getChildList()) atoms[idx++] = a;
        uOld = potentialCompute.computeManyAtomsOld(atoms);
        if (uOld == Double.POSITIVE_INFINITY || Double.isNaN(uOld)) {
            throw new RuntimeException("infinite or NaN uOld "+uOld);
        }

        if(uOld > 1e10){
            throw new ConfigurationOverlapException(box);
        }
        
        groupTransVect.setRandomCube(random);
        groupTransVect.TE(stepSize);
        moveMoleculeAction.actionPerformed(molecule0);
        groupTransVect.TE(-1.0);
        moveMoleculeAction.actionPerformed(molecule1);
        uNew = potentialCompute.computeManyAtoms(atoms);

        return true;
    }

    public double getChi(double temperature) {
        return Math.exp(-(uNew - uOld) / temperature);
    }

    public void acceptNotify() {
        potentialCompute.processAtomU(1);
        // put it back, then compute old contributions to energy
        moveMoleculeAction.actionPerformed(molecule0);
        groupTransVect.TE(-1.0);
        moveMoleculeAction.actionPerformed(molecule1);

        potentialCompute.computeManyAtoms(atoms);

        moveMoleculeAction.actionPerformed(molecule0);
        groupTransVect.TE(-1.0);
        moveMoleculeAction.actionPerformed(molecule1);

        potentialCompute.processAtomU(-1);
    }

    public void rejectNotify() {
        moveMoleculeAction.actionPerformed(molecule0);
        groupTransVect.TE(-1.0);
        moveMoleculeAction.actionPerformed(molecule1);
    }
    
    /**
     * Configures the move to not explicitly calculate the potential for atom
     * pairs that are not neighbors (as determined by the potentialMaster).
     * 
     * Setting this has an effect only if the potentialMaster is an instance of
     * PotentialMasterList
     */
    public void setDoExcludeNonNeighbors(boolean newDoExcludeNonNeighbors) {
        doExcludeNonNeighbors = newDoExcludeNonNeighbors;
    }

    /**
     * Returns true if the move does not explicitly calculate the potential for
     * atom pairs that are not neighbors (as determined by the potentialMaster).
     */
    public boolean getDoExcludeNonNeighbors() {
        return doExcludeNonNeighbors;
    }
    
}
