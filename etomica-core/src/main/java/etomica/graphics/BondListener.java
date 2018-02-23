/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;
import etomica.atom.iterator.AtomsetIteratorDirectable;
import etomica.box.Box;
import etomica.chem.models.Model;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculeSetSinglet;
import etomica.molecule.iterator.MoleculeIteratorMolecule;
import etomica.potential.IPotential;
import etomica.potential.IteratorDirective.Direction;
import etomica.species.ISpecies;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * BondListener listens for Atoms being added to the Simulation, determines
 * what covalent bonds (if any) apply to the new Atom and informs a BondManager
 * of any new bonds.  When an Atom is removed from the Simulation, the
 * BondListener determines what (if any) bonds the atom participates in and
 * informs the BondManager that the bond is going away.
 * 
 * @author Andrew Schultz
 */
public class BondListener implements AtomLeafAgentManager.AgentSource<ArrayList>, Serializable {
    
    /**
     * Creates a new BondListener for the given Box, using the given
     * BondManager to actually create or remove bonds.
     */
    public BondListener(Box box, BondManager bondManager) {
        this.box = box;
        bondIteratorsHash = new HashMap<ISpecies,Model.PotentialAndIterator[]>();
        atomAgentManager = new AtomLeafAgentManager<ArrayList>(this, box);
        this.bondManager = bondManager;
        atomSetSinglet = new MoleculeSetSinglet();
    }

    /**
     * Infroms the BondListener that is should track bonds associated with all
     * molecules of the given Model.
     *
     * This method should not be called until the model has been added to the
     * simulation.
     */
    public void addModel(Model newModel) {
        ISpecies species = newModel.getSpecies();
        Model.PotentialAndIterator[] bondIterators = newModel.getPotentials();
        bondIteratorsHash.put(species, bondIterators);

        // now find bonds for atoms that already exist
        // the ArrayLists (agents) for the Atoms will already exist
        MoleculeIteratorMolecule moleculeIterator = new MoleculeIteratorMolecule(species);
        moleculeIterator.setBox(box);
        moleculeIterator.reset();
        for (IMolecule molecule = moleculeIterator.nextMolecule(); molecule != null;
             molecule = moleculeIterator.nextMolecule()) {
            // we have an molecule, now grab all of its bonds

            for (int i=0; i<bondIterators.length; i++) {
                IPotential bondedPotential = bondIterators[i].getPotential();
                AtomsetIteratorBasisDependent iterator = bondIterators[i].getIterator();
                if (iterator instanceof AtomsetIteratorDirectable) {
                    // these should all be directable, but perhaps not
                    ((AtomsetIteratorDirectable)iterator).setDirection(Direction.UP);
                }
                else {
                    System.err.println("iterator wasn't directable, strange things may happen");
                }
                atomSetSinglet.atom = molecule;
                iterator.setBasis(atomSetSinglet);
                iterator.setTarget(null);
                iterator.reset();
                for  (IAtomList bondedPair = iterator.next(); bondedPair != null;
                      bondedPair = iterator.next()) {
                    
                    Object bond = bondManager.makeBond(bondedPair, bondedPotential);

                    atomAgentManager.getAgent(bondedPair.get(0)).add(bond);
                }
            }
        }
    }

    /**
     * Informs the BondListener that bonds in molecules of the given model
     * should no longer be tracked.
     */
    public void removeModel(Model model) {
        ISpecies species = model.getSpecies();
        MoleculeIteratorMolecule moleculeIterator = new MoleculeIteratorMolecule(species);
        moleculeIterator.setBox(box);
        moleculeIterator.reset();
        for (IMolecule molecule = moleculeIterator.nextMolecule(); molecule != null;
             molecule = moleculeIterator.nextMolecule()) {
            IAtomList childList = molecule.getChildList();
            for (int iChild = 0; iChild < childList.size(); iChild++) {
                ArrayList list = atomAgentManager.getAgent(childList.get(iChild));
                for (int i=0; i<list.size(); i++) {
                    bondManager.releaseBond(list.get(i));
                }
            }
        }
            
            
        bondIteratorsHash.remove(species);
    }

    public Class getAgentClass() {
        return ArrayList.class;
    }
    
    public ArrayList makeAgent(IAtom newAtom, Box agentBox) {
        // we got a leaf atom in a mult-atom molecule
        ArrayList<Object> bondList = new ArrayList<Object>(); 
        Model.PotentialAndIterator[] bondIterators = bondIteratorsHash.
                get(newAtom.getType().getSpecies());
        
        if (bondIterators != null) {
            IMolecule molecule = newAtom.getParentGroup();
            
            for (int i=0; i<bondIterators.length; i++) {
                IPotential bondedPotential = bondIterators[i].getPotential();
                AtomsetIteratorBasisDependent iterator = bondIterators[i].getIterator();
                // We only want bonds where our atom of interest is the "up" atom.
                // We'll pick up bonds where this atom is the "down" atom when
                // makeAgent is called for that Atom.  We're dependent here on 
                // a molecule being added to the system only after it's 
                // completely formed.  Not depending on that would be "hard".
                if (iterator instanceof AtomsetIteratorDirectable) {
                    // these should all be directable, but perhaps not
                    ((AtomsetIteratorDirectable)iterator).setDirection(Direction.UP);
                }
                else {
                    System.err.println("iterator wasn't directable, strange things may happen");
                }
                atomSetSinglet.atom = molecule;
                iterator.setBasis(atomSetSinglet);
                iterator.setTarget(newAtom);
                iterator.reset();
                for (IAtomList bondedAtoms = iterator.next(); bondedAtoms != null;
                     bondedAtoms = iterator.next()) {
                    Object bond = bondManager.makeBond(bondedAtoms, bondedPotential);
                    bondList.add(bond);
                }
            }
        }
        return bondList;
    }
    
    public void releaseAgent(ArrayList agent, IAtom atom, Box agentBox) {
        // we only release a bond when the "up" atom from the bond goes away
        // so if only the "down" atom goes away, we would leave the bond in
        // (bad).  However, you're not allowed to mutate the model, so deleting
        // one leaf atom of a covalent bond and not the other would be illegal.
        ArrayList<Object> list = agent;
        for (int i=0; i<list.size(); i++) {
            bondManager.releaseBond(list.get(i));
        }
        list.clear();
    }
    
    private static final long serialVersionUID = 1L;
    protected final Box box;
    protected final AtomLeafAgentManager<ArrayList> atomAgentManager;
    protected final HashMap<ISpecies,Model.PotentialAndIterator[]> bondIteratorsHash;
    protected BondManager bondManager;
    protected final MoleculeSetSinglet atomSetSinglet;
}
