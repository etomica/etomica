package etomica.graphics;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;

import etomica.atom.AtomAgentManager;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.atom.SpeciesAgent;
import etomica.atom.iterator.AtomIteratorMolecule;
import etomica.atom.iterator.AtomIteratorTreeRoot;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;
import etomica.atom.iterator.AtomsetIteratorDirectable;
import etomica.atom.iterator.IteratorDirective.Direction;
import etomica.chem.models.Model;
import etomica.phase.Phase;
import etomica.potential.Potential;
import etomica.species.Species;

/**
 * BondListener listens for Atoms being added to the Simulation, determines
 * what covalent bonds (if any) apply to the new Atom and informs a BondManager
 * of any new bonds.  When an Atom is removed from the Simulation, the
 * BondListener determines what (if any) bonds the atom participates in and
 * informs the BondManager that the bond is going away.
 * 
 * @author Andrew Schultz
 */
public class BondListener implements AtomAgentManager.AgentSource, Serializable {

    /**
     * Creates a new BondListener for the given Phase, using the given
     * BondManager to actually create or remove bonds.  The BondListener
     * (and its AtomAgentManager) are consdidered to be "backend".
     */
    public BondListener(Phase phase, BondManager bondManager) {
        this(phase, bondManager, true);
    }
    
    /**
     * Creates a new BondListener for the given Phase, using the given
     * BondManager to actually create or remove bonds.  The BondListener
     * (and its AtomAgentManager, BondManager and bonds) are consdidered to be
     * "backend" (they will be serialized along with the simulation) if the
     * given isBackend parameter is true.  
     */
    public BondListener(Phase phase, BondManager bondManager, boolean isBackend) {
        this.phase = phase;
        bondIteratorsHash = new HashMap();
        atomAgentManager = new AtomAgentManager(this, phase, isBackend);
        this.bondManager = bondManager;
    }

    /**
     * Infroms the BondListener that is should track bonds associated with all
     * molecules of the given Model.
     *
     * This method should not be called until the model has been added to the
     * simulation.
     */
    public void addModel(Model newModel) {
        Species species = newModel.getSpecies();
        Model.PotentialAndIterator[] bondIterators = newModel.getPotentials();
        bondIteratorsHash.put(species, bondIterators);

        // now find bonds for atoms that already exist
        // the ArrayLists (agents) for the Atoms will already exist
        AtomIteratorMolecule moleculeIterator = new AtomIteratorMolecule(new Species[]{species});
        moleculeIterator.setPhase(phase);
        moleculeIterator.reset();
        while(moleculeIterator.hasNext()) {
            IAtom molecule = moleculeIterator.nextAtom();
            // we have an molecule, now grab all of its bonds

            for (int i=0; i<bondIterators.length; i++) {
                Potential bondedPotential = bondIterators[i].getPotential();
                AtomsetIteratorBasisDependent iterator = bondIterators[i].getIterator();
                if (iterator instanceof AtomsetIteratorDirectable) {
                    // these should all be directable, but perhaps not
                    ((AtomsetIteratorDirectable)iterator).setDirection(Direction.UP);
                }
                else {
                    System.err.println("iterator wasn't directable, strange things may happen");
                }
                iterator.setBasis(molecule);
                iterator.setTarget(null);
                iterator.reset();
                while (iterator.hasNext()) {
                    AtomPair bondedPair = (AtomPair)iterator.next();
                    
                    Object bond = bondManager.makeBond(bondedPair, bondedPotential);

                    ((ArrayList)atomAgentManager.getAgent(bondedPair.getAtom(0))).add(bond);
                }
            }
        }
    }

    /**
     * Informs the BondListener that bonds in molecules of the given model
     * should no longer be tracked.
     */
    public void removeModel(Model model) {
        Species species = model.getSpecies();
        AtomIteratorMolecule moleculeIterator = new AtomIteratorMolecule(new Species[]{species});
        moleculeIterator.setPhase(phase);
        moleculeIterator.reset();
        AtomIteratorTreeRoot leafIterator = new AtomIteratorTreeRoot();
        while(moleculeIterator.hasNext()) {
            IAtom molecule = moleculeIterator.nextAtom();
            leafIterator.setRootAtom(molecule);
            leafIterator.reset();
            while (leafIterator.hasNext()) {
                IAtom leafAtom = leafIterator.nextAtom();
                ArrayList list = (ArrayList)atomAgentManager.getAgent(leafAtom);
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
    
    public Object makeAgent(IAtom newAtom) {
        if (!(newAtom.getParentGroup() instanceof SpeciesAgent) && newAtom.isLeaf()) {
            // we got a leaf atom in a mult-atom molecule
            ArrayList bondList = new ArrayList(); 
            Model.PotentialAndIterator[] bondIterators = 
                (Model.PotentialAndIterator[])bondIteratorsHash.
                    get(newAtom.getType().getSpecies());
            
            if (bondIterators != null) {
                IAtomGroup molecule = newAtom.getParentGroup();
                while (!(molecule.getParentGroup() instanceof SpeciesAgent)) {
                    molecule = molecule.getParentGroup();
                }
                
                for (int i=0; i<bondIterators.length; i++) {
                    Potential bondedPotential = bondIterators[i].getPotential();
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
                    iterator.setBasis(molecule);
                    iterator.setTarget(newAtom);
                    iterator.reset();
                    while (iterator.hasNext()) {
                        AtomSet bondedAtoms = iterator.next();
                        Object bond = bondManager.makeBond(bondedAtoms, bondedPotential);
                        bondList.add(bond);
                    }
                }
            }
            return bondList;
        }
        return null;
    }
    
    public void releaseAgent(Object agent, IAtom atom) {
        // we only release a bond when the "up" atom from the bond goes away
        // so if only the "down" atom goes away, we would leave the bond in
        // (bad).  However, you're not allowed to mutate the model, so deleting
        // one leaf atom of a covalent bond and not the other would be illegal.
        ArrayList list = (ArrayList)agent;
        for (int i=0; i<list.size(); i++) {
            bondManager.releaseBond(list.get(i));
        }
        ((ArrayList)agent).clear();
    }
    
    private static final long serialVersionUID = 1L;
    protected final Phase phase;
    protected final AtomAgentManager atomAgentManager;
    protected final HashMap bondIteratorsHash;
    protected BondManager bondManager;
}
