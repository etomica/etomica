package etomica.potential;

import etomica.api.IAtom;
import etomica.api.IAtomSet;
import etomica.api.IAtomType;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IPotential;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.atom.AtomPair;
import etomica.atom.AtomSetSinglet;
import etomica.atom.AtomTypeAgentManager;
import etomica.atom.iterator.IteratorDirective;
import etomica.space.ISpace;
import etomica.util.Arrays;

/**
 * PotentialMaster suitable for monatomic molecules without using neighbor
 * lists and avoids PotentialMaster's performance problems.  Potentials should
 * be added via the method taking IAtomTypeLeafs.  This currently only handles
 * 1 and 2 body potentials.  This class assumes that only one 2-body potential
 * applies to any given pair of leaf atoms, but more than one 1-body potential
 * may apply to any atom.  PotentialGroups are not created to hold the leaf
 * potentials.  Calling addPotential with an ISpecies array will fail.
 * 
 * @author Andrew Schultz
 */
public class PotentialMasterMonatomic extends PotentialMaster implements AtomTypeAgentManager.AgentSource {

    public PotentialMasterMonatomic(ISimulation sim, ISpace space) {
        super(space);
        potentialAgentManager = new AtomTypeAgentManager(this, sim.getSpeciesManager(), sim.getEventManager(), true);
        potentialIterator = potentialAgentManager.makeIterator();
        atomSetSinglet = new AtomSetSinglet();
        atomPair = new AtomPair();
    }
    
    public void addPotential(IPotential potential, ISpecies[] species) {
        throw new RuntimeException("Probably not the method you really wanted to call");
    }

    public void addPotential(IPotential potential, IAtomTypeLeaf[] atomTypes) {
        if (potential.nBody() > 2) {
            throw new RuntimeException("I only understand 1 and 2-body potentials");
        }
        allPotentials = (IPotential[])etomica.util.Arrays.addObject(allPotentials, potential);
        for (int i=0; i<atomTypes.length; i++) {
            boolean alreadyDone = false;
            for (int j=0; j<i; j++) {
                if (atomTypes[i] == atomTypes[j]) {
                    alreadyDone = true;
                    break;
                }
            }
            if (alreadyDone) {
                continue;
            }
            PotentialArrayByType potentialArray = (PotentialArrayByType)potentialAgentManager.getAgent(atomTypes[i]);
            IAtomType otherType = null;
            if (potential.nBody() == 2) {
                otherType = atomTypes[1-i];
            }
            potentialArray.addPotential(potential, otherType);
        }
    }

    public void removePotential(IPotential potential) {
        super.removePotential(potential);
        potentialIterator.reset();
        while (potentialIterator.hasNext()) {
            ((PotentialArrayByType)potentialIterator.next()).removePotential(potential);
        }
        allPotentials = (IPotential[])Arrays.removeObject(allPotentials,potential);
    }

    public void calculate(IBox box, IteratorDirective id, PotentialCalculation pc) {
        if(!enabled) return;
        IAtom targetAtom = id.getTargetAtom();

        IAtomSet leafList = box.getLeafList();
        if (targetAtom != null) {
            if (targetAtom instanceof IMolecule) {
                targetAtom = ((IMolecule)targetAtom).getChildList().getAtom(0);
            }
            final int targetIndex = box.getLeafIndex(targetAtom);
            final PotentialArrayByType potentialArray = (PotentialArrayByType)potentialAgentManager.getAgent(targetAtom.getType());
            IPotential[] potentials = potentialArray.getPotentials();
            for(int i=0; i<potentials.length; i++) {
                potentials[i].setBox(box);
            }
            calculate(targetAtom, leafList, targetIndex, potentialArray, id.direction(), pc);
        }
        else {
            // invoke setBox on all potentials
            for (int i=0; i<allPotentials.length; i++) {
                allPotentials[i].setBox(box);
            }
            for (int i=0; i<leafList.getAtomCount(); i++) {
                IAtom atom = leafList.getAtom(i);
                PotentialArrayByType potentialArray = (PotentialArrayByType)potentialAgentManager.getAgent(atom.getType());
                calculate(atom, leafList, i, potentialArray, IteratorDirective.Direction.UP, pc);
            }
        }
        if(lrcMaster != null) {
            lrcMaster.calculate(box, id, pc);
        }
    }
    
    protected void calculate(IAtom leafAtom, IAtomSet leafList, int leafIndex, PotentialArrayByType potentialArray, IteratorDirective.Direction direction, PotentialCalculation pc) {
        
        IPotential[] potentials = potentialArray.getPotentials();
        int leafCount = leafList.getAtomCount();
        for(int i=0; i<potentials.length; i++) {
            if (potentials[i].nBody() == 1) {
                atomSetSinglet.atom = leafAtom;
                pc.doCalculation(atomSetSinglet, potentials[i]);
            }
        }

        IAtomType[] types = potentialArray.getTypes();
        if (direction != IteratorDirective.Direction.DOWN) {
            atomPair.atom0 = leafAtom;
            for (int j=leafIndex+1; j<leafCount; j++) {
                atomPair.atom1 = leafList.getAtom(j);
                IAtomType type1 = atomPair.atom1.getType();
                for (int i=0; i<types.length; i++) {
                    if (types[i] == type1) {
                        pc.doCalculation(atomPair, potentials[i]);
                        break;
                    }
                }
            }
        }
        if (direction != IteratorDirective.Direction.UP) {
            atomPair.atom1 = leafAtom;
            for (int j=0; j<leafIndex; j++) {
                atomPair.atom0 = leafList.getAtom(j);
                IAtomType type0 = atomPair.atom0.getType();
                for (int i=0; i<types.length; i++) {
                    if (types[i] == type0) {
                        pc.doCalculation(atomPair, potentials[i]);
                        break;
                    }
                }
            }
        }
    }
    
    public Class getTypeAgentClass() {
        return PotentialArrayByType.class;
    }
    
    public Object makeAgent(IAtomType type) {
        return new PotentialArrayByType();
    }
    
    public void releaseAgent(Object agent, IAtomType type) {
    }

    protected final AtomTypeAgentManager.AgentIterator potentialIterator;
    protected IPotential[] allPotentials = new IPotential[0];
    protected final AtomTypeAgentManager potentialAgentManager;
    protected final AtomSetSinglet atomSetSinglet;
    protected final AtomPair atomPair;
}
