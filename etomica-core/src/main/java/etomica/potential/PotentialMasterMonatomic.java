/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.*;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.simulation.Simulation;
import etomica.species.ISpecies;
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
public class PotentialMasterMonatomic extends PotentialMaster implements AtomTypeAgentManager.AgentSource<PotentialArrayByType> {

    protected final AtomTypeAgentManager<PotentialArrayByType> potentialAgentManager;
    protected final AtomSetSinglet atomSetSinglet;
    protected final AtomPair atomPair;
    protected IPotential[] allPotentials = new IPotential[0];
    
    public PotentialMasterMonatomic(Simulation sim) {
        super();
        potentialAgentManager = new AtomTypeAgentManager<>(this, sim);
        atomSetSinglet = new AtomSetSinglet();
        atomPair = new AtomPair();
    }
    
    public void addPotential(IPotentialMolecular potential, ISpecies[] species) {
        throw new RuntimeException("Probably not the method you really wanted to call");
    }

    public void addPotential(IPotentialAtomic potential, AtomType[] atomTypes) {
        if (potential.nBody() > 2) {
            throw new RuntimeException("I only understand 1 and 2-body potentials");
        }
        allPotentials = Arrays.addObject(allPotentials, potential);
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
            PotentialArrayByType potentialArray = potentialAgentManager.getAgent(atomTypes[i]);
            AtomType otherType = null;
            if (potential.nBody() == 2) {
                otherType = atomTypes[1-i];
            }
            potentialArray.addPotential(potential, otherType);
        }
        if(potential instanceof PotentialTruncated) {
            Potential0Lrc lrcPotential = ((PotentialTruncated) potential).makeLrcPotential(atomTypes);
            if(lrcPotential != null) {
                lrcMaster().addPotential(lrcPotential);
            }
        }
    }

    public void removePotential(IPotentialAtomic potential) {
        super.removePotential(potential);
        for (PotentialArrayByType potentialArray : this.potentialAgentManager.getAgents().values()) {
            potentialArray.removePotential(potential);
        }
        allPotentials = Arrays.removeObject(allPotentials,potential);
    }

    public void calculate(Box box, IteratorDirective id, PotentialCalculation pc) {
        if(!enabled) return;
        IAtom targetAtom = id.getTargetAtom();
        IMolecule targetMolecule = id.getTargetMolecule();

        IAtomList leafList = box.getLeafList();
        if (targetAtom != null || targetMolecule != null) {
            IAtom leafAtom = targetAtom;
            if (leafAtom == null) {
                leafAtom = targetMolecule.getChildList().getAtom(0);
            }
            final int targetIndex = leafAtom.getLeafIndex();
            final PotentialArrayByType potentialArray = potentialAgentManager.getAgent(leafAtom.getType());
            IPotential[] potentials = potentialArray.getPotentials();
            for(int i=0; i<potentials.length; i++) {
                potentials[i].setBox(box);
            }
            calculate(leafAtom, leafList, targetIndex, potentialArray, id.direction(), pc);
        }
        else {
            // invoke setBox on all potentials
            for (int i=0; i<allPotentials.length; i++) {
                allPotentials[i].setBox(box);
            }
            for (int i=0; i<leafList.getAtomCount(); i++) {
                IAtom atom = leafList.getAtom(i);
                PotentialArrayByType potentialArray = potentialAgentManager.getAgent(atom.getType());
                calculate(atom, leafList, i, potentialArray, IteratorDirective.Direction.UP, pc);
            }
        }
        if(lrcMaster != null) {
            lrcMaster.calculate(box, id, pc);
        }
    }

    protected void calculate(IAtom leafAtom, IAtomList leafList, int leafIndex, PotentialArrayByType potentialArray, IteratorDirective.Direction direction, PotentialCalculation pc) {

        IPotential[] potentials = potentialArray.getPotentials();
        int leafCount = leafList.getAtomCount();
        for(int i=0; i<potentials.length; i++) {
            if (potentials[i].nBody() == 1) {
                atomSetSinglet.atom = leafAtom;
                pc.doCalculation(atomSetSinglet, (IPotentialAtomic)potentials[i]);
            }
        }

        AtomType[] types = potentialArray.getTypes();
        if (direction != IteratorDirective.Direction.DOWN) {
            atomPair.atom0 = leafAtom;
            for (int j=leafIndex+1; j<leafCount; j++) {
                atomPair.atom1 = leafList.getAtom(j);
                AtomType type1 = atomPair.atom1.getType();
                for (int i=0; i<types.length; i++) {
                    if (types[i] == type1) {
                        pc.doCalculation(atomPair, (IPotentialAtomic)potentials[i]);
                        break;
                    }
                }
            }
        }
        if (direction != IteratorDirective.Direction.UP) {
            atomPair.atom1 = leafAtom;
            for (int j=0; j<leafIndex; j++) {
                atomPair.atom0 = leafList.getAtom(j);
                AtomType type0 = atomPair.atom0.getType();
                for (int i=0; i<types.length; i++) {
                    if (types[i] == type0) {
                        pc.doCalculation(atomPair, (IPotentialAtomic)potentials[i]);
                        break;
                    }
                }
            }
        }
    }

    public PotentialArrayByType makeAgent(AtomType type) {
        return new PotentialArrayByType();
    }

    public void releaseAgent(PotentialArrayByType agent, AtomType type) {
    }
}
