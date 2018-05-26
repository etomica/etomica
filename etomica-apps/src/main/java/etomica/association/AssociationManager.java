/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomIterator;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveEvent;
import etomica.integrator.mcmove.MCMoveTrialCompletedEvent;
import etomica.nbr.NeighborIterator;
import etomica.nbr.cell.Api1ACell;
import etomica.nbr.cell.NeighborIteratorCell;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.util.IEvent;
import etomica.util.IListener;

/**
 * Class to define and track atom associations.  Constructed given an iterator
 * that defines the set of atoms that is managed, and an association definition
 * that is used to determine if atoms are associated.  Only pairwise associations
 * are considered.
 * To incorporate an instance of this class in a simulation, add an instance of
 * this EnergySum inner class to the simulation:<br>
 * sim.setEnergySum(associationManager.new EnergySum());
 * <p>
 * If a MC simulation, this should be done before any MCMove classes are instantiated.
 * Also, register this as a listener to IntegratorMC:<br>
 * integrator.addMCMoveListener(associationManager);
 */
public class AssociationManager implements AgentSource<AtomArrayList>, IListener {

    private final Box box;
    private final AtomLeafAgentManager<AtomArrayList> agentManager;
    private final NeighborIterator neighborIterator;
    private final AtomArrayList associatedAtoms = new AtomArrayList();
    private AssociationDefinition associationDefinition;

    public AssociationManager(Box box, PotentialMasterCell potentialMaster, AssociationDefinition definition) {
        this.box = box;
        agentManager = new AtomLeafAgentManager<AtomArrayList>(this, box);
        associationDefinition = definition;
        this.neighborIterator = new NeighborIteratorCell(potentialMaster.getNbrCellManager(box));
    }

    public AssociationDefinition getAssociationDefinition() {
        return associationDefinition;
    }

    public void initialize() {
        IAtomList atomList = box.getLeafList();//list of all atoms in this box
        for (int i = 0; i < atomList.size(); i += 1) {
            IAtom atomi = atomList.get(i);
            agentManager.getAgent(atomi).clear();
        }
        for (int i = 0; i < atomList.size() - 1; i += 1) {
            IAtom atomi = atomList.get(i);//definition of atom i
            for (int j = i + 1; j < atomList.size(); j += 1) {
                IAtom atomj = atomList.get(j);
                if (associationDefinition.isAssociated(atomi, atomj)) {
                    agentManager.getAgent(atomi).add(atomj);//i and j are associated
                    agentManager.getAgent(atomj).add(atomi);
                }
            }
        }
        associatedAtoms.clear();
        for (int i = 0; i < atomList.size(); i += 1) {
            IAtom atomi = atomList.get(i);
            if (agentManager.getAgent(atomi).size() > 0) {
                associatedAtoms.add(atomi);
            }
        }

    }

    public IAtomList getAssociatedAtoms() {
        return associatedAtoms;
    }

    //need also to handle associatedAtoms list
    public void actionPerformed(IEvent evt) {
        MCMoveEvent mcEvent = (MCMoveEvent) evt;
        if (mcEvent instanceof MCMoveTrialCompletedEvent && ((MCMoveTrialCompletedEvent) mcEvent).isAccepted()) {
            return;
        }

        AtomIterator iterator = mcEvent.getMCMove().affectedAtoms(box);
        iterator.reset();
        for (IAtom atomi = iterator.nextAtom(); atomi != null; atomi = iterator.nextAtom()) {
            AtomArrayList listi = agentManager.getAgent(atomi);//list of the old bonds
            if (listi.size() > 0) {
                for (int i = 0; i < listi.size(); i++) {
                    IAtom atomj = listi.get(i);
                    AtomArrayList listj = agentManager.getAgent(atomj);
                    listj.remove(listj.indexOf(atomi));//remove atom i from the listj
                    if (listj.size() == 0) {
                        associatedAtoms.remove(associatedAtoms.indexOf(atomj));
                    }
                }
                listi.clear();//remove all the elements from listi
                associatedAtoms.remove(associatedAtoms.indexOf(atomi));
            }

            IAtom targetAtom = atomi;
            neighborIterator.forEachNeighbor(atomi, null, (atom1, atom2) -> {
                IAtom atomj = atom1 == targetAtom ? atom1 : atom2;
                if (associationDefinition.isAssociated(targetAtom, atomj)) {
                    if (listi.isEmpty()) {
                        associatedAtoms.add(targetAtom);
                    }
                    listi.add(atomj);
                    AtomArrayList listj = agentManager.getAgent(atomj);
                    if (listj.isEmpty()) {
                        associatedAtoms.add(atomj);
                    }
                    listj.add(targetAtom);
                }
            });
        }
//        int numAssociatedAtoms = associatedAtoms.getAtomCount();
//        initialize();
//        if (numAssociatedAtoms != associatedAtoms.getAtomCount()){
//        	throw new RuntimeException ("-_-;");
//        }
        for (IAtom atomi = iterator.nextAtom(); atomi != null; atomi = iterator.nextAtom()) {
            if (getAssociatedAtoms(atomi).size() > 1) {
                System.out.println("atomi" + atomi);
            } else if (getAssociatedAtoms(atomi).size() == 1) {
                IAtom atomj = getAssociatedAtoms(atomi).get(0);
                if (getAssociatedAtoms(atomj).size() == 0) {
                    System.out.println("Wrong");
                } else if (getAssociatedAtoms(atomj).size() > 1) {
                    AtomArrayList listj = agentManager.getAgent(atomj);
                    System.out.println("Wrong:smer");
                    System.out.println("listj = " + listj + " atomi= " + atomi + " atomj= " + atomj);
                }
            }
        }
    }


    /**
     * Returns the number of atoms on the list of associations of the given atom.
     */
    public IAtomList getAssociatedAtoms(IAtom atom) {
        return agentManager.getAgent(atom);
    }

    public AtomArrayList makeAgent(IAtom a, Box agentBox) {
        return new AtomArrayList();
    }

    public void releaseAgent(AtomArrayList agent, IAtom atom, Box agentBox) {
    }
}
