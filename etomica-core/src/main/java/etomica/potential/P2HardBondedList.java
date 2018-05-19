/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.box.BoxAgentManager;
import etomica.box.BoxAgentManager.BoxAgentSource;
import etomica.potential.P2HardBondedList.BondArrayList;
import etomica.simulation.Simulation;
import etomica.space.Tensor;

/**
 * Hard potential that wraps two others.  One potential applies to atoms that are found on each other's
 * list of bonded atoms, and the other applies for pairs not on each other's list. The list of bonded atoms must
 * be configured externally to this class, using the setBonded method.
 * 
 * @author David Kofke
 *
 */
public class P2HardBondedList extends Potential2 implements PotentialHard, AgentSource<BondArrayList> {
        
    public P2HardBondedList(Simulation sim, PotentialHard bondedPotential, 
            PotentialHard nonBondedPotential) {
        super(sim.getSpace());
        this.bondedPotential = bondedPotential;
        this.nonBondedPotential = nonBondedPotential;
        
        //box agent manager is used to handle multiple bonded-atoms lists across boxes
        MyBoxAgentSource bas = new MyBoxAgentSource(this);
        boxAgentManager = new BoxAgentManager<AtomLeafAgentManager<BondArrayList>>(bas, sim);
    }

    public double getRange() {
        return Math.max(bondedPotential.getRange(), nonBondedPotential.getRange());
    }

    public void bump(IAtomList pair, double falseTime) {
        MyBond bond = getBond(pair);
        lastCollisionIsBonded = (bond != null);
        if(lastCollisionIsBonded) {
//            bondedPotential.setBondLengthSquared(bond.bondLengthSquared);
            ((P2PenetrableSquareWell)bondedPotential).setCoreDiameterSquared(bond.bondLengthSquared);
            bondedPotential.bump(pair, falseTime);
        } else {
            nonBondedPotential.bump(pair, falseTime);
        }
    }

    public double collisionTime(IAtomList pair, double falseTime) {
        MyBond bond = getBond(pair);
        if(bond != null) {
//            bondedPotential.setBondLengthSquared(bond.bondLengthSquared);
            ((P2PenetrableSquareWell)bondedPotential).setCoreDiameterSquared(bond.bondLengthSquared);
            return bondedPotential.collisionTime(pair, falseTime);
        }
        return nonBondedPotential.collisionTime(pair, falseTime);
    }

    public double energyChange() {
        return lastCollisionIsBonded ? bondedPotential.energyChange() : nonBondedPotential.energyChange();
    }

    public double lastCollisionVirial() {
        return lastCollisionIsBonded ? bondedPotential.lastCollisionVirial() : nonBondedPotential.lastCollisionVirial();
    }

    public Tensor lastCollisionVirialTensor() {
        return lastCollisionIsBonded ? bondedPotential.lastCollisionVirialTensor() : nonBondedPotential.lastCollisionVirialTensor();
    }
    
    public double energy(IAtomList pair) {
        return isBonded(pair) ? bondedPotential.energy(pair) : nonBondedPotential.energy(pair);
    }

    public void setBox(Box box) {
        agentManager = boxAgentManager.getAgent(box);
        bondedPotential.setBox(box);
        nonBondedPotential.setBox(box);
    }
    
    /**
     * Returns true if the pair of atoms are on each others' bonded list; false otherwise.
     */
    public boolean isBonded(IAtom atom0, IAtom atom1) {
        return getBondedList(atom0).contains(atom1);
    }
    public boolean isBonded(IAtomList pair) {
        return isBonded(pair.get(0), pair.get(1));
    }
    
    private MyBond getBond(IAtomList pair) {
        return getBondedList(pair.get(0)).getBond(pair.get(1));
    }
    
    /**
     * Adds/removes the pair of atoms to/from each other's bonded list.
     * Adding a bond to bonded atoms does not cause an error, and has no effect.
     * Remove a bond from unbonded atoms causes an exception.
     * 
     * @param bonding true for adding to bonded list; false for removing
     * @throws IndexOutOfBoundsException if specifying unbonding, and atoms are not bonded already
     */
    public void setBonded(boolean bonding, IAtom atom0, IAtom atom1, double bondLengthSquared) {
        BondArrayList list0 = getBondedList(atom0);
        BondArrayList list1 = getBondedList(atom1);
        if(bonding) {
            if(isBonded(atom0, atom1)) return;
            list0.add(new MyBond(atom1, bondLengthSquared));
            list1.add(new MyBond(atom0, bondLengthSquared));
        } else {
            list0.remove(list0.indexOf(atom1));//throws exception if not in list
            list1.remove(list1.indexOf(atom0));
        }
    }
    
    /**
     * Returns the list containing the bonds for the given atom.
     */
    public BondArrayList getBondedList(IAtom atom) {
        return agentManager.getAgent(atom);
    }

    public BondArrayList makeAgent(IAtom a, Box agentBox) {
        return new BondArrayList();
    }

    public void releaseAgent(BondArrayList agent, IAtom atom, Box agentBox) {
        agent.clear();
    }
    

    private final PotentialHard bondedPotential, nonBondedPotential;
    private boolean lastCollisionIsBonded;
    private AtomLeafAgentManager<BondArrayList> agentManager;
    private final BoxAgentManager<AtomLeafAgentManager<BondArrayList>> boxAgentManager;
    private static final long serialVersionUID = 1L;
    
    //inner class to define BoxAgentSource.  Need to do this instead of implementing
    //BoxAgentSource in P2HardBondedList, because of conflicting definitions of getAgentClass
    private class MyBoxAgentSource implements BoxAgentSource<AtomLeafAgentManager<BondArrayList>> {

        AgentSource<BondArrayList> as;
        public MyBoxAgentSource(AgentSource<BondArrayList> asource) {
            as = asource;
        }

        public AtomLeafAgentManager<BondArrayList> makeAgent(Box box) {
            return new AtomLeafAgentManager<BondArrayList>(as, box);
        }

        public void releaseAgent(AtomLeafAgentManager<BondArrayList> agent) {
            agent.dispose();
            
        }
    }
    
    private class MyBond {
        final double bondLengthSquared;
        final IAtom partner;
        
        MyBond(IAtom partner, double bondLengthSquared) {
            this.partner = partner;
            this.bondLengthSquared = bondLengthSquared;
        }
    }
    
    public class BondArrayList {
        MyBond[] bondList;
        int itemsInList = 0;
        
        public BondArrayList() {
            bondList = new MyBond[10];
        }
        
        public int size() {
            return itemsInList;
        }

        public void clear() {
            for(int i = 0; i < itemsInList; i++) {
                bondList[i] = null;
            }
            itemsInList = 0;
        }

        public int indexOf(IAtom elem) {
            for(int i = 0; i < itemsInList; i++) {
                if(elem == bondList[i].partner) {
                    return i;
                }
            }
            return -1;
        }
        
        public MyBond getBond(IAtom elem) {
            int index = indexOf(elem);
            if(index == -1) return null;
            return bondList[index];
        }
        
        public boolean contains(IAtom elem) {
            return (indexOf(elem) != -1);
        }

        public boolean add(MyBond bond) {

            if(itemsInList == bondList.length) {
                MyBond[] tempList = new MyBond[(int)((float)itemsInList * (1.0f + 0.3f)+1)];

                for(int i = 0; i < bondList.length; i++) {
                    tempList[i] = bondList[i];
                }
                bondList = tempList;
            }
            bondList[itemsInList] = bond; 
            itemsInList++;

            return true;
        }

        public MyBond remove(int index) {
            MyBond bond = null;

            if(index < 0 || index >= itemsInList) {
                throw new IndexOutOfBoundsException("AtomLeafArrayList.remove invalid index");
            }

            bond = bondList[index];
            for(int i = index; i < itemsInList-1; i++) {
                bondList[i] = bondList[i+1];
            }
            bondList[itemsInList-1] = null;
            itemsInList--;

            return bond;
        }

    }
 
}
