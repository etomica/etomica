package etomica.potential;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.box.BoxAgentManager;
import etomica.box.BoxAgentManager.BoxAgentSource;
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
public class P2HardBondedList extends Potential2 implements PotentialHard, AgentSource {
        
    public P2HardBondedList(Simulation sim, PotentialHard bondedPotential, 
            PotentialHard nonBondedPotential) {
        super(sim.getSpace());
        this.bondedPotential = bondedPotential;
        this.nonBondedPotential = nonBondedPotential;
        
        //box agent manager is used to handle multiple bonded-atoms lists across boxes
        BoxAgentSource bas = new MyBoxAgentSource(this);
        boxAgentManager = new BoxAgentManager(bas, sim);
    }

    public double getRange() {
        return Math.max(bondedPotential.getRange(), nonBondedPotential.getRange());
    }

    public void bump(IAtomList pair, double falseTime) {
        lastCollisionIsBonded = isBonded(pair);
        if(lastCollisionIsBonded) {
            bondedPotential.bump(pair, falseTime);
        } else {
            nonBondedPotential.bump(pair, falseTime);
        }
    }

    public double collisionTime(IAtomList pair, double falseTime) {
        return isBonded(pair) ? bondedPotential.collisionTime(pair, falseTime) : nonBondedPotential.collisionTime(pair, falseTime);
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

    public void setBox(IBox box) {
        agentManager = (AtomLeafAgentManager)boxAgentManager.getAgent(box);
        bondedPotential.setBox(box);
        nonBondedPotential.setBox(box);
    }
    
    /**
     * Returns true if the pair of atoms are on each others' bonded list; false otherwise.
     */
    public boolean isBonded(IAtomList pair) {
        return getBondedList(pair.getAtom(0)).contains(pair.getAtom(1));
    }
    
    /**
     * Adds/removes the pair of atoms to/from each other's bonded list.
     * Adding a bond to bonded atoms does not cause an error, and has no effect.
     * Remove a bond from unbonded atoms causes an exception.
     * 
     * @param bonded true for adding to bonded list; false for removing
     * @param pair pair of atoms to be bonded or unbonded
     * @throws IndexOutOfBoundsException if specifying unbonding, and atoms are not bonded already
     */
    public void setBonded(boolean bonding, IAtomList pair) {
        IAtom atom0 = pair.getAtom(0);
        IAtom atom1 = pair.getAtom(1);
        AtomArrayList list0 = getBondedList(atom0);
        AtomArrayList list1 = getBondedList(atom1);
        if(bonding) {
            if(isBonded(pair)) return;
            list0.add(atom1);
            list1.add(atom0);
        } else {
            list0.remove(list0.indexOf(atom1));//throws exception if not in list
            list1.remove(list1.indexOf(atom0));
        }
    }
    
    /**
     * Returns the list containing the bonded atoms for the given atom.
     */
    private AtomArrayList getBondedList(IAtom atom) {
        return (AtomArrayList)agentManager.getAgent(atom);
    }
    
    public Class getAgentClass() {
        return AtomArrayList.class;
    }

    public Object makeAgent(IAtom a) {
        return new AtomArrayList();
    }

    public void releaseAgent(Object agent, IAtom atom) {
        ((AtomArrayList)agent).clear();
    } 
    

    private final PotentialHard bondedPotential, nonBondedPotential;
    private boolean lastCollisionIsBonded;
    private AtomLeafAgentManager agentManager;
    private final BoxAgentManager boxAgentManager;
    private static final long serialVersionUID = 1L;
    
    //inner class to define BoxAgentSource.  Need to do this instead of implementing
    //BoxAgentSource in P2HardBondedList, because of conflicting definitions of getAgentClass
    private class MyBoxAgentSource implements BoxAgentSource {

        AgentSource as;
        public MyBoxAgentSource(AgentSource asource) {
            as = asource;
        }

        public Class getAgentClass() {
            return AtomLeafAgentManager.class;
        }

        public Object makeAgent(IBox box) {
            return new AtomLeafAgentManager(as, box);
        }

        public void releaseAgent(Object agent) {
            ((AtomLeafAgentManager)agent).dispose();
            
        }
    }
 
}
