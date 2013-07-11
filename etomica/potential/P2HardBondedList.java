package etomica.potential;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
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
 * be configured externally to this class.
 * 
 * @author David Kofke
 *
 */
public class P2HardBondedList extends Potential2 implements PotentialHard, AgentSource {
        
    public P2HardBondedList(Simulation sim, Potential2HardSpherical bondedPotential, 
            Potential2HardSpherical nonBondedPotential) {
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
    
    protected boolean isBonded(IAtomList pair) {
        return ((AtomArrayList)agentManager.getAgent(pair.getAtom(0))).contains(pair.getAtom(1));
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
    

    private final Potential2HardSpherical bondedPotential, nonBondedPotential;
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
