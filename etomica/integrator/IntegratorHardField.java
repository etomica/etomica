package etomica.integrator;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.iterator.AtomsetIterator;
import etomica.atom.iterator.IteratorDirective;
import etomica.exception.ConfigurationOverlapException;
import etomica.potential.Potential;
import etomica.potential.Potential1;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.ICoordinateKinetic;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.util.IRandom;

/**
 * Extension of IntegratorHard for case where a constant external force field is applied.
 *
 * @see IntegratorHard
 * @author David Kofke
 *
 */
public final class IntegratorHardField extends IntegratorHard implements EtomicaElement {

    private static final long serialVersionUID = 1L;
	public final PotentialCalculationForceSum forceSum;
    private final IteratorDirective fieldsOnly = new IteratorDirective();
    private final Space space;
	
    //XXX not serializable
    private final IteratorDirective.PotentialCriterion noFieldsCriterion = new IteratorDirective.PotentialCriterion() {
	    public boolean excludes(Potential candidatePotential) {
	        return (candidatePotential instanceof Potential1);
	    }
    };

    public IntegratorHardField(Simulation sim) {
        this(sim.getPotentialMaster(),sim.getRandom(),sim.getDefaults().timeStep,sim.getDefaults().temperature);
    }
    
    public IntegratorHardField(PotentialMaster potentialMaster, IRandom random,
            double timeStep, double temperature) {
        super(potentialMaster,random,timeStep,temperature);
        space = potentialMaster.getSpace();
        forceSum = new PotentialCalculationForceSum();
        //XXX not serializable
        fieldsOnly.addCriterion(new IteratorDirective.PotentialCriterion() {
            public boolean excludes(Potential candidatePotential) {
                return !(candidatePotential instanceof Potential1);
            }
        });
        upList.addCriterion(noFieldsCriterion);
        downList.addCriterion(noFieldsCriterion);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Collision-based MD simulation in the presence of a constant external field");
        return info;
    }

    /**
    * Advances all atom coordinates by tStep, without any intervening collisions.
    * Uses constant-force kinematics.
    */
    protected void advanceAcrossTimeStep(double tStep) {
        
        calculateForces();
        
        double t2 = 0.5*tStep*tStep;
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            AtomLeaf a = (AtomLeaf)atomIterator.nextAtom();
            HardFieldAgent agent = (HardFieldAgent)agentManager.getAgent(a);
            agent.decrementCollisionTime(tStep);
            a.getCoord().getPosition().PEa1Tv1(tStep,((ICoordinateKinetic)a.getCoord()).getVelocity());
            if(!agent.forceFree) {
//                System.out.println("IntegratorHardField "+agent.force.toString()+" "+a.toString());
                a.getCoord().getPosition().PEa1Tv1(t2*((AtomTypeLeaf)a.getType()).rm(),agent.force);
                ((ICoordinateKinetic)a.getCoord()).getVelocity().PEa1Tv1(tStep*((AtomTypeLeaf)a.getType()).rm(),agent.force);
            }
        }
    }
    
    public void reset() throws ConfigurationOverlapException {
        forceSum.setAgentManager(agentManager);
        calculateForces();
        super.reset();
    }

    private void calculateForces() {
        
        //Compute all forces
        atomIterator.reset();
        while(atomIterator.hasNext()) {   //zero forces on all atoms
            HardFieldAgent iagent = (HardFieldAgent)agentManager.getAgent(atomIterator.nextAtom());
            iagent.force.E(0.0);
            iagent.forceFree = true;
        }
        //Compute forces on each atom
        potential.calculate(phase, fieldsOnly, forceSum);
        
    }//end of calculateForces
    
    /**
    *
    */
    public void scaleMomenta(double s) {
        double rs = 1.0/s;
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            AtomLeaf a = (AtomLeaf)atomIterator.nextAtom();
            ((ICoordinateKinetic)a.getCoord()).getVelocity().TE(s); //scale momentum
            ((Agent)agentManager.getAgent(a)).eventLinker.sortKey *= rs;
        }
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.nextAtom();
 //           System.out.println(a.coord.position().toString()+a.coord.momentum().toString()+"  "+
 //                               a.coord.momentum().squared());
            HardFieldAgent iagent = (HardFieldAgent)agentManager.getAgent(a);
            if(!iagent.forceFree) updateAtom(a);//update because not force free
            Atom partner = iagent.collisionPartner();
            if(partner == null) continue;
            HardFieldAgent jagent = (HardFieldAgent)agentManager.getAgent(partner);
            if(!iagent.forceFree) {
                updateAtom(partner);//update because partner not force free
                continue;
            }
            if(!jagent.forceFree) {
                updateAtom(partner);
                updateAtom(a);
            }
        }
        findNextCollider();
 //       System.out.println();
    }
    
    public Class getAgentClass() {
        return HardFieldAgent.class;
    }
    
    /**
    * Produces the Agent defined by this integrator.
    * One instance of an Agent is placed in each atom controlled by this integrator.
    */
    public Object makeAgent(Atom a) {
        return new HardFieldAgent(a,this);
    }
     
    /**
    * Extends IntegratorHard.Agent to hold a force vector.
    */
    public static class HardFieldAgent extends IntegratorHard.Agent implements IntegratorPhase.Forcible { 
    
        public final IVector force;
        public boolean forceFree = true;
        public HardFieldAgent(Atom a, IntegratorHardField integrator) {
            super(a, integrator);
            force = integrator.space.makeVector();
        }
        public final IVector force() {return force;}
    }//end of Agent
    
    /**
     * Sums the force on each iterated atom and adds it to the integrator agent
     * associated with the atom.
     * Differs from PotentialCalculation.ForceSum in that only 1-body potentials
     * are considered, and also sets forceFree flag of Agent appropriately.
     */
    public static final class PotentialCalculationForceSum extends etomica.potential.PotentialCalculationForceSum {

		public void doCalculation(AtomsetIterator iterator, Potential potential) {
			super.doCalculation(iterator,potential);
            iterator.reset();
            while(iterator.hasNext()) {
                AtomSet atoms = iterator.next();
                ((HardFieldAgent)integratorAgentManager.getAgent(atoms.getAtom(0))).forceFree = false;
            }
		}
    }//end ForceSums

}//end of IntegratorHardField

