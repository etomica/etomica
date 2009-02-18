package etomica.integrator;

import etomica.EtomicaInfo;
import etomica.api.IAtomKinetic;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IPotential;
import etomica.api.IPotentialAtomic;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.api.IVectorMutable;
import etomica.atom.iterator.IteratorDirective;
import etomica.exception.ConfigurationOverlapException;
import etomica.potential.Potential1;
import etomica.space.ISpace;

/**
 * Extension of IntegratorHard for case where a constant external force field is applied.
 *
 * @see IntegratorHard
 * @author David Kofke
 *
 */
public final class IntegratorHardField extends IntegratorHard {

    private static final long serialVersionUID = 1L;
	public final PotentialCalculationForceSum forceSum;
    private final IteratorDirective fieldsOnly = new IteratorDirective();
//    private final Space space;
	
    //XXX not serializable
    private final IteratorDirective.PotentialCriterion noFieldsCriterion = new IteratorDirective.PotentialCriterion() {
	    public boolean excludes(IPotential candidatePotential) {
	        return (candidatePotential instanceof Potential1);
	    }
    };

    public IntegratorHardField(ISimulation sim, IPotentialMaster potentialMaster, ISpace _space) {
        this(sim, potentialMaster, sim.getRandom(),0.05, 1.0, _space);
    }
    
    public IntegratorHardField(ISimulation sim, IPotentialMaster potentialMaster, IRandom random,
            double timeStep, double temperature, ISpace _space) {
        super(sim, potentialMaster,random,timeStep,temperature, _space);
        forceSum = new PotentialCalculationForceSum();
        //XXX not serializable
        fieldsOnly.addCriterion(new IteratorDirective.PotentialCriterion() {
            public boolean excludes(IPotential candidatePotential) {
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
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            HardFieldAgent agent = (HardFieldAgent)agentManager.getAgent((IAtomLeaf)a);
            agent.decrementCollisionTime(tStep);
            a.getPosition().PEa1Tv1(tStep,a.getVelocity());
            if(!agent.forceFree) {
//                System.out.println("IntegratorHardField "+agent.force.toString()+" "+a.toString());
                a.getPosition().PEa1Tv1(t2*((IAtomLeaf)a).getType().rm(),agent.force);
                a.getVelocity().PEa1Tv1(tStep*((IAtomLeaf)a).getType().rm(),agent.force);
            }
        }
    }
    
    public void setBox(IBox newBox) {
        super.setBox(newBox);
        forceSum.setAgentManager(agentManager);
    }
    
    public void reset() {
        super.reset();
        calculateForces();
    }

    private void calculateForces() {
        
        //Compute all forces
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomLeaf atom = leafList.getAtom(iLeaf);
            ((HardFieldAgent)agentManager.getAgent(atom)).forceFree = true;
        }
        //zero forces on all atoms
        forceSum.reset();
        //Compute forces on each atom
        potentialMaster.calculate(box, fieldsOnly, forceSum);
        
    }//end of calculateForces
    
    /**
    *
    */
    public void scaleMomenta(double s) {
        double rs = 1.0/s;
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            a.getVelocity().TE(s); //scale momentum
            ((Agent)agentManager.getAgent((IAtomLeaf)a)).eventLinker.sortKey *= rs;
        }
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
 //           System.out.println(a.coord.position().toString()+a.coord.momentum().toString()+"  "+
 //                               a.coord.momentum().squared());
            HardFieldAgent iagent = (HardFieldAgent)agentManager.getAgent((IAtomLeaf)a);
            if(!iagent.forceFree) updateAtom((IAtomLeaf)a);//update because not force free
            IAtomLeaf partner = iagent.collisionPartner();
            if(partner == null) continue;
            HardFieldAgent jagent = (HardFieldAgent)agentManager.getAgent(partner);
            if(!iagent.forceFree) {
                updateAtom((IAtomLeaf)partner);//update because partner not force free
                continue;
            }
            if(!jagent.forceFree) {
                updateAtom((IAtomLeaf)partner);
                updateAtom((IAtomLeaf)a);
            }
        }
        findNextCollider();
 //       System.out.println();
    }
    
    public Class getSpeciesAgentClass() {
        return HardFieldAgent.class;
    }
    
    /**
    * Produces the Agent defined by this integrator.
    * One instance of an Agent is placed in each atom controlled by this integrator.
    */
    public Object makeAgent(IAtomLeaf a) {
        return new HardFieldAgent(a,this);
    }
     
    /**
    * Extends IntegratorHard.Agent to hold a force vector.
    */
    public static class HardFieldAgent extends IntegratorHard.Agent implements IntegratorBox.Forcible { 
    
        public final IVectorMutable force;
        public boolean forceFree = true;
        public HardFieldAgent(IAtomLeaf a, IntegratorHardField integrator) {
            super(a, integrator);
            force = integrator.space.makeVector();
        }
        public final IVectorMutable force() {return force;}
    }//end of Agent
    
    /**
     * Sums the force on each iterated atom and adds it to the integrator agent
     * associated with the atom.
     * Differs from PotentialCalculation.ForceSum in that only 1-body potentials
     * are considered, and also sets forceFree flag of Agent appropriately.
     */
    public static final class PotentialCalculationForceSum extends etomica.potential.PotentialCalculationForceSum {

		public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
			super.doCalculation(atoms,potential);
			((HardFieldAgent)integratorAgentManager.getAgent(atoms.getAtom(0))).forceFree = false;
		}
    }
}
