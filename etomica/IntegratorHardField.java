package etomica;

/**
 * Extension of IntegratorHard for case where a constant external force field is applied.
 *
 * @see IntegratorHard
 * @author David Kofke
 *
 */
public final class IntegratorHardField extends IntegratorHard implements EtomicaElement {

//    public final IntegratorHardField.ForceSum forceSum;
	public final PotentialCalculationForceSum forceSum;
    private final IteratorDirective fieldsOnly = new IteratorDirective();
	
	private final IteratorDirective.PotentialCriterion noFieldsCriterion = new IteratorDirective.PotentialCriterion() {
			public boolean excludes(Potential potential) {
				return (potential instanceof Potential1);
			}
	};

    public IntegratorHardField() {
        this(Simulation.instance);
    }
    public IntegratorHardField(Simulation sim) {
        super(sim);
        forceSum = new PotentialCalculationForceSum(space);//new IntegratorHardField.ForceSum(sim.space());
        fieldsOnly.addCriterion(new IteratorDirective.PotentialCriterion() {
            public boolean excludes(Potential potential) {
                return !(potential instanceof Potential1);
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
            Atom a = atomIterator.nextAtom();
            Agent agent = (Agent)a.ia;
            agent.decrementCollisionTime(tStep);
            if(a.coord.isStationary()) {continue;}  //skip if atom is stationary
            a.coord.freeFlight(tStep);
            if(!agent.forceFree) {
//                System.out.println("IntegratorHardField "+agent.force.toString()+" "+a.toString());
                a.coord.translateBy(t2*a.coord.rm(),agent.force);
                a.coord.accelerateBy(tStep,agent.force);
            }
        }
    }
    
    public void reset() {
        calculateForces();
        super.reset();
    }

    private void calculateForces() {
        
        //Compute all forces
        atomIterator.reset();
        while(atomIterator.hasNext()) {   //zero forces on all atoms
            Agent iagent = (Agent)atomIterator.nextAtom().ia;
            iagent.force.E(0.0);
            iagent.forceFree = true;
        }
        //Compute forces on each atom
        potential.calculate(firstPhase, fieldsOnly, forceSum);
        
    }//end of calculateForces
    
    /**
    *
    */
    public void scaleMomenta(double s) {
        double rs = 1.0/s;
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.nextAtom();
            if(a.type instanceof AtomType.Wall &&
                !((Agent)a.ia).forceFree) continue;
            a.coord.momentum().TE(s); //scale momentum
            ((Agent)a.ia).collisionTime *= rs;
        }
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.nextAtom();
 //           System.out.println(a.coord.position().toString()+a.coord.momentum().toString()+"  "+
 //                               a.coord.momentum().squared());
            Agent iagent = (Agent)a.ia;
            if(!iagent.forceFree) updateAtom(a);//update because not force free
            Atom partner = iagent.collisionPartner();
            if(partner == null) continue;
            Agent jagent = (Agent)partner.ia;
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
    

    /**
    * Produces the Agent defined by this integrator.
    * One instance of an Agent is placed in each atom controlled by this integrator.
    */
    public final Object makeAgent(Atom a) {
        return new Agent(simulation(),a);
    }
     
    /**
    * Extends IntegratorHard.Agent to hold a force vector.
    */
    public static class Agent extends IntegratorHard.Agent implements Integrator.Forcible { 
    
        public final Space.Vector force;
        public boolean forceFree = true;
        public Agent(Simulation sim, Atom a) {
            super(a);
            force = sim.space().makeVector();
        }
        public final Space.Vector force() {return force;}
    }//end of Agent
    
    /**
     * Sums the force on each iterated atom and adds it to the integrator agent
     * associated with the atom.
     * Differs from PotentialCalculation.ForceSum in that only 1-body potentials
     * are considered, and also sets forceFree flag of Agent appropriately.
     */
    public static final class PotentialCalculationForceSum extends etomica.PotentialCalculationForceSum {
        
        public PotentialCalculationForceSum(Space space) {
             super(space);
        }
        
		public void actionPerformed(Atom atom) {
			super.actionPerformed(atom);
			((Agent)atom.ia).forceFree = false;
		}
//		public void actionPerformed(AtomPair pair) {
//			throw new etomica.exception.MethodNotImplementedException();
//		}
//		public void actionPerformed(Atom3 atom3) {
//			throw new etomica.exception.MethodNotImplementedException();
//		}
    }//end ForceSums

}//end of IntegratorHardField

