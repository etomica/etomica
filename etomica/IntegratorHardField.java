package etomica;

/**
 * Extension of IntegratorHard for case where a constant external force field is applied.
 *
 * @see IntegratorHard
 * @author David Kofke
 *
 */
public final class IntegratorHardField extends IntegratorHard implements EtomicaElement {

    public String getVersion() {return "IntegratorHardField:01.03.17/"+super.getVersion();}
    public final IntegratorHardField.ForceSum forceSum;
    private final IteratorDirective fieldsOnly = new IteratorDirective();

    public IntegratorHardField() {
        this(Simulation.instance);
    }
    public IntegratorHardField(Simulation sim) {
        super(sim);
        forceSum = new IntegratorHardField.ForceSum(sim.space());
        fieldsOnly.addCriterion(new IteratorDirective.PotentialCriterion() {
            public boolean excludes(Potential potential) {
                return !(potential instanceof Potential1);
            }
        });
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
        for(Atom a=firstPhase.firstAtom(); a!=null; a=a.nextAtom()) {
            Agent agent = (Agent)a.ia;
            agent.decrementCollisionTime(tStep);
            if(a.coord.isStationary()) {continue;}  //skip if atom is stationary
            a.coord.freeFlight(tStep);
            if(!agent.forceFree) {
//                System.out.println(agent.force.toString()+" "+a.toString());
                a.coord.translateBy(t2*a.coord.rm(),agent.force);
                a.coord.accelerateBy(tStep,agent.force);
            }
        }
    }
    
    protected void doReset() {
        calculateForces();
        super.doReset();
    }

    private void calculateForces() {
        
        //Compute all forces
        atomIterator.reset();
        while(atomIterator.hasNext()) {   //zero forces on all atoms
            Agent iagent = (Agent)atomIterator.next().ia;
            iagent.force.E(0.0);
            iagent.forceFree = true;
        }
        //Compute forces on each atom
        potential.calculate(fieldsOnly, forceSum);
        
    }//end of calculateForces

    /**
    * Produces the Agent defined by this integrator.
    * One instance of an Agent is placed in each atom controlled by this integrator.
    */
    public final Integrator.Agent makeAgent(Atom a) {
        return new Agent(parentSimulation(),a);
    }
     
    /**
    * Extends IntegratorHard.Agent to hold a force vector.
    */
    public final static class Agent extends IntegratorHard.Agent implements Integrator.Agent.Forcible { 
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
    public static final class ForceSum implements Potential1Calculation {
        
        private final Space.Vector f;
        public ForceSum(Space space) {
             f = space.makeVector();
        }
        
        //atom
        public void calculate(AtomIterator iterator, Potential1 potential) {
            Potential1Soft potentialSoft = (Potential1Soft)potential;
            while(iterator.hasNext()) {
                Atom atom = iterator.next();
                f.E(potentialSoft.gradient(atom));
                Agent iagent = ((Agent)atom.ia);
                iagent.force().ME(f);
                iagent.forceFree = false;
            }//end while
        }//end of calculate
    }//end ForceSums

}//end of IntegratorHardField

