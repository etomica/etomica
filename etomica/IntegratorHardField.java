package etomica;

/**
 * Extension of IntegratorHard for case where a constant external force field is applied.
 *
 * @see IntegratorHard
 * @author David Kofke
 *
 */
 
public final class IntegratorHardField extends IntegratorHard implements EtomicaElement {

    public String getVersion() {return "IntegratorHardField:01.03.17.0/"+super.getVersion();}

    public IntegratorHardField() {
        this(Simulation.instance);
    }
    public IntegratorHardField(Simulation sim) {
        super(sim);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Collision-based MD simulation in the presence of a hard external field");
        return info;
    }

    /**
    * Advances all atom coordinates by tStep, without any intervening collisions.
    * Uses constant-force kinematics.
    * (as defined in the phase)
    */
    protected void advanceAcrossTimeStep(double tStep) {
        
        computeForces();
        
        double t2 = 0.5*tStep*tStep;
        for(Atom a=firstPhase.firstAtom(); a!=null; a=a.nextAtom()) {
            Agent agent = (Agent)a.ia;
            agent.decrementCollisionTime(tStep);
            if(a.isStationary()) {continue;}  //skip if atom is stationary
            a.translateBy(tStep*a.rm(),a.momentum());
            if(!agent.forceFree) {
                a.translateBy(t2*a.rm(),agent.force);
                a.accelerateBy(tStep,agent.force);
            }
        }
    }
    
    private void computeForces() {
        //Compute all forces
        for(Atom a=firstPhase.firstAtom(); a!=null; a=a.nextAtom()) {
            ((Agent)a.ia).force.E(0.0);
        }
        //Add in forces on each atom due to interaction with fields acting in the phase
        for(PotentialField f=firstPhase.firstField(); f!=null; f=f.nextField()) {
            if(!(f instanceof PotentialField.Soft)) continue;
            PotentialField.Soft field = (PotentialField.Soft)f;
            Atom.Iterator iterator = f.getAffectedAtoms();  //iterator for atoms under the influence of this field
            iterator.reset();
            while(iterator.hasNext()) {
                Atom a = iterator.next();
                Agent agent = (Agent)a.ia;
                agent.force.PE(field.force(a));
                agent.forceFree = false;
            }
        }
    }//end of computeForces
    
    protected void doReset() {
        computeForces();
        if(isothermal) scaleMomenta(Math.sqrt(this.temperature/(firstPhase.kineticTemperature())));
        Atom.Iterator iterator = firstPhase.iteratorFactory().makeAtomIteratorUp();
        iterator.reset();
        while(iterator.hasNext()) {upList(iterator.next());}
        findNextCollider();
    }

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
}//end of IntegratorHardField

