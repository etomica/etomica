package etomica;

import java.util.Random;

public final class IntegratorVelocityVerlet extends IntegratorMD {

    AtomPair.Iterator pairIterator;
    Atom.Iterator atomIterator;
    AtomPair.Action forceSum;
    
    //Fields for Andersen thermostat
    Random random = new Random();
    double nu = 0.001;  //heat bath "collision" frequency
                
    public IntegratorVelocityVerlet() {
        this(Simulation.instance);
    }
    public IntegratorVelocityVerlet(Simulation sim) {
        super(sim);
        
    //anonymous class to sum forces
        forceSum = new AtomPair.Action() {
            private final Space.Vector f = parentSimulation().space().makeVector();
            public void action(AtomPair pair) {
                f.E(((Potential.Soft)parentSimulation().getPotential(pair)).force(pair));
                ((Agent)pair.atom1().ia).force.PE(f);
                ((Agent)pair.atom2().ia).force.ME(f);
            }
        };
    }
    
    
	/**
	 * Overrides superclass method to instantiate iterators when iteratorFactory in phase is changed.
	 * Called by Integrator.addPhase and Integrator.iteratorFactorObserver.
	 */
	protected void makeIterators(IteratorFactory factory) {
        pairIterator = factory.makeAtomPairIteratorAll();
        atomIterator = factory.makeAtomIteratorUp();
    }
     
  private double t2;
  public final void setTimeStep(double t) {
    super.setTimeStep(t);
    t2 = timeStep*timeStep;
  }
  

          
//--------------------------------------------------------------
// steps all particles across time interval tStep

    public void doStep() {

        atomIterator.reset();              //reset iterator of atoms
        while(atomIterator.hasNext()) {    //loop over all atoms
            Atom a = atomIterator.next();  //  advancing positions full step
            Agent agent = (Agent)a.ia;     //  and momenta half step
            Space.Vector r = a.position();
            Space.Vector p = a.momentum();
            p.PEa1Tv1(0.5*timeStep,agent.force);  // p += f(old)*dt/2
            r.PEa1Tv1(timeStep*a.rm(),p);         // r += p*dt/m
            agent.force.E(0.0);
        }
        
        //Add in forces on each atom due to interaction with fields acting in the phase
        for(PotentialField f=firstPhase.firstField(); f!=null; f=f.nextField()) {
            if(!(f instanceof PotentialField.Soft)) continue;
            PotentialField.Soft field = (PotentialField.Soft)f;
            Atom.Iterator iterator = f.getAffectedAtoms();  //iterator for atoms under the influence of this field
            iterator.reset();
            while(iterator.hasNext()) {
                Atom a = iterator.next();
                ((Agent)a.ia).force.PE(field.force(a));
            }
        }
        
        //Add in forces on each atom due to interaction with other atoms in phase
        pairIterator.allPairs(forceSum);    //compute forces on each atom
        
        //Finish integration step
        atomIterator.reset();
        while(atomIterator.hasNext()) {     //loop over atoms again
            Atom a = atomIterator.next();   //  finishing the momentum step
            a.momentum().PEa1Tv1(0.5*timeStep,((Agent)a.ia).force);  //p += f(new)*dt/2
        }
        if(isothermal) {  //Andersen thermostat
            atomIterator.reset();
            double nut = nu*timeStep;
            while(atomIterator.hasNext()) {
                Atom a = atomIterator.next();
                if(random.nextDouble() < nut) a.randomizeMomentum(temperature);  //this method in Atom needs some work
            }
        }
        return;
    }
    
    public void setAndersenNu(double n) {nu = n;}
    public double getAndersenNu() {return nu;}
            
//--------------------------------------------------------------

    protected void doReset() {
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.next();
            Agent agent = (Agent)a.ia;
            agent.force.E(0.0);
        }
        pairIterator.allPairs(forceSum);
    }
              
//--------------------------------------------------------------

    public final Integrator.Agent makeAgent(Atom a) {
        return new Agent(parentSimulation(),a);
    }
            
    public final static class Agent implements Integrator.Agent {  //need public so to use with instanceof
        public Atom atom;
        public Space.Vector force;

        public Agent(Simulation sim, Atom a) {
            atom = a;
            force = sim.space().makeVector();
        }
    }
}

