package etomica;

import java.util.Random;

public final class IntegratorVelocityVerlet extends IntegratorMD implements EtomicaElement {

    AtomIterator atomIterator;
    
    public final PotentialCalculation.ForceSum forceSum;
    private final IteratorDirective allAtoms = new IteratorDirective();
    protected PotentialAgent phasePotential;
    
    //Fields for Andersen thermostat
    Random random = new Random();
    double nu = 0.001;  //heat bath "collision" frequency
                
    public IntegratorVelocityVerlet() {
        this(Simulation.instance);
    }
    public IntegratorVelocityVerlet(Simulation sim) {
        super(sim);
        forceSum = new PotentialCalculation.ForceSum(sim.space());
        
        setTimeStep(etomica.units.LennardJones.Time.UNIT.toSim(2.0));
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Molecular dynamics using velocity Verlet integration algorithm");
        return info;
    }

          //need to modify to handle multiple-phase issues
    public boolean addPhase(Phase p) {
        if(!super.addPhase(p)) return false;
        phasePotential = p.potential();
        return true;
    }
    
	/**
	 * Overrides superclass method to instantiate iterators when iteratorFactory in phase is changed.
	 * Called by Integrator.addPhase and Integrator.iteratorFactorObserver.
	 */
	protected void makeIterators(IteratorFactory factory) {
        atomIterator = factory.makeAtomIterator();
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
                
        //Add in forces on each atom due to interaction with other atoms in phase
        phasePotential.calculate(allAtoms, forceSum);
        
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
        phasePotential.calculate(allAtoms, forceSum);
    }
              
//--------------------------------------------------------------

    public final Integrator.Agent makeAgent(Atom a) {
        return new Agent(parentSimulation(),a);
    }
            
    public final static class Agent implements Integrator.Agent.Forcible {  //need public so to use with instanceof
        public Atom atom;
        public Space.Vector force;

        public Agent(Simulation sim, Atom a) {
            atom = a;
            force = sim.space().makeVector();
        }
        
        public Space.Vector force() {return force;}
    }
}

