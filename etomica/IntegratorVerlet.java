package etomica;

public final class IntegratorVerlet extends IntegratorMD implements EtomicaElement {

    AtomPair.Iterator pairIterator;
    Atom.Iterator atomIterator;
    AtomPair.Action forceSum;
    AtomAction verletStep;
    Space.Vector work;
                
    public IntegratorVerlet() {
        this(Simulation.instance);
    }
    public IntegratorVerlet(final Simulation sim) {
        super(sim);
        
    //anonymous class
        forceSum = new AtomPair.Action() {
            private Space.Vector f = sim.space().makeVector();
            public void action(AtomPair pair) {
                f.E(((Potential.Soft)parentSimulation().getPotential(pair)).force(pair));
                ((Agent)pair.atom1().ia).force.PE(f);
                ((Agent)pair.atom2().ia).force.ME(f);
            }
        };
        work = sim.space().makeVector();
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Molecular dynamics using basic Verlet algorithm");
        return info;
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

        atomIterator.reset();
        while(atomIterator.hasNext()) {   //zero forces on all atoms
            ((Agent)atomIterator.next().ia).force.E(0.0);
        }
        //Add in forces on each atom due to interaction with fields acting in the phase
        for(PotentialField f=firstPhase.firstField(); f!=null; f=f.nextField()) {
            PotentialField.Soft field = (PotentialField.Soft)f;
            Atom.Iterator iterator = f.getAffectedAtoms();  //iterator for atoms under the influence of this field
            iterator.reset();
            while(iterator.hasNext()) {
                Atom a = iterator.next();
                ((Agent)a.ia).force.PE(field.force(a));
            }
        }
        
        //Add in forces on each atom due to interaction with other atoms in phase
        pairIterator.allPairs(forceSum);

        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.next();
            Agent agent = (Agent)a.ia;
            Space.Vector r = a.position();
            work.E(r);
            r.PE(agent.rMrLast);
            agent.force.TE(a.rm()*t2);
            r.PE(agent.force);
            agent.rMrLast.E(r);
            agent.rMrLast.ME(work);
        }
        return;
    }
    
//    private static class forceSum implements AtomPair.Action {
//        public void action(AtomPair pair) {
            
        
//--------------------------------------------------------------

    protected void doReset() {
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.next();
            Agent agent = (Agent)a.ia;
            agent.rMrLast.Ea1Tv1(-timeStep*a.rm(),a.momentum());
        }
    }
              
//--------------------------------------------------------------

    public final Integrator.Agent makeAgent(Atom a) {
        return new Agent(parentSimulation(),a);
    }
            
    public final static class Agent implements Integrator.Agent {  //need public so to use with instanceof
        public Atom atom;
        public Space.Vector force;
        public Space.Vector rMrLast;  //r - rLast

        public Agent(Simulation sim, Atom a) {
            atom = a;
            force = sim.space().makeVector();
            rMrLast = sim.space().makeVector();
        }
    }
}

