package etomica.integrator;

import etomica.Atom;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Integrator;
import etomica.IteratorDirective;
import etomica.Phase;
import etomica.PotentialMaster;
import etomica.Space;
import etomica.atom.iterator.AtomIteratorList;
import etomica.potential.PotentialCalculationForceSum;
import etomica.space.ICoordinateKinetic;
import etomica.space.Vector;

/* History
 * 
 * 06/18/03 (DAK) changed doReset so that rMrLast is given by dt*p/m instead of
 * -dt*p/m
 */
public final class IntegratorVerlet extends IntegratorMD implements EtomicaElement {

    private final AtomIteratorList atomIterator = new AtomIteratorList();
    
    public final PotentialCalculationForceSum forceSum;
    private final IteratorDirective allAtoms = new IteratorDirective();
    private final Space space;
    
    Vector work;
                
    public IntegratorVerlet(PotentialMaster potentialMaster, Space space) {
        super(potentialMaster);
        this.space = space;
        forceSum = new PotentialCalculationForceSum(space);
        work = space.makeVector();
        setTimeStep(etomica.units.systems.LJ.SYSTEM.time().toSim(2.0));
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Molecular dynamics using basic Verlet algorithm");
        return info;
    }

	/**
	 * Overrides superclass method to instantiate iterators when iteratorFactory in phase is changed.
	 * Called by Integrator.addPhase and Integrator.iteratorFactorObserver.
	 */
	public boolean addPhase(Phase p) {
	    if(!super.addPhase(p)) return false;
        atomIterator.setList(p.speciesMaster.atomList);
        return true;
    }
    
        
	private double t2;

    public final void setTimeStep(double t) {
        super.setTimeStep(t);
        t2 = timeStep * timeStep;
    }
          
//--------------------------------------------------------------
// steps all particles across time interval tStep

    public void doStep() {

        //Compute forces on each atom
        atomIterator.reset();
        while(atomIterator.hasNext()) {   //zero forces on all atoms
            ((Agent)atomIterator.nextAtom().ia).force.E(0.0);
        }
        potential.calculate(firstPhase, allAtoms, forceSum);

        //take step
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.nextAtom();
            Agent agent = (Agent)a.ia;
            Vector r = a.coord.position();
            work.E(r);
            r.PE(agent.rMrLast);
            agent.force.TE(a.type.rm()*t2);
            r.PE(agent.force);
            agent.rMrLast.E(r);
            agent.rMrLast.ME(work);
        }
    }//end of doStep
    

    public void reset() {
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.nextAtom();
            Agent agent = (Agent)a.ia;
            agent.rMrLast.Ea1Tv1(timeStep,((ICoordinateKinetic)a.coord).velocity());//06/13/03 removed minus sign before timeStep
        }
        super.reset();
    }
              
//--------------------------------------------------------------

    public final Object makeAgent(Atom a) {
        return new Agent(space,a);
    }
            
	public final static class Agent implements Integrator.Forcible {  //need public so to use with instanceof
        public Atom atom;
        public Vector force;
        public Vector rMrLast;  //r - rLast

        public Agent(Space space, Atom a) {
            atom = a;
            force = space.makeVector();
            rMrLast = space.makeVector();
        }
        
        public Vector force() {return force;}
    }//end of Agent
    
/*    public static void main(String[] args) {
        
	    IntegratorVerlet integrator = new IntegratorVerlet();
	    SpeciesSpheres species = new SpeciesSpheres();
	    Phase phase = new Phase();
	    P2LennardJones potential = new P2LennardJones();
	    Controller controller = new Controller();
	    DisplayPhase display = new DisplayPhase();
	    IntegratorMD.Timer timer = integrator.new Timer(integrator.chronoMeter());
	    timer.setUpdateInterval(10);

		MeterEnergy energy = new MeterEnergy();
		energy.setPhase(phase);
		energy.setHistorying(true);
		energy.setActive(true);		
		energy.getHistory().setNValues(500);		
		DisplayPlot plot = new DisplayPlot();
		plot.setLabel("Energy");
		plot.setDataSource(energy.getHistory());
		
		integrator.setSleepPeriod(2);
		
		Simulation.instance.elementCoordinator.go();
        Simulation.makeAndDisplayFrame(Simulation.instance);

    }//end of main 
    */
}//end of IntegratorVerlet

