// includes a main method

package etomica;
import etomica.units.*;

/**
 * Gear 4th-order predictor-corrector integrator.
 *
 * @author Ed Maginn
 * @author David Kofke
 */
public class IntegratorGear4 extends IntegratorMD implements EtomicaElement {

    public String getVersion() {return "IntegratorGear4:01.07.05/"+IntegratorMD.VERSION;}

    protected AtomIterator atomIterator;
    protected final PotentialCalculationForceSum forceSum;
    private final IteratorDirective allAtoms = new IteratorDirective();
    final Space.Vector work1, work2;
    double zeta = 0.0;
    double chi = 0.0;
    double p1, p2, p3, p4;
    double c0, c2, c3, c4;
    
    static final double GEAR0 = 251./720.;
    static final double GEAR2 = 11./12.;
    static final double GEAR3 = 1./3.;
    static final double GEAR4 = 1./24.;
                
    public IntegratorGear4() {
        this(Simulation.instance);
    }
    public IntegratorGear4(final Simulation sim) {
        super(sim);
        forceSum = new PotentialCalculationForceSum(sim.space());
        work1 = sim.space().makeVector();
        work2 = sim.space().makeVector();
        setTimeStep(etomica.units.LennardJones.Time.UNIT.toSim(2.0));
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Molecular dynamics using Gear 4th-order predictor/corrector algorithm");
        return info;
    }

    public void setTimeStep(double dt) {
        super.setTimeStep(dt);
        p1 = dt;
        p2 = p1 * dt / 2.0;
        p3 = p2 * dt / 3.0;
        p4 = p3 * dt / 4.0;
        c0 = GEAR0 * p1;
        c2 = GEAR2 * p1 / p2;
        c3 = GEAR3 * p1 / p3;
        c4 = GEAR4 * p1 / p4;
    }
        

	/**
	 * Overrides superclass method to instantiate iterators when iteratorFactory in phase is changed.
	 * Called by Integrator.addPhase and Integrator.iteratorFactorObserver.
	 */
	protected void makeIterators(IteratorFactory factory) {
        atomIterator = factory.makeAtomIterator();
    }
    
//--------------------------------------------------------------
// steps all particles across time interval tStep

    public void doStep() {
        
        predictor();
        calculateForces();
        corrector();
        
    }//end of doStep
    
    private void calculateForces() {
        
        //Compute all forces
        atomIterator.reset();
        while(atomIterator.hasNext()) {   //zero forces on all atoms
            ((IntegratorGear4.Agent)atomIterator.next().ia).force.E(0.0);
        }
        //Compute forces on each atom
        potential.calculate(allAtoms, forceSum);
        
    }//end of calculateForces
    
    protected void corrector() {
        
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.next();
            Agent agent = (IntegratorGear4.Agent)a.ia;
            Space.Vector r = a.coord.position();
            Space.Vector p = a.coord.momentum();
            work1.E(p);
            work1.PEa1Tv1(chi*a.coord.mass(),r);
            work2.E(work1);
            work2.ME(agent.dr1);
            r.PEa1Tv1(c0*a.coord.rm(), work2);
            agent.dr1.E(work1);
            agent.dr2.PEa1Tv1(c2,work2);
            agent.dr3.PEa1Tv1(c3,work2);
            agent.dr4.PEa1Tv1(c4,work2);
            
            work1.E(agent.force);
            work1.PEa1Tv1(-(zeta+chi),p);
            work2.E(work1);
            work2.ME(agent.dp1);
            p.PEa1Tv1(c0,work2);
            agent.dp1.E(work1);
            agent.dp2.PEa1Tv1(c2,work2);
            agent.dp3.PEa1Tv1(c3,work2);
            agent.dp4.PEa1Tv1(c4,work2);
        }
    }//end of corrector
        
    protected void predictor() {
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.next();
            Agent agent = (Agent)a.ia;
            Space.Vector r = a.coord.position();
            Space.Vector p = a.coord.momentum();
            double rm = a.coord.rm();
            r.PEa1Tv1(p1*rm, agent.dr1);
            r.PEa1Tv1(p2*rm, agent.dr2);
            r.PEa1Tv1(p3*rm, agent.dr3);
            r.PEa1Tv1(p4*rm, agent.dr4);
            
            agent.dr1.PEa1Tv1(p1, agent.dr2);
            agent.dr1.PEa1Tv1(p2, agent.dr3);
            agent.dr1.PEa1Tv1(p3, agent.dr4);
            
            agent.dr2.PEa1Tv1(p1, agent.dr3);
            agent.dr2.PEa1Tv1(p2, agent.dr4);
            
            agent.dr3.PEa1Tv1(p1, agent.dr4);
            
            p.PEa1Tv1(p1, agent.dp1);
            p.PEa1Tv1(p2, agent.dp2);
            p.PEa1Tv1(p3, agent.dp3);
            p.PEa1Tv1(p4, agent.dp4);
            
            agent.dp1.PEa1Tv1(p1, agent.dp2);
            agent.dp1.PEa1Tv1(p2, agent.dp3);
            agent.dp1.PEa1Tv1(p3, agent.dp4);
            
            agent.dp2.PEa1Tv1(p1, agent.dp3);
            agent.dp2.PEa1Tv1(p2, agent.dp4);
            
            agent.dp3.PEa1Tv1(p1, agent.dp4);
        }
    }
//--------------------------------------------------------------


    protected void doReset() {
        potential.set(firstPhase);//assumes only one phase
        calculateForces();
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.next();
            Agent agent = (IntegratorGear4.Agent)a.ia;
            agent.dr1.E(a.coord.momentum());
            agent.dr2.E(agent.force);
            agent.dr3.E(0.0);
            agent.dr4.E(0.0);
            agent.dp1.E(agent.force);
            agent.dp2.E(0.0);
            agent.dp3.E(0.0);
            agent.dp4.E(0.0);
        }
    }
              
//--------------------------------------------------------------

    public Integrator.Agent makeAgent(Atom a) {
        return new Agent(parentSimulation(),a);
    }
            
    public static class Agent implements Integrator.Agent.Forcible {  //need public so to use with instanceof
        public Atom atom;
        public Space.Vector force;
        public Space.Vector dr1, dr2, dr3, dr4;
        public Space.Vector dp1, dp2, dp3, dp4;

        public Agent(Simulation sim, Atom a) {
            atom = a;
            force = sim.space().makeVector();
            dr1 = sim.space().makeVector();
            dr2 = sim.space().makeVector();
            dr3 = sim.space().makeVector();
            dr4 = sim.space().makeVector();
            dp1 = sim.space().makeVector();
            dp2 = sim.space().makeVector();
            dp3 = sim.space().makeVector();
            dp4 = sim.space().makeVector();
        }
        
        public Space.Vector force() {return force;}
    }

    public static void main(String[] args) {
        
	    IntegratorGear4 integratorGear4 = new IntegratorGear4();
	    SpeciesSpheres speciesSpheres1 = new SpeciesSpheres();
//	    speciesSpheres1.setMass(1.0);
	    Phase phase1 = new Phase();
	    P2LennardJones P2LennardJones1 = new P2LennardJones();
	    Controller controller1 = new Controller();
	    DisplayPhase displayPhase1 = new DisplayPhase();
	    DisplayPlot plot = new DisplayPlot();
	    IntegratorMD.Timer timer = integratorGear4.new Timer(integratorGear4.chronoMeter());
	    timer.setUpdateInterval(10);
	//    integratorGear4.setTimeStep(0.005);
	    integratorGear4.setSleepPeriod(2);
		Simulation.instance.panel().setBackground(java.awt.Color.yellow);

        Meter ke = new MeterKineticEnergy();
        Meter temp = new MeterTemperature();
        Meter energy = new MeterEnergy();
        energy.setHistorying(true);
        plot.setDataSources(energy.getHistory());
        Phase phase = Simulation.instance.phase(0);
        ke.setPhase(phase);
        temp.setPhase(phase);
        energy.setPhase(phase);
        DisplayBox box1 = new DisplayBox();
        box1.setMeter(ke);
        box1.setUpdateInterval(10);
        DisplayBox box2 = new DisplayBox();
        box2.setMeter(temp);
        box2.setUnit(new Unit(Kelvin.UNIT));
        DisplayBox box3 = new DisplayBox();
        box3.setMeter(energy);
        box3.setPrecision(7);
                                            
		Simulation.instance.elementCoordinator.go();

   //     P2LennardJones1.setIterator(new AtomPairIterator(phase));
   //     P2LennardJones1.set(speciesSpheres1.getAgent(phase));
				
		Simulation.makeAndDisplayFrame(Simulation.instance);
    }//end of main
    
}//end of IntegratorGear4

