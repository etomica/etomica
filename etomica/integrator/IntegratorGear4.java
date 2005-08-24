// includes a main method

package etomica.integrator;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Integrator;
import etomica.IteratorDirective;
import etomica.Phase;
import etomica.PotentialMaster;
import etomica.Space;
import etomica.atom.Atom;
import etomica.atom.AtomTypeLeaf;
import etomica.potential.PotentialCalculationForceSum;
import etomica.space.ICoordinateKinetic;
import etomica.space.Vector;
import etomica.units.systems.LJ;

/**
 * Gear 4th-order predictor-corrector integrator.
 *
 * @author Ed Maginn
 * @author David Kofke
 */
public class IntegratorGear4 extends IntegratorMD implements EtomicaElement {

    private final PotentialCalculationForceSum forceSum;
    private final IteratorDirective allAtoms = new IteratorDirective();
    protected final Space space;
    final Vector work1, work2;
    double zeta = 0.0;
    double chi = 0.0;
    double p1, p2, p3, p4;
    double c0, c2, c3, c4;
    
    static final double GEAR0 = 251./720.;
    static final double GEAR2 = 11./12.;
    static final double GEAR3 = 1./3.;
    static final double GEAR4 = 1./24.;
                
    public IntegratorGear4(PotentialMaster potentialMaster, Space space) {
        super(potentialMaster);
        this.space = space;
        forceSum = new PotentialCalculationForceSum(space);
        work1 = space.makeVector();
        work2 = space.makeVector();
        //XXX this is totally wrong!  This should be based on the actual temperature and
        //potentials (steepness and depth) used.
        setTimeStep(new LJ().time().toSim(2.0));
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Molecular dynamics using Gear 4th-order predictor/corrector algorithm");
        return info;
    }

    public boolean addPhase(Phase p) {
        if(!super.addPhase(p)) return false;
        return true;
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
        
    
    public void doStep() {
        
        predictor();
        calculateForces();
        corrector();
        
    }//end of doStep
    
    protected void calculateForces() {
        
        //Compute all forces
        atomIterator.reset();
        while(atomIterator.hasNext()) {   //zero forces on all atoms
            ((IntegratorGear4.Agent)atomIterator.nextAtom().ia).force.E(0.0);
        }
        //Compute forces on each atom
        potential.calculate(firstPhase, allAtoms, forceSum);
        
    }//end of calculateForces
    
    protected void corrector() {
        
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.nextAtom();
            Agent agent = (IntegratorGear4.Agent)a.ia;
            Vector r = a.coord.position();
            Vector v = ((ICoordinateKinetic)a.coord).velocity();
            work1.E(v);
            work1.PEa1Tv1(chi,r);
            work2.E(work1);
            work2.ME(agent.dr1);
            r.PEa1Tv1(c0, work2);
            agent.dr1.E(work1);
            agent.dr2.PEa1Tv1(c2,work2);
            agent.dr3.PEa1Tv1(c3,work2);
            agent.dr4.PEa1Tv1(c4,work2);
            
            work1.Ea1Tv1(((AtomTypeLeaf)a.type).rm(),agent.force);
            work1.PEa1Tv1(-(zeta+chi),v);
            work2.E(work1);
            work2.ME(agent.dv1);
            v.PEa1Tv1(c0,work2);
            agent.dv1.E(work1);
            agent.dv2.PEa1Tv1(c2,work2);
            agent.dv3.PEa1Tv1(c3,work2);
            agent.dv4.PEa1Tv1(c4,work2);
        }
    }//end of corrector
        
    protected void predictor() {
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.nextAtom();
            Agent agent = (Agent)a.ia;
            Vector r = a.coord.position();
            Vector v = ((ICoordinateKinetic)a.coord).velocity();
            r.PEa1Tv1(p1, agent.dr1);
            r.PEa1Tv1(p2, agent.dr2);
            r.PEa1Tv1(p3, agent.dr3);
            r.PEa1Tv1(p4, agent.dr4);
            
            agent.dr1.PEa1Tv1(p1, agent.dr2);
            agent.dr1.PEa1Tv1(p2, agent.dr3);
            agent.dr1.PEa1Tv1(p3, agent.dr4);
            
            agent.dr2.PEa1Tv1(p1, agent.dr3);
            agent.dr2.PEa1Tv1(p2, agent.dr4);
            
            agent.dr3.PEa1Tv1(p1, agent.dr4);
            
            v.PEa1Tv1(p1, agent.dv1);
            v.PEa1Tv1(p2, agent.dv2);
            v.PEa1Tv1(p3, agent.dv3);
            v.PEa1Tv1(p4, agent.dv4);
            
            agent.dv1.PEa1Tv1(p1, agent.dv2);
            agent.dv1.PEa1Tv1(p2, agent.dv3);
            agent.dv1.PEa1Tv1(p3, agent.dv4);
            
            agent.dv2.PEa1Tv1(p1, agent.dv3);
            agent.dv2.PEa1Tv1(p2, agent.dv4);
            
            agent.dv3.PEa1Tv1(p1, agent.dv4);
        }
    }

    public void reset() {
        //XXX is this check really necessary?
        if(potential == null || firstPhase == null) return;
        calculateForces();
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.nextAtom();
            Agent agent = (IntegratorGear4.Agent)a.ia;
            agent.dr1.E(((ICoordinateKinetic)a.coord).velocity());
            agent.dr2.Ea1Tv1(((AtomTypeLeaf)a.type).rm(),agent.force);
            agent.dr3.E(0.0);
            agent.dr4.E(0.0);
            agent.dv1.Ea1Tv1(((AtomTypeLeaf)a.type).rm(),agent.force);
            agent.dv2.E(0.0);
            agent.dv3.E(0.0);
            agent.dv4.E(0.0);
        }
        super.reset();
    }
              
    public Object makeAgent(Atom a) {
        return new Agent(space,a);
    }
            
    public static class Agent implements Integrator.Forcible {  //need public so to use with instanceof
        public Atom atom;
        public Vector force;
        public Vector dr1, dr2, dr3, dr4;
        public Vector dv1, dv2, dv3, dv4;

        public Agent(Space space, Atom a) {
            atom = a;
            force = space.makeVector();
            dr1 = space.makeVector();
            dr2 = space.makeVector();
            dr3 = space.makeVector();
            dr4 = space.makeVector();
            dv1 = space.makeVector();
            dv2 = space.makeVector();
            dv3 = space.makeVector();
            dv4 = space.makeVector();
        }
        
        public Vector force() {return force;}
    }

/*    public static void main(String[] args) {
        
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
   */ 
}//end of IntegratorGear4

