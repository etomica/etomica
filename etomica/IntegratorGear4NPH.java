// includes a main method

package etomica;
import etomica.units.*;

/**
 * Gear 4th-order predictor-corrector integrator for constant enthalphy, pressure.
 *
 * @author Ed Maginn
 * @author David Kofke
 */
public final class IntegratorGear4NPH extends IntegratorGear4 implements EtomicaElement {

    public String getVersion() {return "IntegratorGear4NPH:01.10.19/"+IntegratorGear4.VERSION;}

    double vol1, vol2, vol3, vol4;
    protected final ForceSumNPH forceSum;
    private final IteratorDirective allAtoms = new IteratorDirective();
    protected PhaseAction.Inflate inflate;
    double targetH;
    double targetP;
    private double kineticT;
    private double rrp = 300.;
    private double rrh = 300.;
    private double kp, kh;
    private int D;
    
    public IntegratorGear4NPH() {
        this(Simulation.instance);
    }
    public IntegratorGear4NPH(final Simulation sim) {
        super(sim);
        kp = 1.0/rrp/timeStep();
        kh = 1.0/rrp/timeStep();
        D = sim.space().D();
        setIsothermal(true);
        forceSum = new ForceSumNPH(sim.space());
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Gear4 MD at constant pressure and enthalpy");
        return info;
    }
    
    public void setTimeStep(double t) {
        super.setTimeStep(t);
        kp = 1.0/(rrp*t);
        kh = 1.0/(rrp*t);
    }
    
    public void setRelaxationRate(double value) {
        rrp = value;
        setTimeStep(timeStep);
    }
    public double getRelaxationRate() {return rrp;}
    public Dimension getRelaxationRateDimension() {return Dimension.NULL;}
    
    public void setTargetH(double value) {targetH = value;}
    public double getTargetH() {return targetH;}
    
    public void setTargetP(double value) {targetP = value;}
    public double getTargetP() {return targetP;}
    public Dimension getTargetPDimension() {return Dimension.PRESSURE;}

    public boolean addPhase(Phase p) {
        boolean b = super.addPhase(p);
        inflate = new PhaseAction.Inflate(firstPhase);
        return b;
    }
    
    
//--------------------------------------------------------------
// steps all particles across time interval tStep

    public void doStep() {
        
        predictor();
        calculateForces();
        if(isothermal) driveP();
        else drivePH();
        corrector();
        if(isothermal) scaleMomenta(Math.sqrt(this.temperature/kineticT));
        
    }//end of doStep
    
    public void driveP() {
        kineticT = firstPhase.kineticTemperature();
        double volume = firstPhase.volume();
        double pCurrent = firstPhase.getDensity()*kineticT - forceSum.w/(D*volume);
        double dpdt = kp*(targetP - pCurrent);
        chi = (2.0*forceSum.vf - forceSum.rvx - D*dpdt*volume)/
                    (2.*kineticT*D*firstPhase.atomCount() + forceSum.x + D*D*pCurrent*volume);
    }
    
    public void drivePH() {
        
    }
    
    public void scaleMomenta(double s) {
        double rs = 1.0/s;
        for(Atom a=firstPhase.firstAtom(); a!=null; a=a.nextAtom()) {
            a.coord.momentum().TE(s); //scale momentum
        }
    }
    private void calculateForces() {
        
        //Compute all forces
        atomIterator.reset();
        while(atomIterator.hasNext()) {   //zero forces on all atoms
            ((Agent)atomIterator.next().ia).force.E(0.0);
        }
        forceSum.w = 0.0;
        forceSum.x = 0.0;
        forceSum.rvx = 0.0;
        forceSum.vf = 0.0;

        potential.calculate(allAtoms,forceSum);
    }//end of calculateForces
    
    protected void corrector() {
        super.corrector();
        double volOld = firstPhase.boundary().volume();
        double voi = 3.0*volOld*chi;
        double corvol = voi - vol1;
        double volNew = volOld + c0*corvol;
        vol1 = voi;
        vol2 += c2*corvol;
        vol3 += c3*corvol;
        vol4 += c4*corvol;
        double rScale = Math.pow(volNew/volOld,1.0/(double)parentSimulation().space().D());
        inflate.actionPerformed(rScale);
    }//end of corrector
        
    protected void predictor() {
        super.predictor();
        double volOld = firstPhase.boundary().volume();
        double volNew = volOld + p1*vol1 + p2*vol2 + p3*vol3 + p4*vol4;
        vol1 += p1*vol2 + p2*vol3 + p3*vol4;
        vol2 += p1*vol3 + p2*vol4;
        vol3 += p1*vol4;
        double rScale = Math.pow(volNew/volOld,1.0/(double)parentSimulation().space().D());
        inflate.actionPerformed(rScale);
    }
//--------------------------------------------------------------


    protected void doReset() {
        super.doReset();
    }
              
//--------------------------------------------------------------

    public Integrator.Agent makeAgent(Atom a) {
        return new Agent(parentSimulation(),a);
    }
            
    public static class Agent extends IntegratorGear4.Agent {  //need public so to use with instanceof
        public Agent(Simulation sim, Atom a) {
            super(sim, a);
        }
    }
    
    public final class ForceSumNPH implements Potential1Calculation, Potential2Calculation {
            
        double w; //virial sum
        double x; //hypervirial sum
        double rvx; 
        double vf;

        private final Space.Vector f;
        public ForceSumNPH(Space space) {
            f = space.makeVector();
        }
            
        //atom
        public void calculate(AtomIterator iterator, Potential1 potential) {
            Potential1Soft potentialSoft = (Potential1Soft)potential;
            while(iterator.hasNext()) {
                Atom atom = iterator.next();
                f.E(potentialSoft.gradient(atom));
                ((Integrator.Agent.Forcible)atom.ia).force().ME(f);
            }//end while
        }//end of calculate

        //pair
        public void calculate(AtomPairIterator iterator, Potential2 potential) {
            Potential2Soft potentialSoft = (Potential2Soft)potential;
            while(iterator.hasNext()) {
                AtomPair pair = iterator.next();
                double r2 = pair.r2();
                w += potentialSoft.virial(pair);
                double hv = potentialSoft.hyperVirial(pair);
                x += hv;
                rvx += hv * pair.vDotr()/r2;
                f.E(potentialSoft.gradient(pair));
                vf += pair.cPair.vDot(f);
                ((Integrator.Agent.Forcible)pair.atom1().ia).force().PE(f);
                ((Integrator.Agent.Forcible)pair.atom2().ia).force().ME(f);
            }//end while
        }//end of calculate
    }//end ForceSums
        

    public static void main(String[] args) {
	    IntegratorGear4 integratorGear4 = new IntegratorGear4();
	    SpeciesSpheres speciesDisks1 = new SpeciesSpheres();
//	    speciesDisks1.setMass(1.0);
	    Phase phase1 = new Phase();
	    P2LennardJones P2LennardJones1 = new P2LennardJones();
	    Controller controller1 = new Controller();
	    DisplayPhase displayPhase1 = new DisplayPhase();
	    IntegratorMD.Timer timer = integratorGear4.new Timer(integratorGear4.chronoMeter());
	    timer.setUpdateInterval(10);
	    integratorGear4.setTimeStep(0.005);
		Simulation.instance.panel().setBackground(java.awt.Color.yellow);

        Meter ke = new MeterKineticEnergy();
        Meter temp = new MeterTemperature();
        Meter energy = new MeterEnergy();
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
                                            
		Simulation.instance.elementCoordinator.go(); //invoke this method only after all elements are in place
		                                    //calling it a second time has no effect
		
	    Simulation.makeAndDisplayFrame(Simulation.instance);
    }//end of main
}

