// includes a main method

package etomica.integrator;
import etomica.Atom;
import etomica.AtomPair;
import etomica.AtomsetIterator;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Integrator;
import etomica.IteratorDirective;
import etomica.Phase;
import etomica.Potential;
import etomica.PotentialMaster;
import etomica.Space;
import etomica.action.PhaseInflate;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataSource;
import etomica.data.meter.MeterTemperature;
import etomica.data.types.DataDoubleArray;
import etomica.modifier.ModifierBoolean;
import etomica.potential.Potential2Soft;
import etomica.potential.PotentialCalculation;
import etomica.space.CoordinatePairKinetic;
import etomica.space.Vector;
import etomica.units.Dimension;
import etomica.units.Kelvin;

/**
 * Gear 4th-order predictor-corrector integrator for constant enthalphy, pressure.
 *
 * @author Ed Maginn
 * @author David Kofke
 */
public final class IntegratorGear4NPH extends IntegratorGear4 implements EtomicaElement {

    double vol1, vol2, vol3, vol4;
    protected /*final*/ ForceSumNPH forceSumNPH;//MeterTPH won't permit this to be final (?)
    private final IteratorDirective allAtoms = new IteratorDirective();
    protected final PhaseInflate inflate;
    double targetH;
    double targetP;
    double targetT = Kelvin.UNIT.toSim(300.);
    private double rrp = 300.;
    private double rrh = 300.;
    private double kp, kh;
    protected int D;
    protected final MeterTemperature meterTemperature = new MeterTemperature();
    
    public IntegratorGear4NPH(PotentialMaster potentialMaster, Space space) {
        super(potentialMaster, space);
        kp = 1.0/rrp/getTimeStep();
        kh = 1.0/rrh/getTimeStep();
        D = space.D();
        setIsothermal(true);
        forceSumNPH = new ForceSumNPH(space);
        inflate = new PhaseInflate(space);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Gear4 MD at constant pressure and enthalpy");
        return info;
    }
    
    public void setTimeStep(double t) {
        super.setTimeStep(t);
        kp = 1.0/(rrp*t);
        kh = 1.0/(rrh*t);
    }
    
    public void setRelaxationRateP(double value) {
        rrp = value;
        setTimeStep(timeStep);
    }
    public double getRelaxationRateP() {return rrp;}
    public Dimension getRelaxationRatePDimension() {return Dimension.NULL;}
    
    public void setRelaxationRateH(double value) {
        rrh = value;
        setTimeStep(timeStep);
    }
    public double getRelaxationRateH() {return rrh;}
    public Dimension getRelaxationRateHDimension() {return Dimension.NULL;}
    
    public synchronized void setTargetH(double value) {targetH = value;}
    public double getTargetH() {return targetH;}
    public Dimension getTargetHDimension() {return Dimension.ENERGY;}
    
    public synchronized void setTargetP(double value) {targetP = value;}
    public double getTargetP() {return targetP;}
    public Dimension getTargetPDimension() {return Dimension.PRESSURE;}
    
    public synchronized void setTargetT(double value) {targetT = value;}
    public double getTargetT() {return targetT;}
    public Dimension getTargetTDimension() {return Dimension.TEMPERATURE;}

    public boolean addPhase(Phase p) {
        boolean b = super.addPhase(p);
        inflate.setPhase(firstPhase);
        meterTemperature.setPhase(firstPhase);
        return b;
    }
    
    
//--------------------------------------------------------------
// steps all particles across time interval tStep

    public void doStep() {
        
        predictor();
        calculateForces();
        if(isothermal) drivePT();
        else drivePH();
        corrector();
    }//end of doStep
    
    public void drivePT() {
        double kineticT = meterTemperature.getDataAsScalar();
        double mvsq = kineticT * D * firstPhase.atomCount();
        double volume = firstPhase.volume();
        double pCurrent = firstPhase.getDensity()*kineticT - forceSumNPH.w/(D*volume);
        double pDot = kp*(targetP - pCurrent);
        double kDot = kp*(targetT - kineticT)*firstPhase.atomCount();
        chi = ( - forceSumNPH.rvx - D*pDot*volume)/
                    ( forceSumNPH.x + D*D*pCurrent*volume);
        zeta = (forceSumNPH.vf - kDot)/mvsq - chi;
    }
    
    public void drivePH() {
        double kineticT = meterTemperature.getDataAsScalar();
        double mvsq = kineticT * D * firstPhase.atomCount();
        double volume = firstPhase.volume();
        double pCurrent = firstPhase.getDensity()*kineticT - forceSumNPH.w/(D*volume);
        double hCurrent = 0.5*mvsq + forceSumNPH.u + pCurrent*volume;
        double pDot = kp*(targetP - pCurrent);
        double hDot = kh*(targetH - hCurrent);
        zeta = (pDot*volume - hDot)/mvsq;
        chi = (2.0*forceSumNPH.vf -forceSumNPH.rvx - D*pDot*volume - 2*zeta*mvsq)/
                    (2.0*mvsq + forceSumNPH.x + D*D*pCurrent*volume);
    }
    
    protected void calculateForces() {
        
        //Compute all forces
        atomIterator.reset();
        while(atomIterator.hasNext()) {   //zero forces on all atoms
            ((Agent)atomIterator.nextAtom().ia).force.E(0.0);
        }
        forceSumNPH.u = 0.0;
        forceSumNPH.w = 0.0;
        forceSumNPH.x = 0.0;
        forceSumNPH.rvx = 0.0;
        forceSumNPH.vf = 0.0;

        potential.calculate(firstPhase, allAtoms, forceSumNPH);
    }//end of calculateForces
    
    protected void corrector() {
        super.corrector();
        double volOld = firstPhase.boundary().volume();
        double voi = D*volOld*chi;
        double corvol = voi - vol1;
        double volNew = volOld + c0*corvol;
        vol1 = voi;
        vol2 += c2*corvol;
        vol3 += c3*corvol;
        vol4 += c4*corvol;
        double rScale = Math.pow(volNew/volOld,1.0/D);
        inflate.setScale(rScale);
        inflate.actionPerformed();
    }//end of corrector
        
    protected void predictor() {
        super.predictor();
        double volOld = firstPhase.boundary().volume();
        double volNew = volOld + p1*vol1 + p2*vol2 + p3*vol3 + p4*vol4;
        vol1 += p1*vol2 + p2*vol3 + p3*vol4;
        vol2 += p1*vol3 + p2*vol4;
        vol3 += p1*vol4;
        double rScale = Math.pow(volNew/volOld,1.0/D);
        inflate.setScale(rScale);
        inflate.actionPerformed();
    }
//--------------------------------------------------------------


    public void reset() {
        super.reset();
        vol1 = 0.0;
        vol2 = 0.0;
        vol3 = 0.0;
        vol4 = 0.0;
    }
              
//--------------------------------------------------------------

    public Object makeAgent(Atom a) {
        return new Agent(space,a);
    }
            
    public static class Agent extends IntegratorGear4.Agent {  //need public so to use with instanceof
        public Agent(Space space, Atom a) {
            super(space, a);
        }
    }
    
    //inner class used to toggle between NPT and NPH ensembles
    public static class EnsembleToggler implements ModifierBoolean, java.io.Serializable {
        public EnsembleToggler(IntegratorGear4NPH integrator) {
            this.integrator = integrator;
        }
        public void setBoolean(boolean isothermal) {
            integrator.setIsothermal(isothermal);
            if(!isothermal) {
                integrator.calculateForces();
                Phase phase = integrator.getPhase()[0];
                double kineticT = integrator.getMeterTemperature().getDataAsScalar();
                double mvsq = kineticT * phase.space().D() * phase.atomCount();
                double volume = phase.volume();
                double pCurrent = phase.getDensity()*kineticT - integrator.forceSumNPH.w/(phase.space().D()*volume);
                double hCurrent = 0.5*mvsq + integrator.forceSumNPH.u + pCurrent*volume;
                integrator.setTargetH(hCurrent);
            }
        }
        public boolean getBoolean() {return integrator.isIsothermal();}
        private IntegratorGear4NPH integrator;
    }

    //meter group for temperature, pressure, enthalpy, obtaining values from
    //most recent call to the ForceSumNPH instance
    public static final class MeterTPH implements DataSource, java.io.Serializable {
        
        public MeterTPH(IntegratorGear4NPH integrator) {
            data = new DataDoubleArray("TPH",Dimension.UNDEFINED,3);
            this.integrator = integrator;
        }
        
        public DataInfo getDataInfo() {
            return data.getDataInfo();
        }
        
        public Data getData() {
            Phase phase = integrator.getPhase()[0];
            double kineticT = integrator.getMeterTemperature().getDataAsScalar();
            double mvsq = kineticT* phase.space().D() * phase.atomCount();
            double volume = phase.volume();
            double pCurrent = phase.getDensity()*kineticT - integrator.forceSumNPH.w/(phase.space().D()*volume);
            double hCurrent = 0.5*mvsq + integrator.forceSumNPH.u + pCurrent*volume;
            data.getData()[0] = kineticT;
            data.getData()[1] = pCurrent;
            data.getData()[2] = hCurrent/phase.moleculeCount();
            return data;
        }
        
        private DataDoubleArray data;
        private IntegratorGear4NPH integrator;
    }
        
    public final class ForceSumNPH extends PotentialCalculation {
        
        double u; //energy sum
        double w; //virial sum
        double x; //hypervirial sum
        double rvx; 
        double vf;
        private CoordinatePairKinetic cPair;

        private final Vector f;
        public ForceSumNPH(Space space) {
            f = space.makeVector();
            cPair = new CoordinatePairKinetic(space);
        }
        
        //pair
        public void doCalculation(AtomsetIterator iterator, Potential potential2) {
            Potential2Soft potentialSoft = (Potential2Soft)potential2;
            while(iterator.hasNext()) {
                AtomPair pair = (AtomPair)iterator.next();
                cPair.reset(pair);
                double r2 = cPair.r2();
                u += potentialSoft.energy(pair);
                w += potentialSoft.virial(pair);
                double hv = potentialSoft.hyperVirial(pair);
                x += hv;
                rvx += hv * cPair.vDotr()/r2;
                f.E(potentialSoft.gradient(pair));
                vf -= cPair.dv().dot(f); //maybe should be (-)?
                ((Integrator.Forcible)pair.atom0.ia).force().PE(f);
                ((Integrator.Forcible)pair.atom1.ia).force().ME(f);
            }//end while
        }//end of calculate
    }//end ForceSums

/*    public static void main(String[] args) {
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
    */
}

