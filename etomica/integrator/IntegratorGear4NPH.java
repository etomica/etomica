// includes a main method

package etomica.integrator;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.action.PhaseInflate;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.iterator.AtomsetIterator;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataSource;
import etomica.data.meter.MeterTemperature;
import etomica.data.types.DataDoubleArray;
import etomica.exception.ConfigurationOverlapException;
import etomica.modifier.ModifierBoolean;
import etomica.phase.Phase;
import etomica.potential.Potential;
import etomica.potential.Potential2Soft;
import etomica.potential.PotentialCalculation;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.ICoordinateKinetic;
import etomica.space.NearestImageTransformer;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Kelvin;
import etomica.units.Null;
import etomica.units.Pressure;
import etomica.units.Temperature;
import etomica.units.Undefined;

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
    
    public IntegratorGear4NPH(Simulation sim) {
        this(sim.potentialMaster,sim.space,sim.getDefaults().timeStep,
                sim.getDefaults().temperature);
    }
    
    public IntegratorGear4NPH(PotentialMaster potentialMaster, Space space, 
            double timeStep, double temperature) {
        super(potentialMaster, space, timeStep, temperature);
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
    public Dimension getRelaxationRatePDimension() {return Null.DIMENSION;}
    
    public void setRelaxationRateH(double value) {
        rrh = value;
        setTimeStep(timeStep);
    }
    public double getRelaxationRateH() {return rrh;}
    public Dimension getRelaxationRateHDimension() {return Null.DIMENSION;}
    
    public synchronized void setTargetH(double value) {targetH = value;}
    public double getTargetH() {return targetH;}
    public Dimension getTargetHDimension() {return Energy.DIMENSION;}
    
    public synchronized void setTargetP(double value) {targetP = value;}
    public double getTargetP() {return targetP;}
    public Dimension getTargetPDimension() {return Pressure.DIMENSION;}
    
    public synchronized void setTargetT(double value) {targetT = value;}
    public double getTargetT() {return targetT;}
    public Dimension getTargetTDimension() {return Temperature.DIMENSION;}

    public void setPhase(Phase p) {
        super.setPhase(p);
        inflate.setPhase(phase);
        meterTemperature.setPhase(phase);
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
        double mvsq = kineticT * D * phase.atomCount();
        double volume = phase.volume();
        double pCurrent = phase.getDensity()*kineticT - forceSumNPH.w/(D*volume);
        double pDot = kp*(targetP - pCurrent);
        double kDot = kp*(targetT - kineticT)*phase.atomCount();
        chi = ( - forceSumNPH.rvx - D*pDot*volume)/
                    ( forceSumNPH.x + D*D*pCurrent*volume);
        zeta = (forceSumNPH.vf - kDot)/mvsq - chi;
    }
    
    public void drivePH() {
        double kineticT = meterTemperature.getDataAsScalar();
        double mvsq = kineticT * D * phase.atomCount();
        double volume = phase.volume();
        double pCurrent = phase.getDensity()*kineticT - forceSumNPH.w/(D*volume);
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
            agents[atomIterator.nextAtom().getGlobalIndex()].force.E(0.0);
        }
        forceSumNPH.u = 0.0;
        forceSumNPH.w = 0.0;
        forceSumNPH.x = 0.0;
        forceSumNPH.rvx = 0.0;
        forceSumNPH.vf = 0.0;

        potential.calculate(phase, allAtoms, forceSumNPH);
    }//end of calculateForces
    
    protected void corrector() {
        super.corrector();
        double volOld = phase.getBoundary().volume();
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
        double volOld = phase.getBoundary().volume();
        double volNew = volOld + p1*vol1 + p2*vol2 + p3*vol3 + p4*vol4;
        vol1 += p1*vol2 + p2*vol3 + p3*vol4;
        vol2 += p1*vol3 + p2*vol4;
        vol3 += p1*vol4;
        double rScale = Math.pow(volNew/volOld,1.0/D);
        inflate.setScale(rScale);
        inflate.actionPerformed();
    }
//--------------------------------------------------------------


    public void reset() throws ConfigurationOverlapException {
        vol1 = 0.0;
        vol2 = 0.0;
        vol3 = 0.0;
        vol4 = 0.0;
        super.reset();
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
                Phase phase = integrator.getPhase();
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
            data = new DataDoubleArray("TPH",Undefined.DIMENSION,3);
            this.integrator = integrator;
        }
        
        public DataInfo getDataInfo() {
            return data.getDataInfo();
        }
        
        public Data getData() {
            Phase phase = integrator.getPhase();
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
        private final Vector dr;
        private final Vector dv;
        private NearestImageTransformer nearestImageTransformer;

        private final Vector f;
        public ForceSumNPH(Space space) {
            f = space.makeVector();
            dr = space.makeVector();
            dv = space.makeVector();
        }
        
        public void setPhase(Phase phase) {
            nearestImageTransformer = phase.getBoundary();
        }
        
        //pair
        public void doCalculation(AtomsetIterator iterator, Potential potential2) {
            Potential2Soft potentialSoft = (Potential2Soft)potential2;
            while(iterator.hasNext()) {
                AtomPair pair = (AtomPair)iterator.next();

                AtomLeaf atom0 = (AtomLeaf)pair.atom0;
                AtomLeaf atom1 = (AtomLeaf)pair.atom1;
                ICoordinateKinetic coord0 = (ICoordinateKinetic)atom0.coord;
                ICoordinateKinetic coord1 = (ICoordinateKinetic)atom1.coord;
                dv.Ev1Mv2(coord1.velocity(), coord0.velocity());
                
                dr.Ev1Mv2(coord1.position(), coord0.position());
                nearestImageTransformer.nearestImage(dr);

                double r2 = dr.squared();
                u += potentialSoft.energy(pair);
                w += potentialSoft.virial(pair);
                double hv = potentialSoft.hyperVirial(pair);
                x += hv;
                rvx += hv * dr.dot(dv)/r2;
                f.E(potentialSoft.gradient(pair));
                vf -= dv.dot(f); //maybe should be (-)?
                agents[pair.atom0.getGlobalIndex()].force().PE(f);
                agents[pair.atom1.getGlobalIndex()].force().ME(f);
            }
        }
    }
}
