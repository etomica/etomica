package etomica.integrator.mcmove;

import etomica.action.PhaseInflate;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.units.Dimension;
import etomica.units.Pressure;
import etomica.util.IRandom;

/**
 * Standard Monte Carlo volume-change move for simulations in the NPT ensemble.
 *
 * @author David Kofke
 */
public class MCMoveVolume extends MCMovePhaseStep {
    
    private static final long serialVersionUID = 1L;
    protected double pressure;
    private MeterPotentialEnergy energyMeter;
    protected final PhaseInflate inflate;
    private final int D;
    private IRandom random;

    private transient double uOld, hOld, vNew, vScale;
    private transient double uNew = Double.NaN;

    public MCMoveVolume(Simulation sim) {
        this(sim.getPotentialMaster(), sim.getRandom(), sim.getDefaults().pressure);
    }
    
    /**
     * @param potentialMaster an appropriate PotentialMaster instance for calculating energies
     * @param space the governing space for the simulation
     */
    public MCMoveVolume(PotentialMaster potentialMaster, IRandom random, double pressure) {
        super(potentialMaster);
        this.random = random;
        this.D = potentialMaster.getSpace().D();
        inflate = new PhaseInflate(potentialMaster.getSpace());
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        setPressure(pressure);
        energyMeter.setIncludeLrc(true);
        setName("MCMoveVolume");
        
    }
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        energyMeter.setPhase(p);
        inflate.setPhase(p);
    }
    
    public boolean doTrial() {
        double vOld = phase.volume();
        uOld = energyMeter.getDataAsScalar();
        hOld = uOld + pressure*vOld;
        vScale = (2.*random.nextDouble()-1.)*stepSize;
        vNew = vOld * Math.exp(vScale); //Step in ln(V)
        double rScale = Math.exp(vScale/D);
        inflate.setScale(rScale);
        inflate.actionPerformed();
        uNew = Double.NaN;
        return true;
    }//end of doTrial
    
    public double getA() {
        return Math.exp((phase.moleculeCount()+1)*vScale);
    }
    
    public double getB() {
        uNew = energyMeter.getDataAsScalar();
        double hNew = uNew + pressure*vNew;
        return -(hNew - hOld);
    }
    
    public void acceptNotify() {  /* do nothing */}
    
    public void rejectNotify() {
        inflate.undo();
    }

    public double energyChange() {return uNew - uOld;}
    
    public AtomIterator affectedAtoms() {
        return new AtomIteratorAllMolecules(phase);
    }

    public void setPressure(double p) {pressure = p;}
    public final double getPressure() {return pressure;}
    public final double pressure() {return pressure;}
    public Dimension getPressureDimension() {return Pressure.DIMENSION;}
    public final void setLogPressure(int lp) {setPressure(Math.pow(10.,lp));}
    
}