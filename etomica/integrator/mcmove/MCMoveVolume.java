package etomica.integrator.mcmove;

import etomica.action.BoxInflate;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.units.Dimension;
import etomica.units.Pressure;

/**
 * Standard Monte Carlo volume-change move for simulations in the NPT ensemble.
 *
 * @author David Kofke
 */
public class MCMoveVolume extends MCMoveBoxStep {
    
    private static final long serialVersionUID = 2L;
    protected double pressure;
    private MeterPotentialEnergy energyMeter;
    protected final BoxInflate inflate;
    private final int D;
    private IRandom random;

    private transient double uOld, hOld, vNew, vScale;
    private transient double uNew = Double.NaN;

    public MCMoveVolume(ISimulation sim, IPotentialMaster potentialMaster) {
        this(potentialMaster, sim.getRandom(), 1.0);
    }
    
    /**
     * @param potentialMaster an appropriate PotentialMaster instance for calculating energies
     * @param space the governing space for the simulation
     */
    public MCMoveVolume(IPotentialMaster potentialMaster, IRandom random, double pressure) {
        super(potentialMaster);
        this.random = random;
        this.D = potentialMaster.getSpace().D();
        inflate = new BoxInflate(potentialMaster.getSpace());
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        setPressure(pressure);
        energyMeter.setIncludeLrc(true);
    }
    
    public void setBox(IBox p) {
        super.setBox(p);
        energyMeter.setBox(p);
        inflate.setBox(p);
    }
    
    public boolean doTrial() {
        double vOld = box.volume();
        uOld = energyMeter.getDataAsScalar();
        hOld = uOld + pressure*vOld;
        vScale = (2.*random.nextDouble()-1.)*stepSize;
        vNew = vOld * Math.exp(vScale); //Step in ln(V)
        double rScale = Math.exp(vScale/D);
        inflate.setScale(rScale);
        inflate.actionPerformed();
        uNew = energyMeter.getDataAsScalar();
        hNew = uNew + pressure*vNew;
        return true;
    }//end of doTrial
    
    public double getA() {
        return Math.exp((box.moleculeCount()+1)*vScale);
    }
    
    public double getB() {
        return -(hNew - hOld);
    }
    
    public void acceptNotify() {  /* do nothing */}
    
    public void rejectNotify() {
        inflate.undo();
    }

    public double energyChange() {return uNew - uOld;}
    
    public AtomIterator affectedAtoms() {
        return new AtomIteratorAllMolecules(box);
    }

    public void setPressure(double p) {pressure = p;}
    public final double getPressure() {return pressure;}
    public final double pressure() {return pressure;}
    public Dimension getPressureDimension() {return Pressure.DIMENSION;}
    public final void setLogPressure(int lp) {setPressure(Math.pow(10.,lp));}
    private double hNew;
    
}