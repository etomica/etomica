package etomica.integrator.mcmove;

import etomica.action.BoxInflate;
import etomica.api.IBox;
import etomica.api.IFunction;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.space.ISpace;
import etomica.units.Dimension;
import etomica.units.Pressure;
import etomica.util.Function;

/**
 * Standard Monte Carlo volume-change move for simulations in the NPT ensemble.
 *
 * @author David Kofke
 */
public class MCMoveVolume extends MCMoveBoxStep {

    protected double pressure;
    private MeterPotentialEnergy energyMeter;
    protected BoxInflate inflate;
    private final int D;
    private IRandom random;
    protected final AtomIteratorLeafAtoms affectedAtomIterator;
    protected IFunction vBias;

    private transient double biasOld, uOld, hOld, vNew, vScale, hNew;
    private transient double uNew = Double.NaN;

    public MCMoveVolume(ISimulation sim, IPotentialMaster potentialMaster,
    		            ISpace _space) {
        this(potentialMaster, sim.getRandom(), _space, 1.0);
    }

    /**
     * @param potentialMaster an appropriate PotentialMaster instance for calculating energies
     * @param space the governing space for the simulation
     */
    public MCMoveVolume(IPotentialMaster potentialMaster, IRandom random,
    		            ISpace _space, double pressure) {
        super(potentialMaster);
        this.random = random;
        this.D = _space.D();
        inflate = new BoxInflate(_space);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        setPressure(pressure);
        energyMeter.setIncludeLrc(true);
        affectedAtomIterator = new AtomIteratorLeafAtoms();
        vBias = new Function.Constant(1);
    }

    public void setInflater(BoxInflate newInflate) {
        inflate = newInflate;
        inflate.setBox(box);
    }

    public void setBox(IBox p) {
        super.setBox(p);
        energyMeter.setBox(p);
        inflate.setBox(p);
        affectedAtomIterator.setBox(p);
    }

    public void setVolumeBias(IFunction vBias) {
        this.vBias = vBias;
    }

    public boolean doTrial() {
        double vOld = box.getBoundary().volume();
        uOld = energyMeter.getDataAsScalar();
        hOld = uOld + pressure*vOld;
        biasOld = vBias.f(vOld);
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
        // N, not N+1 here because of the shell volume
        // D. S. Corti, Mol. Phys. 100, 1887 (2002).
        double biasNew = vBias.f(box.getBoundary().volume());
        return biasNew/biasOld*Math.exp(box.getMoleculeList().getMoleculeCount()*vScale);
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
        return affectedAtomIterator;
    }

    public void setPressure(double p) {pressure = p;}
    public final double getPressure() {return pressure;}
    public Dimension getPressureDimension() {return Pressure.DIMENSION;}
}