package etomica.normalmode;

import etomica.action.BoxInflate;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.api.IVectorMutable;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.units.Dimension;
import etomica.units.Pressure;

/**
 * Standard Monte Carlo volume-change move for simulations in the NPT ensemble.
 *  for monoclinic primitive, where unit cell lengths: a, b and c will
 *  fluctuate independently to maintain a constant volume
 * 
 * 
 * @author Tai Boon Tan
 */
public class MCMoveVolumeMonoclinic extends MCMoveBoxStep {
    
    private static final long serialVersionUID = 2L;
    protected double pressure;
    private MeterPotentialEnergy energyMeter;
    protected final BoxInflate inflate;
    private IRandom random;
    protected final AtomIteratorLeafAtoms affectedAtomIterator;

    private transient double uOld, hOld, hNew;
    private transient double uNew = Double.NaN;

    public MCMoveVolumeMonoclinic(ISimulation sim, IPotentialMaster potentialMaster,
    		            ISpace _space) {
        this(potentialMaster, sim.getRandom(), _space, 1.0);
    }
    
    /**
     * @param potentialMaster an appropriate PotentialMaster instance for calculating energies
     * @param space the governing space for the simulation
     */
    public MCMoveVolumeMonoclinic(IPotentialMaster potentialMaster, IRandom random,
    		            ISpace _space, double pressure) {
        super(potentialMaster);
        this.random = random;
        inflate = new BoxInflate(_space);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.01);
        setPressure(pressure);
        energyMeter.setIncludeLrc(true);
        affectedAtomIterator = new AtomIteratorLeafAtoms();
    }
    
    public void setBox(IBox p) {
        super.setBox(p);
        energyMeter.setBox(p);
        inflate.setBox(p);
        affectedAtomIterator.setBox(p);
    }
     
    public boolean doTrial() {
        double vOld = box.getBoundary().volume();
        uOld = energyMeter.getDataAsScalar();
        hOld = uOld + pressure*vOld;
                
        IVectorMutable vecAOld = Space3D.makeVector(3);
        IVectorMutable vecBOld = Space3D.makeVector(3);
        IVectorMutable vecCOld = Space3D.makeVector(3);
        
        IVectorMutable vecANew = Space3D.makeVector(3);
        IVectorMutable vecCNew = Space3D.makeVector(3);
        
        vecAOld.E(box.getBoundary().getEdgeVector(0));
        vecBOld.E(box.getBoundary().getEdgeVector(1));
        vecCOld.E(box.getBoundary().getEdgeVector(2));
        
        /*
         * anisotropic scaling: scaleb and scalec
         * they are scale for b-VECTOR and c-VECTOR
         * 
         */
        double scalea = Math.exp((2.*random.nextDouble()-1.)*stepSize);
        double scalec = Math.exp((2.*random.nextDouble()-1.)*stepSize);
        
        vecANew.Ea1Tv1(scalea, vecAOld);
        vecCNew.Ea1Tv1(scalec, vecCOld);
        
        //System.out.println("\nbefore: " + box.getBoundary().volume());
        IVectorMutable rScale = Space3D.makeVector(3);
        rScale.E(new double[]{scalea, 1/(scalea*scalec), scalec});
        inflate.setVectorScale(rScale);
        inflate.actionPerformed();

        uNew = energyMeter.getDataAsScalar();
        hNew = uNew + pressure*vOld;
        return true;
    }//end of doTrial
    
    public double getA() {
        return Math.exp((box.getMoleculeList().getMoleculeCount()+1));
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