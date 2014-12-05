/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
    protected final IVectorMutable scaleVector;
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
        scaleVector = _space.makeVector();
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
                
        /*
         * anisotropic scaling: scaleb and scalec
         * they are scale for b-VECTOR and c-VECTOR
         * 
         */
        double scalea = Math.exp((2.*random.nextDouble()-1.)*stepSize);
        double scalec = Math.exp((2.*random.nextDouble()-1.)*stepSize);
        
        //System.out.println("\nbefore: " + box.getBoundary().volume());
        scaleVector.setX(0, scalea);
        scaleVector.setX(1, 1/(scalea*scalec));
        scaleVector.setX(2, scalec);
        
        inflate.setVectorScale(scaleVector);
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