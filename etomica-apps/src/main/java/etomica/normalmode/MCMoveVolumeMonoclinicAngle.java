/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.BoxInflateAnisotropic;
import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.util.random.IRandom;
import etomica.space.Vector;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.space.Space;

/**
 * Standard Monte Carlo volume-change move for simulations in the NPT ensemble.
 *  for monoclinic primitive, where the x-component of c-Vector is varied to sample
 *  box with angle fluctuation between a-Vecvtor and c-Vector and also 
 *  maintaining a constant volume
 * 
 * @author Tai Boon Tan
 */
public class MCMoveVolumeMonoclinicAngle extends MCMoveBoxStep {
    
    private static final long serialVersionUID = 2L;
    private MeterPotentialEnergy energyMeter;
    protected final BoxInflateAnisotropic inflate;
    private IRandom random;
    protected final AtomIteratorLeafAtoms affectedAtomIterator;

    private transient double uOld;
    private transient double uNew = Double.NaN;
    protected final Vector cVec;

    /**
     * @param potentialMaster an appropriate PotentialMaster instance for calculating energies
     * @param space the governing space for the simulation
     */
    public MCMoveVolumeMonoclinicAngle(PotentialMaster potentialMaster, IRandom random,
                                       Space _space, Box box) {
        super(potentialMaster);
        this.random = random;
        inflate = new BoxInflateAnisotropic(box, _space);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.01);
        energyMeter.setIncludeLrc(true);
        affectedAtomIterator = new AtomIteratorLeafAtoms();
        cVec = _space.makeVector();
    }
    
    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
        inflate.setBox(p);
        affectedAtomIterator.setBox(p);
    }
     
    public boolean doTrial() {
        uOld = energyMeter.getDataAsScalar();
                
        cVec.E(box.getBoundary().getEdgeVector(2));
        
        double scale = Math.exp((2.*random.nextDouble()-1.)*stepSize);
        cVec.setX(0, scale*cVec.getX(0));

        inflate.setCVector(cVec);
        inflate.actionPerformed();

        uNew = energyMeter.getDataAsScalar();
        return true;
    }//end of doTrial
    
    public double getA() {
        return 1;
    }
    
    public double getB() {
        return -(uNew - uOld);
    }
    
    public void acceptNotify() {
    }
    
    public void rejectNotify() {
        inflate.undo();
    }

    public double energyChange() {return uNew - uOld;}
    
    public AtomIterator affectedAtoms() {
        return affectedAtomIterator;
    }
}
