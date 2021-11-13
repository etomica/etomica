/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.BoxInflateAnisotropic;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Standard Monte Carlo volume-change move for simulations in the NPT ensemble.
 *  for monoclinic primitive, where the x-component of c-Vector is varied to sample
 *  box with angle fluctuation between a-Vecvtor and c-Vector and also 
 *  maintaining a constant volume
 * 
 * @author Tai Boon Tan
 */
public class MCMoveVolumeMonoclinicAngle extends MCMoveBoxStep {

    protected final PotentialCompute potentialCompute;
    protected final IntegratorMC integrator;
    protected final BoxInflateAnisotropic inflate;
    private final IRandom random;
    protected final AtomIteratorLeafAtoms affectedAtomIterator;

    private transient double uOld;
    private transient double uNew = Double.NaN;
    protected final Vector cVec;

    public MCMoveVolumeMonoclinicAngle(PotentialCompute potentialCompute, IntegratorMC integrator, IRandom random) {
        super();
        this.random = random;
        this.potentialCompute = potentialCompute;
        this.integrator = integrator;
        inflate = new BoxInflateAnisotropic(integrator.getBox(), integrator.getBox().getSpace());
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.01);
        affectedAtomIterator = new AtomIteratorLeafAtoms();
        cVec = integrator.getBox().getSpace().makeVector();
        setBox(integrator.getBox());
    }
    
    public void setBox(Box p) {
        super.setBox(p);
        inflate.setBox(p);
        affectedAtomIterator.setBox(p);
    }
     
    public boolean doTrial() {
        uOld = integrator.getPotentialEnergy();
                
        cVec.E(box.getBoundary().getEdgeVector(2));
        
        double scale = Math.exp((2.*random.nextDouble()-1.)*stepSize);
        cVec.setX(0, scale*cVec.getX(0));

        inflate.setCVector(cVec);
        inflate.actionPerformed();

        potentialCompute.init();
        uNew = potentialCompute.computeAll(false);
        return true;
    }//end of doTrial

    public double getChi(double temperature) {
        return Math.exp(-(uNew - uOld) / temperature);
    }
    
    public void acceptNotify() {
    }
    
    public void rejectNotify() {
        inflate.undo();
        potentialCompute.init();
        potentialCompute.computeAll(false);
    }

    public double energyChange() {return uNew - uOld;}
    
    public AtomIterator affectedAtoms() {
        return affectedAtomIterator;
    }
}
