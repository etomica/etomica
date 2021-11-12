/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.BoxInflate;
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
 *  for monoclinic primitive, where unit cell lengths: a, b and c will
 *  fluctuate independently to maintain a constant volume
 * 
 * 
 * @author Tai Boon Tan
 */
public class MCMoveVolumeMonoclinic extends MCMoveBoxStep {

    protected final PotentialCompute potentialCompute;
    protected final IntegratorMC integrator;
    protected BoxInflate inflate;
    private final IRandom random;
    protected final Vector scaleVector;
    protected final AtomIteratorLeafAtoms affectedAtomIterator;

    private transient double uOld;
    private transient double uNew = Double.NaN;

    public MCMoveVolumeMonoclinic(PotentialCompute potentialCompute, IntegratorMC integrator, IRandom random) {
        super();
        this.potentialCompute = potentialCompute;
        this.integrator = integrator;
        this.random = random;
        inflate = new BoxInflate(integrator.getBox().getSpace());
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.01);
        affectedAtomIterator = new AtomIteratorLeafAtoms();
        scaleVector = integrator.getBox().getSpace().makeVector();
    }

    public void setBox(Box p) {
        super.setBox(p);
        inflate.setBox(p);
        affectedAtomIterator.setBox(p);
    }
     
    public void setInflater(BoxInflate newInflate) {
        inflate = newInflate;
    }
    
    public boolean doTrial() {
        uOld = integrator.getPotentialEnergy();
                
        /*
         * anisotropic scaling: scaleb and scalec
         * they are scale for b-VECTOR and c-VECTOR
         * 
         */
        double scalea = Math.exp((2.*random.nextDouble()-1.)*stepSize);
        double scaleb = scalea;
        double scalec = 1.0/(scalea*scaleb);
        
        //System.out.println("\nbefore: " + box.getBoundary().volume());
        scaleVector.setX(0, scalea);
        scaleVector.setX(1, scaleb);
        scaleVector.setX(2, scalec);
        
        inflate.setVectorScale(scaleVector);
        inflate.actionPerformed();

        potentialCompute.init();
        uNew = potentialCompute.computeAll(false);
        return true;
    }//end of doTrial

    public double getChi(double temperature) {
        return Math.exp(-(uNew - uOld));
    }
    
    public void acceptNotify() {  /* do nothing */}
    
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
