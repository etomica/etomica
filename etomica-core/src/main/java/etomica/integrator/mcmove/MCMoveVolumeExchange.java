/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.action.BoxInflate;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorNull;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorMC;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.util.random.IRandom;

/**
 * Elementary Monte Carlo trial that exchanges volume between two boxs.  Trial
 * consists of a volume increase in one box (selected at random) and an equal
 * volume decrease in the other.  Used in Gibbs ensemble simulations.
 *
 * @author David Kofke
 */
public class MCMoveVolumeExchange extends MCMoveStep {
    
    private static final long serialVersionUID = 1L;
    private final MeterPotentialEnergy energyMeter;
    protected final Box firstBox;
    protected final Box secondBox;
    private final IntegratorBox integrator1;
    private final IntegratorBox integrator2;
    private final BoxInflate inflate1;
    private final BoxInflate inflate2;
    private transient double uOld1, uOld2;
    private transient double uNew1 = Double.NaN;
    private transient double uNew2 = Double.NaN;
    private final double ROOT;
    private final AtomIteratorLeafAtoms box1AtomIterator;
    private final AtomIteratorLeafAtoms box2AtomIterator;
    private final IRandom random;
    
    private transient double hOld, v1Scale, v2Scale;

    public MCMoveVolumeExchange(PotentialMaster potentialMaster, IRandom random,
                                Space _space,
                                IntegratorBox integrator1,
                                IntegratorBox integrator2) {
        super(potentialMaster, new MCMoveStepTracker());
        this.random = random;
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        ROOT = 1.0/_space.D();
        setStepSizeMax(Double.MAX_VALUE);
        setStepSizeMin(Double.MIN_VALUE);
        setStepSize(0.1);
        box1AtomIterator = new AtomIteratorLeafAtoms();
        box2AtomIterator = new AtomIteratorLeafAtoms();
        energyMeter.setIncludeLrc(true);
        inflate1 = new BoxInflate(_space);
        inflate2 = new BoxInflate(_space);
        this.integrator1 = integrator1;
        this.integrator2 = integrator2;
        firstBox = integrator1.getBox();
        secondBox = integrator2.getBox();
        inflate1.setBox(firstBox);
        inflate2.setBox(secondBox);
        box1AtomIterator.setBox(firstBox);
        box2AtomIterator.setBox(secondBox);
    }
    
    public boolean doTrial() {
        uOld1 = integrator1.getPotentialEnergy();
        uOld2 = integrator2.getPotentialEnergy();
        hOld = uOld1 + uOld2;
        double v1Old = firstBox.getBoundary().volume();
        double v2Old = secondBox.getBoundary().volume();
        double step = stepSize * (random.nextDouble() - 0.5); 
        double vRatio = v1Old/v2Old * Math.exp(step);
        double v2New = (v1Old + v2Old)/(1 + vRatio);
        double v1New = (v1Old + v2Old - v2New);
        v1Scale = v1New/v1Old;
        v2Scale = v2New/v2Old;
        inflate1.setScale(Math.pow(v1Scale,ROOT));
        inflate2.setScale(Math.pow(v2Scale,ROOT));
        // for cell-listing, this will trigger NeighborCellManager notification
        // (as a box listener)
        inflate1.actionPerformed();
        inflate2.actionPerformed();
        return true;
    }//end of doTrial

    public double getChi(double temperature) {
        energyMeter.setBox(firstBox);
        uNew1 = energyMeter.getDataAsScalar();
        energyMeter.setBox(secondBox);
        uNew2 = energyMeter.getDataAsScalar();
        double hNew = uNew1 + uNew2;
        double B = -(hNew - hOld);
        // assume both integrators have the same temperature
        return Math.exp(B / temperature) * Math.pow(v1Scale, (firstBox.getMoleculeList().getMoleculeCount() + 1))
                * Math.pow(v2Scale,(secondBox.getMoleculeList().getMoleculeCount()+1));
    }
    
    public void acceptNotify() {
        if (integrator1 instanceof IntegratorMC) {
            ((IntegratorMC)integrator1).notifyEnergyChange(uNew1-uOld1);
        }
        else {
            //XXX grossly inefficient
            integrator1.reset();
        }
        if (integrator2 instanceof IntegratorMC) {
            ((IntegratorMC)integrator2).notifyEnergyChange(uNew2-uOld2);
        }
        else {
            //XXX grossly inefficient
            integrator2.reset();
        }
    }
    
    public void rejectNotify() {
        inflate1.undo();
        inflate2.undo();
    }

    public double energyChange(Box box) {
        if(this.firstBox == box) return uNew1 - uOld1;
        else if(this.secondBox == box) return uNew2 - uOld2;
        else return 0.0;
    }
    
    public final AtomIterator affectedAtoms(Box box) {
        if(this.firstBox == box) {
            return box1AtomIterator;
        } else if(this.secondBox == box) {
            return box2AtomIterator;
        } else {
            return AtomIteratorNull.INSTANCE;
        }
    }

}
