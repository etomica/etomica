/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.action.BoxInflate;
import etomica.action.BoxInflateDeformable;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Monte Carlo move for changing the box dimensions at constant volume. 
 *
 * @author Andrew Schultz
 */
public class MCMoveBoxSize extends MCMoveBoxStep {
    
    private static final long serialVersionUID = 2L;
    protected double pressure;
    private MeterPotentialEnergy energyMeter;
    protected BoxInflate inflate;
    protected final Space space;
    private IRandom random;
    protected final AtomIteratorLeafAtoms affectedAtomIterator;
    protected final Vector boxScale;

    private transient double uOld, lScale;
    private transient int dim1, dim2;
    private transient double uNew = Double.NaN;

    public MCMoveBoxSize(Simulation sim, PotentialMaster potentialMaster,
                         Space _space) {
        this(potentialMaster, sim.getRandom(), _space);
    }
    
    /**
     * @param potentialMaster an appropriate PotentialMaster instance for calculating energies
     * @param space the governing space for the simulation
     */
    public MCMoveBoxSize(PotentialMaster potentialMaster, IRandom random,
                         Space space) {
        super(potentialMaster);
        this.space = space;
        this.random = random;
        inflate = new BoxInflate(space);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setStepSizeMax(0.5);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        energyMeter.setIncludeLrc(true);
        affectedAtomIterator = new AtomIteratorLeafAtoms();
        boxScale = space.makeVector();
    }
    
    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
        inflate.setBox(p);
        affectedAtomIterator.setBox(p);
    }
    
    public boolean doTrial() {
        if (box.getBoundary() instanceof BoundaryDeformablePeriodic != inflate instanceof BoxInflateDeformable) {
            if (inflate instanceof BoxInflateDeformable) {
                inflate = new BoxInflate(box, space);
            }
            else {
                inflate = new BoxInflateDeformable(box, space);
            }
        }
        uOld = energyMeter.getDataAsScalar();
        lScale = 1 + (2.*random.nextDouble()-1.)*stepSize;
        dim1 = random.nextInt(space.D());
        do {
            dim2 = random.nextInt(space.D());
        } while (dim1 == dim2);
        boxScale.E(1);
        boxScale.setX(dim1, lScale);
        boxScale.setX(dim2, 1.0/lScale);
        inflate.setVectorScale(boxScale);
        inflate.actionPerformed();
        uNew = energyMeter.getDataAsScalar();
        return true;
    }//end of doTrial

    public double getChi(double temperature) {
        return Math.exp(-(uNew - uOld) / temperature);
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
