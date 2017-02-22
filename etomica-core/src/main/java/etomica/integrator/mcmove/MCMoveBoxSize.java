/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.action.BoxInflate;
import etomica.action.BoxInflateDeformable;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.api.IVectorMutable;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.ISpace;

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
    protected final ISpace space;
    private IRandom random;
    protected final AtomIteratorLeafAtoms affectedAtomIterator;
    protected final IVectorMutable boxScale;

    private transient double uOld, lScale;
    private transient int dim1, dim2;
    private transient double uNew = Double.NaN;

    public MCMoveBoxSize(ISimulation sim, IPotentialMaster potentialMaster,
    		            ISpace _space) {
        this(potentialMaster, sim.getRandom(), _space);
    }
    
    /**
     * @param potentialMaster an appropriate PotentialMaster instance for calculating energies
     * @param space the governing space for the simulation
     */
    public MCMoveBoxSize(IPotentialMaster potentialMaster, IRandom random,
    		            ISpace _space) {
        super(potentialMaster);
        this.space = _space;
        this.random = random;
        inflate = new BoxInflate(_space);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setStepSizeMax(0.5);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        energyMeter.setIncludeLrc(true);
        affectedAtomIterator = new AtomIteratorLeafAtoms();
        boxScale = _space.makeVector();
    }
    
    public void setBox(IBox p) {
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
    
    public double getA() {
        return 1;
    }
    
    public double getB() {
        return -(uNew - uOld);
    }
    
    public void acceptNotify() {
//        if (box.getBoundary() instanceof BoundaryDeformablePeriodic) {
//            System.out.print(" => (");
//            double covera = 1;
//            for (int i=0; i<space.D(); i++) {
//                double l = Math.sqrt(((BoundaryDeformablePeriodic)box.getBoundary()).getEdgeVector(i).squared());
//                if (i==2) {
//                    covera *= l;
//                }
//                else {
//                    covera /= Math.sqrt(l);
//                }
//                System.out.print(l+" ");
//            }
//            System.out.println(") "+covera);
//        }
//        else {
//            System.out.println(" => "+box.getBoundary().getBoxSize());
//        }
    }
    
    public void rejectNotify() {
        inflate.undo();
    }

    public double energyChange() {return uNew - uOld;}
    
    public AtomIterator affectedAtoms() {
        return affectedAtomIterator;
    }
}