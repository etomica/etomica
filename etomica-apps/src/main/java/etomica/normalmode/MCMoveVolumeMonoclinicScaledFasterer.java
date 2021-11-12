/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.BoxInflate;
import etomica.action.MoleculeActionTranslateTo;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
import etomica.integrator.IntegratorMCFasterer;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.math.function.Function;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePositionGeometricCenter;
import etomica.potential.compute.PotentialCompute;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Pressure;
import etomica.util.random.IRandom;

/**
 * Monte Carlo volume-change move for simulations of crystalline solids in the
 * NPT ensemble (molecular model). When changing the volume, the molecular COM's 
 * coordinates are scaled away from or toward their lattice in order to 
 * improve the likelihood of acceptance.  
 *
 * @author Tai Boon Tan, Andrew Schultz
 */
public class MCMoveVolumeMonoclinicScaledFasterer extends MCMoveBoxStep {

    protected final PotentialCompute potentialCompute;
    protected final IntegratorMCFasterer integrator;
    protected double pressure;
    protected BoxInflate inflate;
    protected final int D;
    protected IRandom random;
    protected final AtomIteratorLeafAtoms affectedAtomIterator;
    protected double temperature;
    protected final Vector boxSize;
    protected final Vector dest, comOld;
    protected final Vector scaleVector;
    protected final Vector latticeScale;

    protected final MoleculePositionGeometricCenter moleculeCenter;
    protected final MoleculeActionTranslateTo translateTo;

    protected transient double uOld;
    protected transient double uNew = Double.NaN;

    protected Function uLatFunction = uLat0;

    protected Box latticeBox;

    public MCMoveVolumeMonoclinicScaledFasterer(PotentialCompute potentialCompute, IntegratorMCFasterer integrator,
                                                IRandom random, double pressure) {
        this(potentialCompute, integrator, random, pressure, integrator.getBox().getSpace().D());
    }

    public MCMoveVolumeMonoclinicScaledFasterer(PotentialCompute potentialCompute, IntegratorMCFasterer integrator,
                                                IRandom random, double pressure, int D) {
        super();
        this.potentialCompute = potentialCompute;
        this.integrator = integrator;
        this.random = random;
        this.D = D;
        Space space = integrator.getBox().getSpace();
        inflate = new BoxInflate(space);
        setStepSizeMax(0.1);
        setStepSizeMin(0.0);
        setStepSize(0.01);
        setPressure(pressure);
        affectedAtomIterator = new AtomIteratorLeafAtoms();
        
        dest = space.makeVector();
        comOld = space.makeVector();
        boxSize = space.makeVector();
        scaleVector = space.makeVector();
        latticeScale = space.makeVector();
        moleculeCenter = new MoleculePositionGeometricCenter(space);
        translateTo = new MoleculeActionTranslateTo(space);
    }
    
    public void setInflater(BoxInflate newInflate) {
        inflate = newInflate;
    }
    
    public void setLatticeBox(Box newLatticeBox) {
        latticeBox = newLatticeBox;
    }
    
    public void setBox(Box p) {
        super.setBox(p);
        affectedAtomIterator.setBox(p);
    }

    /**
     * Sets the temperature being sampled.  The temperature is used to
     * determine how to scale the atom coordinates.
     * 
     * In actuality, only P/kT is important, but we'll keep the methods
     * separate. 
     */
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }

    /**
     * Sets a function that returns the lattice energy for a given density.
     * If not set, the default lattice energy is taken to be 0 for all densities.
     */
    public void setULatFunction(Function newULatFunction) {
        uLatFunction = newULatFunction;
    }
    
    public Function getULatFunction() {
        return uLatFunction;
    }
    
    public boolean doTrial() {
        uOld = integrator.getPotentialEnergy();

        double scalea = Math.exp((2.*random.nextDouble()-1.)*stepSize);
        
        //System.out.println("\nbefore: " + box.getBoundary().volume());
        scaleVector.setX(0, scalea);
        scaleVector.setX(1, scalea);
        scaleVector.setX(2, 1/(scalea*scalea));
        
        doTransform();

        potentialCompute.init();
        uNew = potentialCompute.computeAll(false);
        return true;
    }
    

    public void doTransform() {

        IMoleculeList moleculeList = box.getMoleculeList();
        IMoleculeList moleculeLatticeList = latticeBox.getMoleculeList();
        int nMolecules = moleculeList.size();

        
        for (int i = 0; i<moleculeList.size(); i++) {
            IMolecule moleculei = moleculeList.get(i);
            IMolecule moleculeLattice = moleculeLatticeList.get(i);
            comOld.E(moleculeCenter.position(moleculei));
            comOld.ME(moleculeCenter.position(moleculeLattice));
            // comOld is now the deviation of the molecule from its lattice site.
            // we now move the molecule to that position
            translateTo.setDestination(comOld);
            translateTo.actionPerformed(moleculei);
        }

        double vOld = box.getBoundary().volume();
        double deltaVAB = vOld * Math.log(scaleVector.getX(0));

        inflate.setVectorScale(scaleVector);
        inflate.setBox(latticeBox);
        inflate.actionPerformed();

        double latticeScaleAB = Math.exp((pressure*deltaVAB)/((nMolecules-1)*temperature*2));
        latticeScale.setX(0, latticeScaleAB);
        latticeScale.setX(1, latticeScaleAB);
        // latticeScaleC = 1/latticeScaleAB^2
        latticeScale.setX(2, Math.exp((-pressure*deltaVAB)/((nMolecules-1)*temperature)));

        for (int i = 0; i<moleculeList.size(); i++) {
            IMolecule moleculei = moleculeList.get(i);
            IMolecule moleculeLattice = moleculeLatticeList.get(i);
            
            dest.E(moleculeCenter.position(moleculei));
            dest.TE(latticeScale);
            dest.PE(moleculeCenter.position(moleculeLattice));
            translateTo.setDestination(dest);
            translateTo.actionPerformed(moleculei);
        }

        if  (box.getBoundary() instanceof BoundaryDeformablePeriodic) {
            for (int i=0; i<dest.getD(); i++) {
                boxSize.E(box.getBoundary().getEdgeVector(i));
                boxSize.TE(scaleVector.getX(i));
                ((BoundaryDeformablePeriodic)box.getBoundary()).setEdgeVector(i, boxSize);
            }
        }
        else {
            boxSize.E(box.getBoundary().getBoxSize());
            boxSize.TE(scaleVector);
            box.getBoundary().setBoxSize(boxSize);
        }
        
    }

    public double getChi(double temp) {
        return Math.exp(-(uNew - uOld) / temp);
    }
    
    public void acceptNotify() {
        /* do nothing */
    }
    
    public void rejectNotify() {
        // invert scaleVector
        latticeScale.E(scaleVector);
        scaleVector.E(1);
        scaleVector.DE(latticeScale);
        doTransform();
        potentialCompute.init();
        potentialCompute.computeAll(false);
    }

    public double energyChange() {return uNew - uOld;}
    
    public AtomIterator affectedAtoms() {
        return affectedAtomIterator;
    }

    public void setPressure(double p) {pressure = p;}
    public final double getPressure() {return pressure;}
    public Dimension getPressureDimension() {return Pressure.DIMENSION;}

    /**
     * Nominal function for lattice energy
     */
    public static Function uLat0 = new Function() {
        public double f(double x) {return 0;}
    };
}
