/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.BoxInflate;
import etomica.action.MoleculeActionTranslateTo;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
import etomica.integrator.IntegratorMC;
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
public class MCMoveVolumeSolidNPTMolecular extends MCMoveBoxStep {

    protected final IntegratorMC integrator;
    protected final PotentialCompute potentialCompute;
    protected double pressure;
    protected BoxInflate inflate;
    protected final int D;
    protected IRandom random;
    protected final AtomIteratorLeafAtoms affectedAtomIterator;
    protected double temperature;
    protected final Vector boxSize;
    protected final Vector dest, comOld;

    protected final MoleculePositionGeometricCenter moleculeCenter;
    protected final MoleculeActionTranslateTo translateTo;

    protected transient double uOld, vOld, vNew, vScale;
    protected transient double uNew = Double.NaN, latticeScale;

    protected Function uLatFunction = uLat0;

    protected Box latticeBox;

    public MCMoveVolumeSolidNPTMolecular(PotentialCompute potentialCompute, IntegratorMC integrator,
                                         IRandom random, double pressure) {
        this(potentialCompute, integrator, random, pressure, integrator.getBox().getSpace().D());
    }

    public MCMoveVolumeSolidNPTMolecular(PotentialCompute potentialCompute, IntegratorMC integrator,
                                         IRandom random, double pressure, int D) {
        super();
        this.integrator = integrator;
        this.potentialCompute = potentialCompute;
        this.random = random;
        this.D = D;
        Space space = integrator.getBox().getSpace();
        inflate = new BoxInflate(space);
        setStepSizeMax(0.1);
        setStepSizeMin(0.0);
        setStepSize(0.001);
        setPressure(pressure);
        affectedAtomIterator = new AtomIteratorLeafAtoms();
        
        dest = space.makeVector();
        comOld = space.makeVector();
        boxSize = space.makeVector();
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
        if (uOld == Double.POSITIVE_INFINITY) {
            throw new RuntimeException("oops initial inf");
        }
        vOld = box.getBoundary().volume();
        vScale = (2.*random.nextDouble()-1.)*stepSize;
        vNew = vOld * Math.exp(vScale); //Step in ln(V)
        doTransform(vScale);

        potentialCompute.init();
        uNew = potentialCompute.computeAll(false);
        return true;
    }
    

    public void doTransform(double vScaleLocal) {

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

        double vOldLocal = box.getBoundary().volume();
        double vNewLocal = vOldLocal * Math.exp(vScaleLocal); //Step in ln(V)
        double rScale = Math.exp(vScaleLocal/boxSize.getD());
        inflate.setScale(rScale);

        inflate.setBox(latticeBox);
        inflate.actionPerformed();

        double uLatOld = nMolecules*uLatFunction.f(nMolecules/vOldLocal);
        double uLatNew = nMolecules*uLatFunction.f(nMolecules/vNewLocal);

        int transDim = boxSize.getD();
        latticeScale = Math.exp((pressure*(vNewLocal-vOldLocal)+(uLatNew-uLatOld) - vScaleLocal)/((nMolecules*D-transDim)*temperature));

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
                boxSize.TE(rScale);
                ((BoundaryDeformablePeriodic)box.getBoundary()).setEdgeVector(i, boxSize);
            }
        }
        else {
            boxSize.E(box.getBoundary().getBoxSize());
            boxSize.TE(rScale);
            box.getBoundary().setBoxSize(boxSize);
        }
        
    }

    public double getChi(double temperature) {
        int nMolecules = box.getMoleculeList().size();
        double uLatOld = nMolecules*uLatFunction.f(nMolecules/vOld);
        double uLatNew = nMolecules*uLatFunction.f(nMolecules/vNew);
        return Math.exp(-((uNew - uLatNew) - (uOld - uLatOld)) / temperature);
    }
    
    public void acceptNotify() {
        /* do nothing */
    }
    
    public void rejectNotify() {
        doTransform(-vScale);
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
