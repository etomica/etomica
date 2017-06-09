/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.BoxInflate;
import etomica.action.MoleculeActionTranslateTo;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.math.function.Function;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePositionGeometricCenter;
import etomica.potential.PotentialMaster;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Dimension;
import etomica.units.Pressure;
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
    
    protected double pressure;
    protected MeterPotentialEnergy energyMeter;
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

    /**
     * @param potentialMaster an appropriate PotentialMaster instance for calculating energies
     * @param space the governing space for the simulation
     */
    public MCMoveVolumeSolidNPTMolecular(PotentialMaster potentialMaster, IRandom random,
                                         Space space, double pressure) {
        this(potentialMaster, random, space, pressure, space.D());
    }
    
    public MCMoveVolumeSolidNPTMolecular(PotentialMaster potentialMaster, IRandom random,
                                         Space space, double pressure, int D) {
        super(potentialMaster);
        this.random = random;
        this.D = D;
        inflate = new BoxInflate(space);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setStepSizeMax(0.1);
        setStepSizeMin(0.0);
        setStepSize(0.001);
        setPressure(pressure);
        energyMeter.setIncludeLrc(true);
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
        energyMeter.setBox(p);
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
        uOld = energyMeter.getDataAsScalar();
        if (uOld == Double.POSITIVE_INFINITY) {
            throw new RuntimeException("oops initial inf");
        }
        vOld = box.getBoundary().volume();
        vScale = (2.*random.nextDouble()-1.)*stepSize;
        vNew = vOld * Math.exp(vScale); //Step in ln(V)
        doTransform(vScale);
        
        uNew = energyMeter.getDataAsScalar();
        return true;
    }
    

    public void doTransform(double vScaleLocal) {

        IMoleculeList moleculeList = box.getMoleculeList();
        IMoleculeList moleculeLatticeList = latticeBox.getMoleculeList();
        int nMolecules = moleculeList.getMoleculeCount();

        for (int i=0; i<moleculeList.getMoleculeCount(); i++) {
            IMolecule moleculei = moleculeList.getMolecule(i);
            IMolecule moleculeLattice = moleculeLatticeList.getMolecule(i);
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

        for (int i=0; i<moleculeList.getMoleculeCount(); i++) {
            IMolecule moleculei = moleculeList.getMolecule(i);
            IMolecule moleculeLattice = moleculeLatticeList.getMolecule(i);
            
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
    
    public double getA() {
        return 1;
    }
    
    public double getB() {
        int nMolecules = box.getMoleculeList().getMoleculeCount();
        double uLatOld = nMolecules*uLatFunction.f(nMolecules/vOld);
        double uLatNew = nMolecules*uLatFunction.f(nMolecules/vNew);
        return -((uNew-uLatNew) - (uOld-uLatOld));
    }
    
    public void acceptNotify() {
        /* do nothing */
    }
    
    public void rejectNotify() {
        doTransform(-vScale);
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
