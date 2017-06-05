/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.util.random.IRandom;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.normalmode.CoordinateDefinition.BasisCell;

/**
 * MCMove that performs random displacements in the harmonic coordinates for
 * the k=0 wave vector.  The acceptance criteria is the change in target energy
 * for the move.  The class moves molecules in all cells together.  If they
 * are different before the move trial, the trial will homogenize the system.
 *
 * @author Andrew Schultz
 */
public class MCMoveHarmonicStep extends MCMoveBoxStep {

    public MCMoveHarmonicStep(PotentialMaster potentialMaster, IRandom random) {
        super(potentialMaster);
        
        this.random = random;
        iterator = new AtomIteratorLeafAtoms();
        energyMeter = new MeterPotentialEnergy(potentialMaster);
    }
    
    public void setCoordinateDefinition(CoordinateDefinition newCoordinateDefinition) {
        coordinateDefinition = newCoordinateDefinition;
        u = new double[coordinateDefinition.getCoordinateDim()];
        uOld = new double[coordinateDefinition.getCoordinateDim()];
    }
    
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }

    /**
     * Sets the indices of the harmonic modes of the k=0 wave vector which
     * should be sampled.
     */
    public void setModes(int[] newModes) {
        modes = newModes;
    }
    
    /**
     * Informs the move of the eigenvectors for the k=0 wave vector.  The
     * actual eigenvectors used will be those specified via setModes
     */
    public void setEigenVectors(double[][] newEigenVectors) {
        eigenVectors = newEigenVectors;
    }
    
    public void setBox(Box newBox) {
        super.setBox(newBox);
        iterator.setBox(newBox);
        energyMeter.setBox(newBox);
        latticeEnergy = energyMeter.getDataAsScalar();
    }

    public AtomIterator affectedAtoms() {
        return iterator;
    }

    public boolean doTrial() {
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        BasisCell[] cells = coordinateDefinition.getBasisCells();

        BasisCell cell = cells[0];
        energyOld = energyMeter.getDataAsScalar();
        double[] calcedU = coordinateDefinition.calcU(cell.molecules);
        
        double sqrtCells = Math.sqrt(cells.length);
        for (int i=0; i<coordinateDim; i++) {
            uOld[i] = calcedU[i];
            u[i] = calcedU[i];
        }
        for (int i=0; i<modes.length; i++) {
            double delta = (2*random.nextDouble()-1) * stepSize;
            for (int j=0; j<coordinateDim; j++) {
                u[j] += delta * eigenVectors[modes[i]][j] / sqrtCells;
            }
        }

        // Set all the atoms to the new values of u
        for (int iCell = 0; iCell<cells.length; iCell++) {
            coordinateDefinition.setToU(cells[iCell].molecules, u);
        }
        
        energyNew = energyMeter.getDataAsScalar();
        return true;
    }
    
    public double getA() {
        return 1;
    }

    public double getB() {
        return -(energyNew - energyOld);
    }
    
    public void acceptNotify() {
    }

    public double energyChange() {
        return energyNew - energyOld;
    }

    public void rejectNotify() {
        // Set all the atoms back to the old values of u
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        for (int iCell = 0; iCell<cells.length; iCell++) {
            coordinateDefinition.setToU(cells[iCell].molecules, uOld);
        }
    }

    private static final long serialVersionUID = 1L;
    protected CoordinateDefinition coordinateDefinition;
    protected final AtomIteratorLeafAtoms iterator;
    protected double[] uOld, u;
    private double[][] eigenVectors;
    protected final IRandom random;
    protected double energyOld, energyNew, latticeEnergy;
    protected final MeterPotentialEnergy energyMeter;
    protected int[] modes;
}
