/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.data.DataSourceScalar;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space1d.Space1D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.dimensions.Energy;

/**
 * Meter that calculates the harmonic energy of a configuration given
 * eigenvectors and omegas corresponding to wave vectors.
 * @author Andrew Schultz
 */
public class MeterHarmonicEnergy extends DataSourceScalar {

    public MeterHarmonicEnergy(CoordinateDefinition coordinateDefinition, NormalModes normalModes) {
        super("Harmonic Energy", Energy.DIMENSION);
        this.coordinateDefinition = coordinateDefinition;
        this.normalModes = normalModes;

        int coordinateDim = coordinateDefinition.getCoordinateDim();
        realT = new double[coordinateDim];
        imaginaryT = new double[coordinateDim];

        Box box = coordinateDefinition.getBox();
        normalModes.getWaveVectorFactory().makeWaveVectors(box);
        setWaveVectors(normalModes.getWaveVectorFactory().getWaveVectors(),normalModes.getWaveVectorFactory().getCoefficients());
        setEigenvectors(normalModes.getEigenvectors());
        setOmegaSquared(normalModes.getOmegaSquared());
    }
    
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }

    public double getDataAsScalar() {
        double energySum = 0;
        for (int iVector = 0; iVector < waveVectors.length; iVector++) {
            coordinateDefinition.calcT(waveVectors[iVector], realT, imaginaryT);
            // we want to calculate Q = A T
            // where A is made up of eigenvectors as columns
            int coordinateDim = coordinateDefinition.getCoordinateDim();
            for (int i=0; i<coordinateDim; i++) {
                if (Double.isInfinite(omegaSquared[iVector][i])) {
                    continue;
                }
                double realCoord = 0, imaginaryCoord = 0;
                for (int j=0; j<coordinateDim; j++) {
                    realCoord += eigenvectors[iVector][i][j] * realT[j];
                    imaginaryCoord += eigenvectors[iVector][i][j] * imaginaryT[j];
                }
                // coordinates are now actually the normal mode coordinates divided by sqrt(wvc*2)
                // if wvc=0.5, realCoord and imagCoord are the normal mode coordinates
                // if wvc=1, realCoord and imagCoord are the normal mode coordinates/sqrt(2)
                double normalCoordSq = realCoord*realCoord + imaginaryCoord*imaginaryCoord;
                // if wvc=0.5, then we want 0.5*coord^2
                // if wvc=1.0, we want coord^2, since coord^2 implicitly includes a factor of 0.5
                energySum += waveVectorCoefficients[iVector] * normalCoordSq * omegaSquared[iVector][i];
            }
        }
        return energySum;//don't multiply by 1/2 because we're summing over only half of the wave vectors
    }

    public Box getBox() {
        return coordinateDefinition.getBox();
    }

    protected void setWaveVectors(Vector[] newWaveVectors, double[] coefficients) {
        waveVectors = newWaveVectors;
        waveVectorCoefficients = coefficients;
    }
    
    protected void setEigenvectors(double[][][] eigenvectors) {
        this.eigenvectors = eigenvectors.clone();
    }
    
    protected void setOmegaSquared(double[][] omega2) {
        omegaSquared = new double[omega2.length][omega2[0].length];
        for (int i=0; i<omegaSquared.length; i++) {
            for (int j=0; j<omegaSquared[i].length; j++) {
                // omega is sqrt(kT)/eigenvalue
                omegaSquared[i][j] = omega2[i][j];
            }
        }
    }
    
    private static final long serialVersionUID = 1L;
    protected CoordinateDefinition coordinateDefinition;
    protected double[] realT, imaginaryT;
    protected Vector[] waveVectors;
    protected double[] waveVectorCoefficients;
    protected double[][][] eigenvectors;
    protected double[][] omegaSquared;
    protected NormalModes normalModes;
    
    public static void main(String[] args) {
        
        int numAtoms = 8;
        double L = 10;
        Space sp = Space1D.getInstance();
        Simulation sim = new Simulation(sp);

        SpeciesSpheresMono species = new SpeciesSpheresMono(sim, sp);
        sim.addSpecies(species);

        Box box = new Box(new BoundaryRectangularPeriodic(sim.getSpace(), L), sim.getSpace());
        sim.addBox(box);
        box.setNMolecules(species, numAtoms);

        IAtomList atoms = box.getLeafList();
        
        Primitive primitive = new PrimitiveCubic(sim.getSpace());

        CoordinateDefinition coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, sim.getSpace());
        coordinateDefinition.initializeCoordinates(new int[]{numAtoms});

        for(int i=0; i<numAtoms; i++) {
            atoms.get(i).getPosition().E((i+0.5)*L/numAtoms - 0.5*L);
            System.out.println(atoms.get(i).getPosition().getX(0));
        }
        
        NormalModes normalModes = new NormalModes1DHR(box.getBoundary(), numAtoms);

        MeterHarmonicEnergy meter = new MeterHarmonicEnergy(coordinateDefinition, normalModes);
        atoms.get(1).getPosition().PE(0.5);
        atoms.get(6).getPosition().PE(-1.5);
        System.out.println("Harmonic energy: "+meter.getDataAsScalar());
    }
}
