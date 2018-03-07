/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import etomica.util.random.IRandom;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.*;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
import etomica.potential.Potential2HardSpherical;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Vector;
import etomica.species.SpeciesSpheresMono;
import etomica.units.dimensions.Null;


/**
 * Uses a Widom-like insertion of a mode to calculate a probability.
 * Uses a different box than the main simulation, to assume a mode & rod 
 * are removed
 * 
 * @author cribbin
 *
 */
public class MeterDifferentImageSubtract1D extends DataSourceScalar {

    protected final IRandom random;
    public int nInsert, counter;
    protected double temperature;
    WaveVectorFactory waveVectorFactory;
    private MeterPotentialEnergy meterPE;
    private CoordinateDefinition cDef, simCDef;
    private int cDim, simCDim;
    private Vector[] waveVectors, simWaveVectors;
    private double[] simRealT, simImagT;
    private double[] newU;
    private double[] wvCoeff, simWVCoeff, sqrtWVC;
    private double[][][] eigenVectors, simEigenVectors;
    private double[][] simOmegaSquared;
    private Box box;
    private int numAtoms;
    private Boundary bdry;
    private NormalModes nm;
    
    public MeterDifferentImageSubtract1D(String string,
            int numSimAtoms, double density, Simulation sim,
            Primitive simPrimitive, Basis simBasis, CoordinateDefinition simCD,
            NormalModes simNM, double temp) {
        super(string, Null.DIMENSION);
        this.random = sim.getRandom();
        this.temperature = temp;

        simWaveVectors = simNM.getWaveVectorFactory().getWaveVectors();
        this.simCDef = simCD;
        simCDim = simCD.getCoordinateDim();
        simEigenVectors = simNM.getEigenvectors();
        simWVCoeff = simNM.getWaveVectorFactory().getCoefficients();
        simRealT = new double[simCDim];
        simImagT = new double[simCDim];
        simOmegaSquared = simNM.getOmegaSquared();

        numAtoms = numSimAtoms - 1;
        bdry = new BoundaryRectangularPeriodic(sim.getSpace(), numAtoms / density);
        box = sim.makeBox(bdry);
        box.setNMolecules(sim.getSpecies(0), numAtoms);

        int[] nCells = new int[]{numAtoms};
        cDef = new CoordinateDefinitionLeaf(box, simPrimitive,
                simBasis, sim.getSpace());
        cDef.initializeCoordinates(nCells);
        cDim = cDef.getCoordinateDim();

        nm = new NormalModes1DHR(box.getBoundary(), numAtoms);
        nm.setHarmonicFudge(1.0);
        nm.setTemperature(temperature);
        nm.getOmegaSquared();
        waveVectorFactory = nm.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);
        waveVectors = nm.getWaveVectorFactory().getWaveVectors();
        eigenVectors = nm.getEigenvectors();
        wvCoeff = nm.getWaveVectorFactory().getCoefficients();
        sqrtWVC = new double[wvCoeff.length];
        for (int i = 0; i < wvCoeff.length; i++) {
            sqrtWVC[i] = Math.sqrt(2 * wvCoeff[i]);
        }

        PotentialMasterList potentialMaster = new PotentialMasterList(sim, sim.getSpace());
        Potential2 potential = new P2HardSphere(sim.getSpace(), 1.0, true);
        potential = new P2XOrder(sim.getSpace(), (Potential2HardSpherical) potential);
        potential.setBox(box);
        potentialMaster.addPotential(potential, new AtomType[]{((SpeciesSpheresMono) sim.getSpecies(0)).getLeafType(), ((SpeciesSpheresMono) sim.getSpecies(0)).getLeafType()});
        double neighborRange = 1.01 / density;
        potentialMaster.setRange(neighborRange);
        //find neighbors now.  Don't hook up NeighborListManager since the
        //  neighbors won't change
        potentialMaster.getNeighborManager(box).reset();

        meterPE = new MeterPotentialEnergy(potentialMaster, box);

    }
    
    public double getDataAsScalar() {
        BasisCell[] cells = cDef.getBasisCells();
        BasisCell cell = cells[0];
        //nan this makes it 1D
        newU = new double[cDim];
        for(int i=0; i < newU.length; i++){
            newU[i] = 0.0;
        }
        
        //Calculate normal mode coordinates of simulation system.
        double[] realCoord = new double[simWaveVectors.length];
        double[] imagCoord = new double[simWaveVectors.length];
        for (int wvcount = 0; wvcount < simWaveVectors.length; wvcount++){
            simCDef.calcT(simWaveVectors[wvcount], simRealT, simImagT);
            realCoord[wvcount] = 0.0;
            imagCoord[wvcount] = 0.0;
            for (int i = 0; i < simCDim; i++){
                for (int j = 0; j < simCDim; j++){
                    realCoord[wvcount] += simEigenVectors[wvcount][i][j] * simRealT[j];
                    imagCoord[wvcount] += simEigenVectors[wvcount][i][j] * simImagT[j];
                }
            }
            if(simWVCoeff[wvcount] == 1.0){
                realCoord[wvcount] *= Math.sqrt(2);
                imagCoord[wvcount] *= Math.sqrt(2);
            } else {
                imagCoord[wvcount] = 0.0;
            }
        }
        
        //nan this will not work for more than 1D
        //Calculation of harmonic energy
        //The if statement determines whether the real coordinate or the 
        // imaginary coordinate is used, and zeroes out the last coordinate
        // calculated.
        double harmonic = 0.0;
        int index = simWVCoeff.length-1;
        if(simWVCoeff[index] == 1.0){
            harmonic = 0.5 * simOmegaSquared[index][0] * 
                    imagCoord[index] * imagCoord[index];
            imagCoord[index] = 0.0;
        } else {
            harmonic = 0.5 * simOmegaSquared[index][0] * 
                    realCoord[index] * realCoord[index];
            realCoord[index] = 0.0;
        }
        
        //Calculate the positions for the meter's system
        for (int iCell = 0; iCell < cells.length; iCell++){
            cell = cells[iCell];
            for (int j = 0; j < cDim; j++) {
                newU[j] = 0.0;
            }
            for (int wvcount = 0; wvcount < waveVectors.length; wvcount++){
                //Calculate the change in positions.
                double kR = waveVectors[wvcount].dot(cell.cellPosition);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                for (int i = 0; i < cDim; i++){
                    for (int j = 0; j < cDim; j++){
                       newU[j] += sqrtWVC[wvcount] * eigenVectors[wvcount][i][j] 
                            * (realCoord[wvcount] * coskR - imagCoord[wvcount] * sinkR);
                    }
                }
            }
            
            double normalization = 1/Math.sqrt(cells.length);
            for (int i=0; i<cDim; i++) {
                newU[i] *= normalization;
            }
            cDef.setToU(cells[iCell].molecules, newU);
        }
        
//        System.out.println("1d harmonic " + harmonic+ "  "+ meterPE.getDataAsScalar());
//        return Math.exp(-1*(meterPE.getDataAsScalar()+harmonic));
        return meterPE.getDataAsScalar() + harmonic;
    }

    
}
