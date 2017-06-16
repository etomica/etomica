/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.space.*;
import etomica.space.Vector;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.BasisOrthorhombicHexagonal;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.simulation.Simulation;
import etomica.species.SpeciesSpheresMono;
import etomica.math.DoubleRange;
import etomica.data.histogram.HistogramSimple;

/**
 * Simulation class of hard spheres in 1D or 3D that calculates the dq/dx
 * Jacobian.
 */
public class SimModesJ extends Simulation {

    public SimModesJ(Space _space, int numAtoms) {
        super(_space);

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);

        Basis basis;
        int[] nCells;
        primitive = new PrimitiveCubic(space, 1);
        if (space.D() == 1) {
            basis = new BasisMonatomic(space);
            bdry = new BoundaryRectangularPeriodic(space, numAtoms);
            nCells = new int[]{numAtoms};
        }
        else if (space.D() == 2) {
            basis = new BasisOrthorhombicHexagonal();
            int n = (int)Math.round(Math.pow(numAtoms/basis.getScaledCoordinates().length, 1.0/2.0));
            nCells = new int[]{n,n};
            bdry = new BoundaryDeformableLattice(primitive, nCells);
        }
        else {
            basis = new BasisCubicFcc();
            int n = (int)Math.round(Math.pow(numAtoms/basis.getScaledCoordinates().length, 1.0/3.0));
            nCells = new int[]{n,n,n};
            bdry = new BoundaryDeformableLattice(primitive, nCells);
        }
        box.setBoundary(bdry);

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(nCells);
        
        normalModes = new NormalModesVariable(space, space.D()*numAtoms, coordinateDefinition);
        Vector[] waveVectors = normalModes.getWaveVectors();
        double[][] eigenVectors = normalModes.getEigenVectors();
        double[] phaseAngles = normalModes.getPhaseAngles();
        WaveVectorFactory1D waveVectorFactory = new WaveVectorFactory1D();
        waveVectorFactory.makeWaveVectors(box);
        int v = 0;
        boolean useNominalModes = true;
        for (int i=0; i<numAtoms*space.D(); i++) {
            if (space.D() == 1) {
                waveVectors[i] = space.makeVector();
                // density = 1
                waveVectors[i].E(useNominalModes ? waveVectorFactory.getWaveVectors()[v].getX(0) : Math.PI*random.nextDouble());
                eigenVectors[i][0] = 1;
                phaseAngles[i] = useNominalModes ? 0 : random.nextDouble()*2*Math.PI;
                if (useNominalModes && waveVectorFactory.getCoefficients()[v] == 1) {
                    i++;
                    waveVectors[i] = space.makeVector();
                    // density = 1
                    waveVectors[i].E(waveVectors[i-1]);
                    eigenVectors[i][0] = 1;
                    phaseAngles[i] = 0.5*Math.PI;
                }
                v++;
            }
            else {
                waveVectors[i] = space.makeVector();
//                waveVectors[i].E(Math.PI*random.nextDouble());
                eigenVectors[i][i%coordinateDefinition.getCoordinateDim()] = 1;
                phaseAngles[i] = random.nextDouble()*2*Math.PI;
                throw new RuntimeException("please fix my wave vector");
            }
        }

        integrator = new IntegratorMC(null, random, 1.0);
        integrator.setBox(box);
        MCMoveWV moveWV = new MCMoveWV(space, normalModes, coordinateDefinition, random);
        ((MCMoveStepTracker)moveWV.getTracker()).setNoisyAdjustment(true);
        integrator.getMoveManager().addMCMove(moveWV);
        
        MCMovePhaseAngle movePhaseAngle = new MCMovePhaseAngle(space, normalModes, coordinateDefinition, random);
        ((MCMoveStepTracker)movePhaseAngle.getTracker()).setNoisyAdjustment(true);
        integrator.getMoveManager().addMCMove(movePhaseAngle);

        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
    }

    /**
     * @param args
     */
    public static void main(String[] args) {

        // defaults
        int D = 1;
        int nA = 10;

        // parse arguments
        if (args.length > 0) {
            nA = Integer.parseInt(args[0]);
        }

        System.out.println("Calculating "
                + (D == 1 ? "1D" : (D == 3 ? "FCC" : "2D hexagonal"))
                + " hard sphere Jacobian for N="+nA);

        // construct simulation
        final SimModesJ sim = new SimModesJ(Space.getInstance(D), nA);

        final HistogramSimple totalHistogram = new HistogramSimple(nA*10, new DoubleRange(0, 2.0*Math.PI*(nA/2)/nA));
        final HistogramSimple[] histograms = new HistogramSimple[nA-1];
        for (int i=0; i<nA-1; i++) {
            histograms[i] = new HistogramSimple(nA*10, new DoubleRange(0, 2.0*Math.PI*(nA/2)/nA));
        }

        final IntegratorListener dumpHist = new IntegratorListener() {
            public void integratorStepStarted(IntegratorEvent e) {}
            public void integratorStepFinished(IntegratorEvent e) {
                if (++count != 1000) return;
                try {
                    FileWriter hw = new FileWriter("wv_total_histogram.dat");
                    double[] hist = totalHistogram.getHistogram();
                    double[] xval = totalHistogram.xValues();
                    for (int i=0; i<hist.length; i++) {
                        hw.write(xval[i]+" "+hist[i]+"\n");
                    }
                    hw.close();
                    
                    for (int i=0; i<histograms.length; i++) {
                        hw = new FileWriter("wv_histogram"+i+".dat");
                        hist = histograms[i].getHistogram();
                        xval = histograms[i].xValues();
                        for (int j=0; j<hist.length; j++) {
                            hw.write(xval[j]+" "+hist[j]+"\n");
                        }
                        hw.close();
                    }                
                }
                catch (IOException ex) {
                    throw new RuntimeException(ex);
                }
            }
            public void integratorInitialized(IntegratorEvent e) {}
            int count;
        };
        sim.integrator.getEventManager().addListener(dumpHist);
        
        sim.integrator.getEventManager().addListener(new IntegratorListener() {
            
            public void integratorStepStarted(IntegratorEvent e) {}
            
            public void integratorStepFinished(IntegratorEvent e) {
                Vector[] wv = sim.normalModes.getWaveVectors();
                try {
                    FileWriter fw = new FileWriter("wv.out", true);
                    fw.write(sim.integrator.getStepCount()+" ");
                    for (int i=1; i<wv.length; i++) {
                        fw.write(wv[i].getX(0)+" ");
                        wvs[i-1] = wv[i].getX(0);
                        totalHistogram.addValue(wv[i].getX(0));
                    }
                    fw.write("\n");
                    fw.close();
                    Arrays.sort(wvs);
                    for (int i=0; i<wvs.length; i++) {
                        histograms[i].addValue(wvs[i]);
                    }
                }
                catch (IOException ex) {
                    throw new RuntimeException(ex);
                }
            }
            
            public void integratorInitialized(IntegratorEvent e) {
                try {
                    FileWriter fw = new FileWriter("wv.out", false);
                    fw.close();
                }
                catch (IOException ex) {}
                wvs = new double[sim.normalModes.getWaveVectors().length-1];
                integratorStepFinished(e);
            }
            double[] wvs;
        });
        sim.activityIntegrate.setMaxSteps(10000);
        sim.getController().actionPerformed();
        
        dumpHist.integratorStepFinished(null);
    }

    private static final long serialVersionUID = 1L;
    public final ActivityIntegrate activityIntegrate;
    public final IntegratorMC integrator;
    public final Box box;
    public final Boundary bdry;
    public final Primitive primitive;
    public final CoordinateDefinition coordinateDefinition;
    public final NormalModesVariable normalModes;
}
