/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.ensembles;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.DataSplitter;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.*;
import etomica.potential.IPotentialAtomic;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesGeneral;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class LJMC1D extends Simulation {
    // The final keyword tells use that each of these variables are created once.
    public final SpeciesGeneral species;
    public final Box box;
    public final IntegratorMC integrator;
    public final MCMoveAtomCoupled mcMoveAtom;
    public final MCMoveHarmonicStep mcMoveHarmonicStep;
    public PotentialMasterList potentialMaster;
    public final CoordinateDefinitionLeaf coordinates;


    public LJMC1D(Space _space) {
        this(_space, 3.0, 150, 4.1, 2);
    }

    public LJMC1D(Space _space, double temperature, int N, double truncationRadius, double rho) {
        this(_space, temperature, N, truncationRadius, rho, new int[0]);
    }

    public LJMC1D(Space _space, double temperature, int N, double truncationRadius, double rho, int[] modes) {
        super(_space);

        // Species.
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        potentialMaster = new PotentialMasterList(this, 6.5, space);        // This holds the various potentials to be calculated during the simulation.
        potentialMaster.setCellRange(2);

        // Controller and integrator.
        box = this.makeBox();
        integrator = new IntegratorMC(potentialMaster, random, temperature, box);

        P2LennardJones potential = new P2LennardJones(space);
        P2SoftSphericalTruncated p2Truncated = new P2SoftSphericalTruncated(space, potential, truncationRadius);
        potentialMaster.addPotential(p2Truncated, new AtomType[]{species.getLeafType(), species.getLeafType()});

        //construct box
        Vector dim = space.makeVector(1);
        double L = N / rho;
        dim.E(L);
        box.getBoundary().setBoxSize(dim);
        box.setNMolecules(species, N);
        BasisBigCell basis = new BasisBigCell(space, new BasisMonatomic(space), new int[]{N});
        coordinates = new CoordinateDefinitionLeaf(box, new PrimitiveCubic(space, L), basis, space);
        coordinates.initializeCoordinates(new int[]{1});
        potentialMaster.reset();
        p2Truncated.setTruncationRadius(100);       // returns the strength, irrespective of how far the atoms are.

        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster, box);
        if (modes.length == 0) {
            mcMoveAtom = new MCMoveAtomCoupled(potentialMaster, meterPE, random, space);
            mcMoveAtom.setPotential(new IPotentialAtomic[]{p2Truncated});
            mcMoveAtom.setDoExcludeNonNeighbors(true);
            integrator.getMoveManager().addMCMove(mcMoveAtom);
            mcMoveHarmonicStep = null;
        } else {
            mcMoveAtom = null;
            mcMoveHarmonicStep = new MCMoveHarmonicStep(potentialMaster, random);
            double[][] eigenvectors = new double[modes.length][N];
            int[] moveModes = new int[modes.length];
            int maxmode = N / 2;
            for (int i = 0; i < modes.length; i++) {
                moveModes[i] = i;
                boolean doCos = modes[i] <= N / 2;
                int k = doCos ? modes[i] : (modes[i] - maxmode);
                for (int j = 0; j < N; j++) {
                    double arg = -2 * Math.PI / N * k * j;
                    eigenvectors[i][j] = 2 * (doCos ? Math.cos(arg) : Math.sin(arg)) / Math.sqrt(N);
                }
            }
            mcMoveHarmonicStep.setCoordinateDefinition(coordinates);
            mcMoveHarmonicStep.setModes(moveModes);
            mcMoveHarmonicStep.setEigenVectors(eigenvectors);
            integrator.getMoveManager().addMCMove(mcMoveHarmonicStep);
        }
    }

    public static void main(String[] args) throws IOException {
        Space space = Space.getInstance(1);

//        String path = "/Users/sykherebrown/Masters Research/1D Harmonic Crystal Output/";

        double[] temperatureList = new double[] { 0.01 };
        int [] systemSizeList = new int[] { 10 };
        double [] truncationRadiusList = new double[] { 2.1 };
        double [] rhoList = new double[] { 1.25 };
        int [] modesList = new int[] { 1 };

//        double[] temperatureList = new double[] { 1.0 };
//        int [] systemSizeList = new int[] { 75 };
//        double [] truncationRadiusList = new double[] { 2.1 };
//        double [] rhoList = new double[] { 1.25 };
//        int [] modesList = new int[] { 0, 10, 20, 30, 40, 50, 60, 70 };

        long steps = 10000;

        // Creating File and FileWriter objects.
//        String filename = "1DHarmonicCrystalPressureNormalModes.txt";
//        // String filename = "1DHarmonicCrystalPressureData.txt";
//        File file = new File(path + filename);
//        FileWriter fileWriter = new FileWriter(file);
//
//        fileWriter.write("Temperature" + "," + "N" + "," + "TruncationRadius" + "," + "Rho" +
//                "," + "NormalMode" + "," + "ConventionalPressure" + "," + "ConventionalPressureError" +
//                "," + "pressureHMANormalMode" + "," + "pressureHMANormalModeError" + ","
//                + "pressureHMARealSpace" + "," + "pressureHMARealSpaceError" + "\n");


        // The code below is deeply nested, is there a better way to configure the code?
        for(int mode : modesList) {
            for(double temperature : temperatureList) {
                for(int N : systemSizeList) {
                    for(double truncationRadius : truncationRadiusList) {
                        for(double rho : rhoList) {

                            double latticeConstant = 1 / rho;

                            if (truncationRadius / latticeConstant <= 1.001)
                                continue;

                            LJMC1D sim = new LJMC1D(space, temperature, N, truncationRadius, rho, new int[] {1, 2, 3, 4, 5, 6, 7, 8, 9});

                            MeterPressureHMA pMeterHMA = new MeterPressureHMA(space, sim.potentialMaster, sim.coordinates, true, mode);
                            pMeterHMA.setTemperature(sim.integrator.getTemperature());
                            pMeterHMA.setTruncationRadius(truncationRadius);
                            pMeterHMA.setPRes();        // shouldn't this be done automatically?
                            pMeterHMA.calculateGij();

                            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps/10));

                            // NOTE: The pressures that are being written to the file need to be obtained from the AccumulatorAverageFixed().
                            final AccumulatorAverageFixed pAccumulatorConventional = new AccumulatorAverageFixed();
                            final AccumulatorAverageFixed pAccumulatorHMARealSpace = new AccumulatorAverageFixed();
                            final AccumulatorAverageFixed pAccumulatorHMANormalMode = new AccumulatorAverageFixed();

                            // TODO: Pump directly into the Accumulator from the Meter.
                            DataSplitter splitter = new DataSplitter();
                            final DataPumpListener pPump = new DataPumpListener(pMeterHMA, splitter, N);
                            splitter.setDataSink(1, pAccumulatorConventional);
                            splitter.setDataSink(4, pAccumulatorHMANormalMode);
                            splitter.setDataSink(5, pAccumulatorHMARealSpace);

                            sim.integrator.getEventManager().addListener(pPump);

                            // Production code.
                            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));

                            // TODO: Is the standard deviation of the pressure more useful, compared to the error?
                            int pressureIndex = pAccumulatorConventional.AVERAGE.index;      // index 1 represents the pressure's average.
                            int errorIndex = pAccumulatorConventional.STANDARD_DEVIATION.index;          // index 2 represents the error in the pressure's average.

                            double pressureConventional = pAccumulatorConventional.getData().getValue(pressureIndex);
                            double pressureConventionalError = pAccumulatorConventional.getData().getValue(errorIndex);

                            double pressureHMANormalMode = pAccumulatorHMANormalMode.getData().getValue(pressureIndex);
                            double pressureHMANormalModeError = pAccumulatorHMANormalMode.getData().getValue(errorIndex);

                            double pressureHMARealSpace = pAccumulatorHMARealSpace.getData().getValue(pressureIndex);
                            double pressureHMARealSpaceError = pAccumulatorHMARealSpace.getData().getValue(errorIndex);

                            // Does this run? Is this only being called once? Is it being called multiple times?
                            double[] q = sim.mcMoveHarmonicStep.q;

//                            fileWriter.write(temperature + "," + N + "," + truncationRadius + "," + rho +
//                                    "," + q[0] + "," + pressureConventional + "," + pressureConventionalError +
//                                    "," + pressureHMANormalMode + "," + pressureHMANormalModeError + ","
//                                    + pressureHMARealSpace + "," + pressureHMARealSpaceError + "\n");
//
//                            fileWriter.flush();

                            pMeterHMA.closeFile();
                            System.out.println(temperature + "," + N + "," + truncationRadius + "," + rho + "," + pressureConventional + "," + pressureConventionalError +
                                    "," + pressureHMANormalMode + "," + pressureHMANormalModeError + ","
                                    + pressureHMARealSpace + "," + pressureHMARealSpaceError + "\n");
                        }
                    }
                }
            }

        }


        // fileWriter.close();
    }
}
