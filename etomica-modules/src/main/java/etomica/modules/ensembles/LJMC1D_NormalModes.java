/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.ensembles;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveInsertDelete;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.integrator.mcmove.MCMoveVolume;
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
import java.io.FileWriter;
import java.io.File;

import java.io.IOException;

public class LJMC1D_NormalModes extends Simulation {
    private static final long serialVersionUID = 1L;
    public final SpeciesGeneral species;
    public final Box box;
    public final ActivityIntegrate activityIntegrate;
    public final IntegratorMC integrator;
    public final MCMoveAtomCoupled mcMoveAtom;
    public final MCMoveVolume mcMoveVolume;
    public final MCMoveInsertDelete mcMoveID;
    public PotentialMasterList potentialMaster;
    public final CoordinateDefinitionLeaf coordinates;
    public final MCMoveHarmonicStep mcMoveMode;


    public LJMC1D_NormalModes(Space _space) {
        super(_space);

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);

        potentialMaster = new PotentialMasterList(this, 6.5, space);
        potentialMaster.setCellRange(2);
        int N = 150;

        // Controller and integrator
        box = this.makeBox();
        integrator = new IntegratorMC(potentialMaster, random, 3.0, box);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addActivity(activityIntegrate);

        P2LennardJones potential = new P2LennardJones(space);
        P2SoftSphericalTruncated p2Truncated = new P2SoftSphericalTruncated(space, potential, 4.1);
        potentialMaster.addPotential(p2Truncated, new AtomType[]{species.getLeafType(), species.getLeafType()});
        p2Truncated.setTruncationRadius(100);

        //construct box
        Vector dim = space.makeVector(1);
        double L = 75;
        dim.E(L);
        box.getBoundary().setBoxSize(dim);
        box.setNMolecules(species, N);
        coordinates = new CoordinateDefinitionLeaf(box, new PrimitiveCubic(space, L / N), space);
        coordinates.initializeCoordinates(new int[]{N});
        potentialMaster.reset();

        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster, box);

        // Below is the block of code concerned with normal-mode sampling.
        mcMoveMode = new MCMoveHarmonicStep(potentialMaster, random);
        mcMoveMode.setCoordinateDefinition(coordinates);
        mcMoveMode.setBox(box);
        mcMoveMode.setModes(new int [N]);
        integrator.getMoveManager().addMCMove(mcMoveMode);
        ((MCMoveStepTracker) mcMoveMode.getTracker()).setMaxAdjustInterval(50000);
        ((MCMoveStepTracker) mcMoveMode.getTracker()).setMinAdjustStep(1.05);

        mcMoveAtom = new MCMoveAtomCoupled(potentialMaster, meterPE, random, space);
        mcMoveAtom.setPotential(new IPotentialAtomic[]{p2Truncated});
        mcMoveAtom.setDoExcludeNonNeighbors(true);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        ((MCMoveStepTracker) mcMoveAtom.getTracker()).setMaxAdjustInterval(50000);
        ((MCMoveStepTracker) mcMoveAtom.getTracker()).setMinAdjustStep(1.05);

        mcMoveVolume = new MCMoveVolume(potentialMaster, random, space, 1) {
            public boolean doTrial() {
                double vOld = box.getBoundary().volume();
                uOld = energyMeter.getDataAsScalar();
                hOld = uOld + pressure * vOld;
                biasOld = vBias.f(vOld);
                vScale = (2. * random.nextDouble() - 1.) * stepSize;
                vNew = vOld * Math.exp(vScale); //Step in ln(V)
                if (vNew < Math.pow(6, space.D())) return false;
                double rScale = Math.exp(vScale / D);
                inflate.setScale(rScale);
                inflate.actionPerformed();
                uNew = energyMeter.getDataAsScalar();
                hNew = uNew + pressure * vNew;
                return true;
            }
        };

        mcMoveID = new MCMoveInsertDelete(potentialMaster, random, space) {
            public void acceptNotify() {
                if (moleculeList.size() > 999 && insert) {
                    rejectNotify();
                    uNew = 0;
                    return;
                }
                super.acceptNotify();
            }
        };
        mcMoveID.setSpecies(species);
    }

    public LJMC1D_NormalModes(Space _space, double temperature, int N, double truncationRadius, double rho) {
        super(_space);

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);

        potentialMaster = new PotentialMasterList(this, 6.5, space);
        potentialMaster.setCellRange(2);

        // Controller and integrator.
        box = this.makeBox();
        integrator = new IntegratorMC(potentialMaster, random, temperature, box);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addActivity(activityIntegrate);

        P2LennardJones potential = new P2LennardJones(space);
        P2SoftSphericalTruncated p2Truncated = new P2SoftSphericalTruncated(space, potential, truncationRadius);
        potentialMaster.addPotential(p2Truncated, new AtomType[]{species.getLeafType(), species.getLeafType()});

        //construct box
        Vector dim = space.makeVector(1);
        double L = N / rho;
        dim.E(L);
        box.getBoundary().setBoxSize(dim);
        box.setNMolecules(species, N);
        coordinates = new CoordinateDefinitionLeaf(box, new PrimitiveCubic(space, 1 / rho), space);
        coordinates.initializeCoordinates(new int[]{N});
        potentialMaster.reset();
        p2Truncated.setTruncationRadius(100);       // returns the strength, irrespective of how far the atoms are.

        // Below is the block of code concerned with normal-mode sampling.
        mcMoveMode = new MCMoveHarmonicStep(potentialMaster, random);
        mcMoveMode.setCoordinateDefinition(coordinates);
        mcMoveMode.setBox(box);
        mcMoveMode.setModes(new int [N]);
        integrator.getMoveManager().addMCMove(mcMoveMode);
        ((MCMoveStepTracker) mcMoveMode.getTracker()).setMaxAdjustInterval(50000);
        ((MCMoveStepTracker) mcMoveMode.getTracker()).setMinAdjustStep(1.05);

        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster, box);
        mcMoveAtom = new MCMoveAtomCoupled(potentialMaster, meterPE, random, space);
        mcMoveAtom.setPotential(new IPotentialAtomic[]{p2Truncated});
        mcMoveAtom.setDoExcludeNonNeighbors(true);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        ((MCMoveStepTracker) mcMoveAtom.getTracker()).setMaxAdjustInterval(50000);
        ((MCMoveStepTracker) mcMoveAtom.getTracker()).setMinAdjustStep(1.05);

        mcMoveVolume = new MCMoveVolume(potentialMaster, random, space, 1) {
            public boolean doTrial() {
                double vOld = box.getBoundary().volume();
                uOld = energyMeter.getDataAsScalar();
                hOld = uOld + pressure * vOld;
                biasOld = vBias.f(vOld);
                vScale = (2. * random.nextDouble() - 1.) * stepSize;
                vNew = vOld * Math.exp(vScale); //Step in ln(V)
                if (vNew < Math.pow(6, space.D())) return false;
                double rScale = Math.exp(vScale / D);
                inflate.setScale(rScale);
                inflate.actionPerformed();
                uNew = energyMeter.getDataAsScalar();
                hNew = uNew + pressure * vNew;
                return true;
            }
        };

        mcMoveID = new MCMoveInsertDelete(potentialMaster, random, space) {
            public void acceptNotify() {
                if (moleculeList.size() > 999 && insert) {
                    rejectNotify();
                    uNew = 0;
                    return;
                }
                super.acceptNotify();
            }
        };
        mcMoveID.setSpecies(species);
    }

    public static void main(String[] args) throws IOException {
        Space space = Space.getInstance(1);

        String path = "/Users/sykherebrown/Masters Research/1D Harmonic Crystal Output/";

        // Below are list a lists of parameters that we're varying.
        double[] temperatureList = new double[] { 0.1, 0.5, 1.0, 2.0, 3.0 };
        int [] systemSizeList = new int[] { 75, 150 };
        double [] truncationRadiusList = new double[] { 1.1, 2.1, 3.1, 4.1 };
        double [] rhoList = new double[] { 0.2, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0 };

//        double[] temperatureList = new double[] { 1.0 };
//        int [] systemSizeList = new int[] { 150 };
//        double [] truncationRadiusList = new double[] { 2.1 };
//        double [] rhoList = new double[] { 1.25 };

        long steps = 1000000;

        // Creating File and FileWriter objects.
        String filename = "1DHarmonicCrystalNormalModesSamplingAnharmonicPressure.txt";
//        String filename = "1DHarmonicCrystalNormalModesSamplingData.txt";
//        String filename = "1DHarmonicCrystalPressureDataAnharmonicDividedByTemp.txt";
//        String filename = "1DHarmonicCrystalPressureData.txt";

        File file = new File(path + filename);
        FileWriter fileWriter = new FileWriter(file);

        // This block of code below was relevant when calculating pressures.
        fileWriter.write("temperature" + "," + "N" + "," + "truncationRadius" + "," + "rho" +
                "," + "Q0" + "," + "Q1" + "pressureConventional" + "," + "pressureConventionalError" + "," +
                "pressureHMANormalMode" + "," + "pressureHMANormalModeError" + ","
                + "pressureHMARealSpace" + "," + "pressureHMARealSpaceError" + "\n");

        // This block of code is relevant when calculating the energy.
//        fileWriter.write("temperature" + "," + "N" + "," + "truncationRadius" + "," + "rho" +
//                "," + "Q1" + "," + "Q2" + "," +
//                "Energy" + "," + "Energy - Harmonic Contribution" + ","
//                + "Energy Error" + "," + "Energy - Harmonic Contribution Error" + "\n");

        for(double temperature : temperatureList) {
            for(int N : systemSizeList) {
                for(double truncationRadius : truncationRadiusList) {
                    for(double rho : rhoList) {

                        double latticeConstant = 1 / rho;

                        // This 'if' statement will for nearest-neighbors less than 1. Be sure to cast to 'int'.
                        double nearestNeighbors = truncationRadius / latticeConstant;
                        if (nearestNeighbors <= 1.001)
                            continue;

                        LJMC1D_NormalModes sim = new LJMC1D_NormalModes(space, temperature, N, truncationRadius, rho);

                        // Equilibration first.
                        sim.activityIntegrate.setMaxSteps(steps / 10);
                        sim.getController().isRunningActivityStep();
                        sim.getController().completeActivities();
                        sim.getController().restartCurrentActivity();
                        // sim.getController().actionPerformed();
                        // sim.getController().reset();

                        MeterPressureHMA pMeterHMA = new MeterPressureHMA(space, sim.potentialMaster, sim.coordinates, true, 0);
                        pMeterHMA.setTemperature(sim.integrator.getTemperature());
                        pMeterHMA.setTruncationRadius(truncationRadius);
                        pMeterHMA.setPRes();        // This should probably be done automatically.
                        pMeterHMA.calculateGij();

                        final AccumulatorAverageFixed pAccumulatorConventional = new AccumulatorAverageFixed();
                        final AccumulatorAverageFixed pAccumulatorHMARealSpace = new AccumulatorAverageFixed();
                        final AccumulatorAverageFixed pAccumulatorHMANormalMode = new AccumulatorAverageFixed();

                        DataSplitter splitter = new DataSplitter();
                        final DataPumpListener pPump = new DataPumpListener(pMeterHMA, splitter, N);
                        splitter.setDataSink(1, pAccumulatorConventional);
                        splitter.setDataSink(4, pAccumulatorHMANormalMode);
                        splitter.setDataSink(5, pAccumulatorHMARealSpace);

                        sim.integrator.getEventManager().addListener(pPump);

                        // Production code.
                        sim.activityIntegrate.setMaxSteps(steps);
                        sim.getController().isRunningActivityStep();
                        sim.getController().completeActivities();
                        // sim.getController().actionPerformed();

                        // TODO: Is the standard deviation of the pressure more useful, compared to the error?
                        int pressureIndex = pAccumulatorConventional.AVERAGE.index;      // index 1 represents the pressure's average.
                        int errorIndex = pAccumulatorConventional.ERROR.index;          // index 2 represents the error in the pressure's average.

                        double pressureConventional = pAccumulatorConventional.getData().getValue(pressureIndex);
                        double pressureConventionalError = pAccumulatorConventional.getData().getValue(errorIndex);

                        double pressureHMANormalMode = pAccumulatorHMANormalMode.getData().getValue(pressureIndex);
                        double pressureHMANormalModeError = pAccumulatorHMANormalMode.getData().getValue(errorIndex);

                        double pressureHMARealSpace = pAccumulatorHMARealSpace.getData().getValue(pressureIndex);
                        double pressureHMARealSpaceError = pAccumulatorHMARealSpace.getData().getValue(errorIndex);

//                        double energy = sim.mcMoveMode.getEnergy();
                        double[] q = sim.mcMoveMode.getModes();

                        // This will write the normal-mode coordinates
                        fileWriter.write(temperature + "," + N + "," + truncationRadius + "," + rho +
                                "," + q[1] + "," + q[2] + "," + pressureConventional + "," + pressureConventionalError + "," +
                                pressureHMANormalMode + "," + pressureHMANormalModeError + ","
                                + pressureHMARealSpace + "," + pressureHMARealSpaceError + "\n");

                        System.out.println(pressureConventional + ", " +
                                pressureHMANormalMode + ", " + pressureHMARealSpace);

                        fileWriter.flush();
                    }
                }
            }
        }
        fileWriter.close();
    }
}
