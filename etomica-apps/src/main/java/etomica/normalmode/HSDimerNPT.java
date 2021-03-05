/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.BoxInflateAnisotropic;
import etomica.action.BoxInflateDeformable;
import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.histogram.HistogramExpanding;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.meter.MeterVolume;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.*;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcpBaseCentered;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveMonoclinic;
import etomica.nbr.list.NeighborListManagerSlanty;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.CoordinateDefinitionHSDimer.IntegerFunction;
import etomica.potential.P2HardSphere;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.dimensions.Null;
import etomica.util.IEvent;
import etomica.util.IListener;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * NPT simulation for hard sphere solid using an MCMove that does coordinate
 * scaling with volume changes.
 */
public class HSDimerNPT extends Simulation {

    public static boolean doFancyMonoclinic = false;
    public static boolean doFancyOriented = true;
    public static boolean doRot = true;
    public static boolean doShape = true;
    public static boolean doGraphics = true;
    public final PotentialMasterList potentialMaster;
    public final IntegratorMC integrator;
    public final SpeciesGeneral species;
    public final Box box, latticeBox;
    public final CoordinateDefinitionHSDimer coordinateDefinition;
    public final double theta;

    public HSDimerNPT(final Space space, int numMolecules, boolean fancyMove, double pSet, double rho, int[] nC, int cp, double L, double thetaFrac, double targetAcc) {
        super(space);

        species = SpeciesHSDimer.create(space, true, L);
        addSpecies(species);

        potentialMaster = new PotentialMasterList(this, 2, new NeighborListManagerSlanty.NeighborListSlantyAgentSource(2), space);

        double tol = 1e-8;
        double a = Math.sqrt(3.0) + tol;
        double b = 1.0 + tol;
        double cx = (cp == 1 ? (L * L + 1) / Math.sqrt(3) : (L * L - 1) / Math.sqrt(3)) + tol;
        double cz = L * Math.sqrt(1 - L * L / 3.0) + Math.sqrt(2.0 / 3.0) + tol;
        double c = Math.sqrt(cx * cx + cz * cz);
        double boxAngle = Math.atan2(cz, cx);

        System.out.println("a: " + a);
        System.out.println("b: " + b);
        System.out.println("c: " + c + " (" + cx + ",0," + cz + ")");

        theta = Math.asin(L / Math.sqrt(3.0));
        System.out.println("close-packed theta: " + theta);

        double sigma = 1.0;
        Vector[] boxDim = new Vector[3];

        Basis unitBasis;
        if (cp == 1 || cp == 2) {
            boxDim[0] = Vector.of(new double[]{nC[0] * a, 0.0, 0.0});
            boxDim[1] = Vector.of(new double[]{0.0, nC[1] * b, 0.0});
            boxDim[2] = Vector.of(new double[]{nC[2] * cx, 0.0, nC[2] * cz});
            unitBasis = new BasisHcpBaseCentered();
        } else {
            throw new RuntimeException("not yet");
        }
        latticeBox = this.makeBox(new BoundaryDeformablePeriodic(space, boxDim));
        integrator = new IntegratorMC(potentialMaster, getRandom(), 1.0, latticeBox);
        this.getController().addActivity(new ActivityIntegrate(integrator));

        P2HardSphere p2 = new P2HardSphere(space, sigma, false);
        potentialMaster.addPotential(p2, new AtomType[]{species.getLeafType(), species.getLeafType()});

        Boundary boundary = new BoundaryDeformablePeriodic(space, boxDim);
        box = this.makeBox(boundary);
        box.setNMolecules(species, numMolecules);

        Primitive primitive = new PrimitiveMonoclinic(space, nC[0] * a, nC[1] * b, nC[2] * c,
                boxAngle);

        Basis basis = new BasisBigCell(space, unitBasis, new int[]{nC[0], nC[1], nC[2]});

        coordinateDefinition = new CoordinateDefinitionHSDimer(getSpeciesManager(), box, primitive, basis, space);
        Vector[][] axes = new Vector[1][3];
        int[] iaxis = new int[]{2, 0, 1};
        for (int i = 0; i < 3; i++) {
            axes[0][i] = space.makeVector();
            axes[0][i].setX(iaxis[i], 1);
        }
        IntegerFunction selector = new IntegerFunction() {
            public int f(int i) {
                return 0;
            }
        };
        coordinateDefinition.setOrientations(axes, new double[]{theta}, selector);
        coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});
        latticeBox.setNMolecules(species, numMolecules);
        CoordinateDefinitionHSDimer coordinateDefinitionLattice = new CoordinateDefinitionHSDimer(getSpeciesManager(), latticeBox, primitive, basis, space);
        coordinateDefinitionLattice.setOrientations(axes, new double[]{theta}, selector);
        coordinateDefinitionLattice.initializeCoordinates(new int[]{1, 1, 1});

        potentialMaster.getNeighborManager(box).reset();

        MCMoveMoleculeCoupled mcMove = new MCMoveMoleculeCoupled(potentialMaster, getRandom(), space);
        mcMove.setBox(box);
        mcMove.setDoExcludeNonNeighbors(true);
        integrator.getMoveManager().addMCMove(mcMove);

        Vector[] drSum = new Vector[nC[2]];
        for (int i = 0; i < drSum.length; i++) {
            drSum[i] = space.makeVector();
        }
        double maxPhi = Math.PI / 6;
        if (doRot) {
            MCMoveRotateMolecule3D rotate = new MCMoveRotateMolecule3D(potentialMaster, getRandom(), space);
            rotate.setBox(box);
//            integrator.getMoveManager().addMCMove(rotate);
            MCMoveRotateMoleculePhiTheta rotateTheta = new MCMoveRotateMoleculePhiTheta(potentialMaster, random, space, drSum, false);
            rotateTheta.setMaxPhi(maxPhi);
            integrator.getMoveManager().addMCMove(rotateTheta);
//            ((MCMoveStepTracker)rotateTheta.getTracker()).setNoisyAdjustment(true);
            MCMoveRotateMoleculePhiTheta rotatePhi = new MCMoveRotateMoleculePhiTheta(potentialMaster, random, space, drSum, true);
            rotatePhi.setMaxPhi(maxPhi);
            integrator.getMoveManager().addMCMove(rotatePhi);
            integrator.getMoveManager().setFrequency(rotateTheta, 0.5);
            integrator.getMoveManager().setFrequency(rotatePhi, 0.5);
//            ((MCMoveStepTracker)rotatePhi.getTracker()).setNoisyAdjustment(true);
        }

        double d3 = 1.0 + L * (1.5 - 0.5 * L * L); //1.792 for L=0.6;
        System.out.println("d3 " + d3);

        if (rho != 0) {
            double initRho = Math.abs(rho) / d3;
            BoxInflateDeformable inflater = new BoxInflateDeformable(box, space);
            inflater.setTargetDensity(initRho);
            inflater.actionPerformed();
            inflater.setBox(latticeBox);
            inflater.actionPerformed();
        }

        if (rho > 0) {
//            MCMoveVolumeMonoclinic mcMoveVolMonoclinic = new MCMoveVolumeMonoclinic(potentialMaster, getRandom(), space);
//            mcMoveVolMonoclinic.setBox(box);
//            mcMoveVolMonoclinic.setStepSize(0.001);
//            ((MCMoveStepTracker)mcMoveVolMonoclinic.getTracker()).setNoisyAdjustment(true);
//            integrator.getMoveManager().addMCMove(mcMoveVolMonoclinic);
//
//            MCMoveVolumeMonoclinicAngle mcMoveVolMonoclinicAngle = new MCMoveVolumeMonoclinicAngle(potentialMaster, getRandom(), space, box);
//            mcMoveVolMonoclinicAngle.setBox(box);
//            mcMoveVolMonoclinicAngle.setStepSize(0.001);
//            ((MCMoveStepTracker)mcMoveVolMonoclinicAngle.getTracker()).setNoisyAdjustment(true);
//            integrator.getMoveManager().addMCMove(mcMoveVolMonoclinicAngle);

        } else {
            double p = pSet / d3;
            final BoxInflateDeformable inflater = new BoxInflateDeformable(box, space);
            final BoxInflateDeformable inflaterLat = new BoxInflateDeformable(latticeBox, space);
            MCMoveBoxStep mcMoveVolume;
            if (fancyMove) {
                // fancy move
                if (false) {
                    mcMoveVolume = new MCMoveVolumeSolidNPTMolecular(potentialMaster, getRandom(), space, p);
                } else {
                    if (!doRot) throw new RuntimeException("oops");
                    mcMoveVolume = new MCMoveVolumeSolidNPTMolecularOriented(potentialMaster, getRandom(), space, p, drSum);
                    ((MCMoveVolumeSolidNPTMolecularOriented) mcMoveVolume).setNominalTheta(theta);
                    ((MCMoveVolumeSolidNPTMolecularOriented) mcMoveVolume).setMaxPhi(maxPhi);
                    ((MCMoveVolumeSolidNPTMolecularOriented) mcMoveVolume).setThetaFrac(thetaFrac);
                }
                ((MCMoveVolumeSolidNPTMolecular) mcMoveVolume).setTemperature(1.0);
                ((MCMoveVolumeSolidNPTMolecular) mcMoveVolume).setInflater(inflaterLat);
                ((MCMoveVolumeSolidNPTMolecular) mcMoveVolume).setLatticeBox(latticeBox);

                if (doShape) {
                    if (doFancyMonoclinic) {

                        MCMoveVolumeMonoclinicScaled mcMoveVolMonoclinic = new MCMoveVolumeMonoclinicScaled(potentialMaster, getRandom(), space, p);
                        mcMoveVolMonoclinic.setTemperature(1.0);
                        mcMoveVolMonoclinic.setInflater(inflaterLat);
                        mcMoveVolMonoclinic.setLatticeBox(latticeBox);
                        mcMoveVolMonoclinic.setStepSize(0.001);
                        ((MCMoveStepTracker) mcMoveVolMonoclinic.getTracker()).setNoisyAdjustment(true);
                        integrator.getMoveManager().addMCMove(mcMoveVolMonoclinic);
                    } else {
                        MCMoveVolumeMonoclinic mcMoveVolMonoclinic = new MCMoveVolumeMonoclinic(potentialMaster, getRandom(), space);
                        mcMoveVolMonoclinic.setInflater(inflater);
                        mcMoveVolMonoclinic.setStepSize(0.001);
                        ((MCMoveStepTracker) mcMoveVolMonoclinic.getTracker()).setNoisyAdjustment(true);
                        integrator.getMoveManager().addMCMove(mcMoveVolMonoclinic);
                    }


                    final BoxInflateAnisotropic inflateAngle = new BoxInflateAnisotropic(latticeBox, space);
                    integrator.getMoveEventManager().addListener(new IListener() {
                        public void actionPerformed(IEvent event) {
                            if (event instanceof MCMoveTrialCompletedEvent && ((MCMoveTrialCompletedEvent) event).isAccepted()) {
                                Vector scaleVec = space.makeVector();
                                MCMove move = ((MCMoveTrialCompletedEvent) event).getMCMove();
                                if (move instanceof MCMoveVolumeMonoclinic) {
                                    for (int i = 0; i < 3; i++) {
                                        Vector lbv = latticeBox.getBoundary().getEdgeVector(i);
                                        Vector bv = box.getBoundary().getEdgeVector(i);
                                        scaleVec.setX(i, bv.getX(i) / lbv.getX(i));
                                    }
                                    inflaterLat.setVectorScale(scaleVec);
                                    inflaterLat.actionPerformed();
                                } else if (move instanceof MCMoveVolumeMonoclinicAngle) {
                                    scaleVec.E(box.getBoundary().getEdgeVector(2));
                                    inflateAngle.setCVector(scaleVec);
                                    inflateAngle.actionPerformed();
                                }
                            }
                        }
                    });
                }

            } else {
                // standard move
                mcMoveVolume = new MCMoveVolume(potentialMaster, getRandom(), space, p);
                mcMoveVolume.setStepSize(1e-4);
                ((MCMoveVolume) mcMoveVolume).setInflater(inflater);

                MCMoveVolumeMonoclinic mcMoveVolMonoclinic = new MCMoveVolumeMonoclinic(potentialMaster, getRandom(), space);
                mcMoveVolMonoclinic.setInflater(inflater);
                mcMoveVolMonoclinic.setStepSize(0.0001);
                ((MCMoveStepTracker) mcMoveVolMonoclinic.getTracker()).setNoisyAdjustment(true);
                if (doShape) {
                    integrator.getMoveManager().addMCMove(mcMoveVolMonoclinic);
                }
            }
            ((MCMoveStepTracker) mcMoveVolume.getTracker()).setAcceptanceTarget(targetAcc);
            ((MCMoveStepTracker) mcMoveVolume.getTracker()).setNoisyAdjustment(true);
            integrator.getMoveManager().addMCMove(mcMoveVolume);

            MCMoveVolumeMonoclinicAngle mcMoveVolMonoclinicAngle = new MCMoveVolumeMonoclinicAngle(potentialMaster, getRandom(), space, box);
            mcMoveVolMonoclinicAngle.setBox(box);
            mcMoveVolMonoclinicAngle.setStepSize(0.001);
            ((MCMoveStepTracker) mcMoveVolMonoclinicAngle.getTracker()).setNoisyAdjustment(true);
            if (doShape) {
                integrator.getMoveManager().addMCMove(mcMoveVolMonoclinicAngle);
            }
        }

    }

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        final HSMD3DParameters params = new HSMD3DParameters();
        if (args.length > 0) {
            ParseArgs parseArgs = new ParseArgs(params);
            parseArgs.parseArgs(args);
        }

        int[] nC = params.nC;
        int numMolecules = nC[0] * nC[1] * nC[2] * 2;
        double rho0 = params.rho;
        if (params.pressure > 0) System.out.println("pressure: " + params.pressure);
        if (rho0 < 0) {
            if (params.pressure == 45) {
                rho0 = -1.3;
            } else if (params.pressure == 100) {
                rho0 = -1.4;
            }
            System.out.println("initial density: " + (-rho0));
        } else {
            System.out.println("density: " + rho0);
        }
        final HSDimerNPT sim = new HSDimerNPT(Space3D.getInstance(), numMolecules, params.fancyMove,
                params.pressure, rho0, params.nC, params.cp, params.L, params.thetaFrac, params.targetAcc);

        final MeterVolume meterVolume = new MeterVolume();
        meterVolume.setBox(sim.box);
        DataFork volumeFork = new DataFork();
        DataPumpListener volumePump = new DataPumpListener(meterVolume, volumeFork, numMolecules);
        sim.integrator.getEventManager().addListener(volumePump);
        long nData = params.numSteps / numMolecules;
        long blockSize = (nData + 99) / 100;
        final AccumulatorAverage volumeAvg = doGraphics ? new AccumulatorAverageCollapsing() : new AccumulatorAverageFixed(blockSize);
        if (doGraphics) {
            volumeAvg.setIncludeACInError(true);
        }
        volumeFork.addDataSink(volumeAvg);

        Vector a0 = sim.box.getBoundary().getEdgeVector(0);
        Vector b0 = sim.box.getBoundary().getEdgeVector(1);
        Vector c0 = sim.box.getBoundary().getEdgeVector(2);
        final double yx0 = Math.sqrt(b0.squared() / a0.squared());
        final double zx0 = Math.sqrt(c0.squared() / a0.squared());
        DataSourceScalar meterXY = new DataSourceScalar("xy", Null.DIMENSION) {
            public double getDataAsScalar() {
                Vector a = sim.box.getBoundary().getEdgeVector(0);
                Vector b = sim.box.getBoundary().getEdgeVector(1);
                return (Math.sqrt(b.squared() / a.squared()) / yx0) - 1;
            }
        };
        DataPumpListener xyPump = new DataPumpListener(meterXY, null, numMolecules);
        sim.integrator.getEventManager().addListener(xyPump);

        DataSourceScalar meterXZ = new DataSourceScalar("xy", Null.DIMENSION) {
            public double getDataAsScalar() {
                Vector a = sim.box.getBoundary().getEdgeVector(0);
                Vector c = sim.box.getBoundary().getEdgeVector(2);
                return (Math.sqrt(c.squared() / a.squared()) / zx0) - 1;
            }
        };
        DataPumpListener xzPump = new DataPumpListener(meterXZ, null, numMolecules);
        sim.integrator.getEventManager().addListener(xzPump);

        final double cxcz0 = c0.getX(0) / c0.getX(2);
        final DataSourceScalar meterCXCZ = new DataSourceScalar("cx/cz", Null.DIMENSION) {
            public double getDataAsScalar() {
                Vector c = sim.box.getBoundary().getEdgeVector(2);
                return (c.getX(0) / c.getX(2) / cxcz0) - 1;
            }
        };
        DataPumpListener cxczPump = new DataPumpListener(meterCXCZ, null, numMolecules);
        sim.integrator.getEventManager().addListener(cxczPump);

        MeterTilt meterTilt = new MeterTilt(sim.space, sim.species, nC[2]);
        meterTilt.setBox(sim.box);
        DataFork tiltFork = new DataFork();
        DataPumpListener tiltPump = new DataPumpListener(meterTilt, tiltFork, numMolecules);
        sim.integrator.getEventManager().addListener(tiltPump);
        final AccumulatorAverage tiltAvg = new AccumulatorAverageFixed(blockSize);
        tiltFork.addDataSink(tiltAvg);

        MeterPhiDeviation meterPhiDeviation = new MeterPhiDeviation(sim.space);
        meterPhiDeviation.setBox(sim.box);
        DataFork phiDeviationFork = new DataFork();
        DataPumpListener phiDeviationPump = new DataPumpListener(meterPhiDeviation, phiDeviationFork, numMolecules);
        if (params.rho > 0) {
            sim.integrator.getEventManager().addListener(phiDeviationPump);
        }
        final AccumulatorAverage phiDeviationAvg = doGraphics ? new AccumulatorAverageCollapsing() : new AccumulatorAverageFixed(blockSize);
        phiDeviationFork.addDataSink(phiDeviationAvg);

        MeterThetaDeviation meterThetaDeviation = new MeterThetaDeviation(sim.space);
        meterThetaDeviation.setBox(sim.box);
        meterThetaDeviation.setNominalTheta(sim.theta);
        DataFork thetaDeviationFork = new DataFork();
        DataPumpListener thetaDeviationPump = new DataPumpListener(meterThetaDeviation, thetaDeviationFork, numMolecules);
        if (params.rho > 0) {
            sim.integrator.getEventManager().addListener(thetaDeviationPump);
        }
        final AccumulatorAverage thetaDeviationAvg = doGraphics ? new AccumulatorAverageCollapsing() : new AccumulatorAverageFixed(blockSize);
        thetaDeviationFork.addDataSink(thetaDeviationAvg);


//      MeterDisplacement meterDisplacement = new MeterDisplacement(sim.space, sim.coordinateDefinition, 0.001);
//      sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterDisplacement, params.numAtoms));

        MeterDisplacementMoleculeRMS meterDisplacementRMS = new MeterDisplacementMoleculeRMS(sim.space, sim.coordinateDefinition);
        final AccumulatorAverage displacementAvg = doGraphics ? new AccumulatorAverageCollapsing() : new AccumulatorAverageFixed(blockSize);
        DataPumpListener displacementAvgPump = new DataPumpListener(meterDisplacementRMS, displacementAvg, numMolecules);
        sim.integrator.getEventManager().addListener(displacementAvgPump);
//
//      MeterDisplacementMax meterDisplacementMax = new MeterDisplacementMax(sim.space, sim.coordinateDefinition, 0.001);
//      AccumulatorAverageFixed displacementMax = new AccumulatorAverageFixed();
//      DataPumpListener displacementMaxPump = new DataPumpListener(meterDisplacementMax, displacementMax, params.numAtoms);
//      sim.integrator.getEventManager().addListener(displacementMaxPump);
//
//      MeterMaxExpansion meterMaxExpansion = new MeterMaxExpansion(sim.space, sim.box, sim.potentialMaster.getNeighborManager(sim.box));
//      AccumulatorAverageFixed maxExpansionAvg = new AccumulatorAverageFixed();
//      DataPumpListener maxExpansionPump = new DataPumpListener(meterMaxExpansion, maxExpansionAvg, params.numAtoms);
//      sim.integrator.getEventManager().addListener(maxExpansionPump);

        if (doGraphics) {

            final MeterTiltRotation meterTiltRotation = new MeterTiltRotation(sim.space, sim.species, nC[2]);
            meterTiltRotation.setBox(sim.box);
            DataFork tiltRotationFork = new DataFork();
            DataPumpListener tiltRotationPump = new DataPumpListener(meterTiltRotation, tiltRotationFork, numMolecules);
            sim.integrator.getEventManager().addListener(tiltRotationPump);
            AccumulatorAverageFixed tiltRotationAvg = new AccumulatorAverageFixed();
            tiltRotationFork.addDataSink(tiltRotationAvg);

            final MeterTiltHistogram meterTiltHistogram = new MeterTiltHistogram(sim.space, sim.species);
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterTiltHistogram, numMolecules));
            meterTiltHistogram.setBox(sim.box);
            DataFork tiltHistogramFork = new DataFork();
            DataPumpListener tiltHistogramPump = new DataPumpListener(meterTiltHistogram, tiltHistogramFork, numMolecules);
            sim.integrator.getEventManager().addListener(tiltHistogramPump);

            final MeterTiltRotationHistogram meterRotationHistogram = new MeterTiltRotationHistogram(sim.space, sim.species, nC[2]);
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterRotationHistogram, numMolecules));
            meterRotationHistogram.setBox(sim.box);
            DataFork tiltRotationHistogramFork = new DataFork();
            DataPumpListener tiltRotationHistogramPump = new DataPumpListener(meterRotationHistogram, tiltRotationHistogramFork, 100 * numMolecules);
            sim.integrator.getEventManager().addListener(tiltRotationHistogramPump);

            MeterTiltRotationStdev meterTiltRotationStdev = new MeterTiltRotationStdev(sim.space, sim.species, nC[2]);
            meterTiltRotationStdev.setBox(sim.box);
            DataFork tiltRotationStdevFork = new DataFork();
            DataPumpListener tiltRotationStdevPump = new DataPumpListener(meterTiltRotationStdev, tiltRotationStdevFork, numMolecules);
            sim.integrator.getEventManager().addListener(tiltRotationStdevPump);
            AccumulatorAverageFixed tiltRotationStdevAvg = new AccumulatorAverageFixed();
            tiltRotationStdevFork.addDataSink(tiltRotationStdevAvg);


            MeterPlaneSlip meterPlaneSlip = new MeterPlaneSlip(sim.space, sim.species, nC[2], nC[0], nC[1]);
            meterPlaneSlip.setBox(sim.box);
            DataFork planeSlipFork = new DataFork();
            DataPumpListener planeSlipPump = new DataPumpListener(meterPlaneSlip, planeSlipFork, numMolecules);
            sim.integrator.getEventManager().addListener(planeSlipPump);


            if (false) {
                String filename = "vol";
                filename += params.fancyMove ? "_fancy" : "_simple";
                if (doFancyMonoclinic) filename += "mono";
                if (doFancyOriented) filename += "orient";
                filename += "foo2";
                filename += ".dat";
                final File wvFile = new File(filename);
                if (wvFile.exists()) wvFile.delete();
                IAction wv = new IAction() {
                    public void actionPerformed() {
                        try {
                            FileWriter fw = new FileWriter(wvFile, true);
                            double v = meterVolume.getDataAsScalar();
                            double cxcz = meterCXCZ.getDataAsScalar();
                            IData planeRot = meterTiltRotation.getData();
                            long step = sim.integrator.getStepCount();
                            fw.write(step + "  " + v + "  " + cxcz + " ");
                            for (int i = 0; i < planeRot.getLength(); i++) {
                                fw.write(" " + planeRot.getValue(i));
                            }
                            fw.write("\n");
                            fw.close();
                        } catch (IOException ex) {
                            throw new RuntimeException(ex);
                        }
                    }
                };
                int wvInterval = numMolecules;
                wvInterval *= 10;
                if (!params.fancyMove) {
                    wvInterval *= 10;
                }

                sim.integrator.getEventManager().addListener(new IntegratorListenerAction(wv, wvInterval));
            }


            SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DiameterHashByType diameter = new DiameterHashByType();
            diameter.setDiameter(sim.species.getLeafType(), 1.0);
            graphic.getDisplayBox(sim.box).setDiameterHash(diameter);

//			ColorSchemeNeighbor colorScheme = new ColorSchemeNeighbor(sim, sim.potentialMaster, sim.box);
//			colorScheme.setAtom(sim.box.getLeafList().getAtom(122));
//			graphic.getDisplayBox(sim.box).setColorScheme(colorScheme);

            DataSourceCountSteps stepCounter = new DataSourceCountSteps(sim.integrator);

            AccumulatorHistory volumeHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
            volumeHistory.setPushInterval(100);
            volumeHistory.setTimeDataSource(stepCounter);
            volumeFork.addDataSink(volumeHistory);
            DisplayPlot volumePlot = new DisplayPlot();
            volumeHistory.setDataSink(volumePlot.getDataSet().makeDataSink());
            volumePlot.setDoLegend(false);
            volumePlot.setLabel("volume");
            graphic.add(volumePlot);

            final AccumulatorHistogram volumeHistogram = new AccumulatorHistogram(new HistogramExpanding(0.05));
            if (params.rho <= 0) {
                volumeFork.addDataSink(volumeHistogram);
                volumeHistogram.setPushInterval(100);
                DisplayPlot volumeHistogramPlot = new DisplayPlot();
                volumeHistogram.setDataSink(volumeHistogramPlot.getDataSet().makeDataSink());
                volumeHistogramPlot.setLabel("v histogram");
                volumeHistogramPlot.setDoLegend(false);
                graphic.add(volumeHistogramPlot);

                AccumulatorHistory volumeAvgHistory = new AccumulatorHistory(new HistoryCollapsingDiscard());
                volumeAvgHistory.setTimeDataSource(stepCounter);
                volumeAvg.addDataSink(volumeAvgHistory, new StatType[]{volumeAvg.AVERAGE});
                volumeAvgHistory.setDataSink(volumePlot.getDataSet().makeDataSink());
                AccumulatorHistory volumeAvgPHistory = new AccumulatorHistory(new HistoryCollapsingDiscard());
                volumeAvgPHistory.setTimeDataSource(stepCounter);
                DataProcessor volumeAvgP = new DataProcessor() {

                    DataInfoDouble dataInfo = new DataInfoDouble("foo", Null.DIMENSION);
                    DataDouble data = new DataDouble();

                    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
                        return dataInfo;
                    }

                    protected IData processData(IData inputData) {
                        data.x = inputData.getValue(0) + inputData.getValue(1);
                        return data;
                    }
                };
                volumeAvg.addDataSink(volumeAvgP, new StatType[]{volumeAvg.AVERAGE, volumeAvg.ERROR});
                volumeAvgP.setDataSink(volumeAvgPHistory);
                volumeAvgPHistory.setDataSink(volumePlot.getDataSet().makeDataSink());

                AccumulatorHistory volumeStdevHistory = new AccumulatorHistory(new HistoryCollapsingDiscard());
                volumeStdevHistory.setTimeDataSource(stepCounter);
                volumeAvg.addDataSink(volumeStdevHistory, new StatType[]{volumeAvg.STANDARD_DEVIATION});
                DisplayPlot volumeStdevPlot = new DisplayPlot();
                volumeStdevHistory.setDataSink(volumeStdevPlot.getDataSet().makeDataSink());
                volumeStdevPlot.setLabel("volume stdev");
                volumeStdevPlot.setDoLegend(false);
                graphic.add(volumeStdevPlot);
            }


            DataSplitter tiltSplitter = new DataSplitter();
            tiltFork.addDataSink(tiltSplitter);
            DisplayPlot tiltPlot = new DisplayPlot();
            for (int i = 0; i < nC[2] + 1; i++) {
                AccumulatorHistory tiltHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
                tiltHistory.setPushInterval(100);
                tiltHistory.setTimeDataSource(stepCounter);
                tiltSplitter.setDataSink(i, tiltHistory);
                tiltHistory.setDataSink(tiltPlot.getDataSet().makeDataSink());
                tiltPlot.setLegend(new DataTag[]{tiltHistory.getTag()}, (i == 0 ? "total" : "" + (i - 1)));
            }
            tiltPlot.setLabel("tilt");
            graphic.add(tiltPlot);

            DisplayPlot tiltHistogramPlot = new DisplayPlot();
            tiltHistogramFork.addDataSink(tiltHistogramPlot.getDataSet().makeDataSink());
            tiltHistogramPlot.setDoLegend(false);
            tiltHistogramPlot.setLabel("tilt histogram");
            graphic.add(tiltHistogramPlot);


            DataSplitter tiltRotationSplitter = new DataSplitter();
            tiltRotationFork.addDataSink(tiltRotationSplitter);
            DisplayPlot tiltRotationPlot = new DisplayPlot();
            for (int i = 0; i < nC[2] + 1; i++) {
                AccumulatorHistory tiltRotationHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
                tiltRotationHistory.setPushInterval(100);
                tiltRotationHistory.setTimeDataSource(stepCounter);
                tiltRotationSplitter.setDataSink(i, tiltRotationHistory);
                tiltRotationHistory.setDataSink(tiltRotationPlot.getDataSet().makeDataSink());
                tiltRotationPlot.setLegend(new DataTag[]{tiltRotationHistory.getTag()}, (i == 0 ? "total" : "" + (i - 1)));
            }
            tiltRotationPlot.setLabel("rotation");
            graphic.add(tiltRotationPlot);

            DataSplitter tiltRotationStdevSplitter = new DataSplitter();
            tiltRotationStdevFork.addDataSink(tiltRotationStdevSplitter);
            DisplayPlot tiltRotationStdevPlot = new DisplayPlot();
            for (int i = 0; i < nC[2] + 1; i++) {
                AccumulatorHistory tiltRotationStdevHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
                tiltRotationStdevHistory.setPushInterval(1000);
                tiltRotationStdevHistory.setTimeDataSource(stepCounter);
                tiltRotationStdevSplitter.setDataSink(i, tiltRotationStdevHistory);
                tiltRotationStdevHistory.setDataSink(tiltRotationStdevPlot.getDataSet().makeDataSink());
                tiltRotationStdevPlot.setLegend(new DataTag[]{tiltRotationStdevHistory.getTag()}, (i == 0 ? "total" : "" + (i - 1)));
            }
            tiltRotationStdevPlot.setLabel("rot stdev");
            graphic.add(tiltRotationStdevPlot);

            AccumulatorHistory xyHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
            xyHistory.setPushInterval(100);
            xyHistory.setTimeDataSource(stepCounter);
            xyPump.setDataSink(xyHistory);
            AccumulatorHistory xzHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
            xzHistory.setPushInterval(100);
            xzHistory.setTimeDataSource(stepCounter);
            xzPump.setDataSink(xzHistory);
            AccumulatorHistory cxczHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
            cxczHistory.setPushInterval(100);
            cxczHistory.setTimeDataSource(stepCounter);
            cxczPump.setDataSink(cxczHistory);
            DisplayPlot xyzPlot = new DisplayPlot();
            if (false) {
                xyHistory.setDataSink(xyzPlot.getDataSet().makeDataSink());
                xyzPlot.setLegend(new DataTag[]{xyHistory.getTag()}, "b/a");
            }
            xzHistory.setDataSink(xyzPlot.getDataSet().makeDataSink());
            xyzPlot.setLegend(new DataTag[]{xzHistory.getTag()}, "c/a");
            cxczHistory.setDataSink(xyzPlot.getDataSet().makeDataSink());
            xyzPlot.setLegend(new DataTag[]{cxczHistory.getTag()}, "cx/cz");
            xyzPlot.setLabel("abc");
            graphic.add(xyzPlot);

            DisplayPlot rotationHistogramPlot = new DisplayPlot();
            DataGroupSplitter histogramSplitter = new DataGroupSplitter();
            tiltRotationHistogramFork.addDataSink(histogramSplitter);
            for (int i = 0; i < 2; i++) {
                histogramSplitter.setDataSink(i, rotationHistogramPlot.getDataSet().makeDataSink());
                rotationHistogramPlot.setLegend(new DataTag[]{meterRotationHistogram.getTag(i)}, i == 0 ? "total" : "plane");
            }
            rotationHistogramPlot.setLabel("rot histogram");
            graphic.add(rotationHistogramPlot);

            DataSplitter planeSlipSplitter = new DataSplitter();
            planeSlipFork.addDataSink(planeSlipSplitter);
            DisplayPlot planeSlipPlot = new DisplayPlot();
            for (int i = 0; i < nC[2] + 1; i++) {
                AccumulatorHistory tiltRotationStdevHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
                tiltRotationStdevHistory.setTimeDataSource(stepCounter);
                tiltRotationStdevSplitter.setDataSink(i, tiltRotationStdevHistory);
                tiltRotationStdevHistory.setDataSink(tiltRotationStdevPlot.getDataSet().makeDataSink());
                tiltRotationStdevPlot.setLegend(new DataTag[]{tiltRotationStdevHistory.getTag()}, (i == 0 ? "total" : "" + (i - 1)));
            }
            for (int i = 0; i < nC[2] * 2; i++) {
                AccumulatorHistory planeSlipHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
                planeSlipHistory.setPushInterval(100);
                planeSlipHistory.setTimeDataSource(stepCounter);
                planeSlipSplitter.setDataSink(i, planeSlipHistory);
                planeSlipHistory.setDataSink(planeSlipPlot.getDataSet().makeDataSink());
                planeSlipPlot.setLegend(new DataTag[]{planeSlipHistory.getTag()}, ((i % 2 == 0) ? "x" : "y") + (i / 2));
            }
            planeSlipPlot.setLabel("slip");
            graphic.add(planeSlipPlot);

            graphic.getController().getResetAveragesButton().setPostAction(new IAction() {
                public void actionPerformed() {
                    meterRotationHistogram.reset();
                    meterTiltHistogram.reset();
                    volumeAvg.reset();
                    volumeHistogram.reset();
                }
            });

            graphic.makeAndDisplayFrame();

            return;
        }
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps / 10));
        volumeAvg.reset();
        displacementAvg.reset();
        thetaDeviationAvg.reset();
        phiDeviationAvg.reset();
        System.out.println("equilibration finished");
        sim.integrator.getMoveManager().setEquilibrating(false);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps));

        if (params.rho <= 0) {
            double vavg = volumeAvg.getData().getValue(volumeAvg.AVERAGE.index);
            double verr = volumeAvg.getData().getValue(volumeAvg.ERROR.index);
            double vstdev = volumeAvg.getData().getValue(volumeAvg.STANDARD_DEVIATION.index);
            double vcorr = volumeAvg.getData().getValue(volumeAvg.BLOCK_CORRELATION.index);
            System.out.println("avg volume " + vavg + "  err " + verr + "  stdev " + vstdev + "  correlation " + vcorr);
            System.out.println("avg density " + numMolecules / vavg + " " + numMolecules / (vavg * vavg) * verr);
        } else {
            double davg = displacementAvg.getData().getValue(displacementAvg.AVERAGE.index);
            double dstdev = displacementAvg.getData().getValue(displacementAvg.STANDARD_DEVIATION.index);
            double derr = displacementAvg.getData().getValue(displacementAvg.ERROR.index);
            System.out.println("displacement avg " + davg + " stdev " + dstdev + " err " + derr);

            double thetaavg = thetaDeviationAvg.getData().getValue(thetaDeviationAvg.AVERAGE.index);
            double thetastdev = thetaDeviationAvg.getData().getValue(thetaDeviationAvg.STANDARD_DEVIATION.index);
            double thetaerr = thetaDeviationAvg.getData().getValue(thetaDeviationAvg.ERROR.index);
            System.out.println("cos theta avg " + thetaavg + " stdev " + thetastdev + " err " + thetaerr);

            double phiavg = phiDeviationAvg.getData().getValue(phiDeviationAvg.AVERAGE.index);
            double phistdev = phiDeviationAvg.getData().getValue(phiDeviationAvg.STANDARD_DEVIATION.index);
            double phierr = phiDeviationAvg.getData().getValue(phiDeviationAvg.ERROR.index);
            System.out.println("phi/sintheta avg " + phiavg + " stdev " + phistdev + " err " + phierr);
        }
    }

    public static class HSMD3DParameters extends ParameterBase {
        public int[] nC = new int[]{3, 6, 4};
        public long numSteps = 10000000;
        public double pressure = 45;
        public int cp = 2;
        public double L = 0.6;
        public double rho = -1.3;
        public boolean fancyMove = false;
        public double thetaFrac = 1.0;
        public double targetAcc = 0.5;
    }
}
