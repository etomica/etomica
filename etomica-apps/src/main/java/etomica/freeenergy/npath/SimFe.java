/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.freeenergy.npath;

import etomica.action.BoxInflate;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHash;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.chem.elements.Iron;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.history.HistoryScrolling;
import etomica.data.meter.*;
import etomica.eam.P2EAM;
import etomica.eam.PotentialCalculationEnergySumEAM;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicBcc;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeHcp4;
import etomica.lattice.SpaceLattice;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.lattice.crystal.PrimitiveHCP4;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.PotentialCalculationForceSum;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.*;
import etomica.units.dimensions.Null;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;
import java.util.Arrays;

import static etomica.freeenergy.npath.SimFe.Crystal.HCP;

public class SimFe extends Simulation {
    
    public final PotentialMasterList potentialMaster;
    public IntegratorVelocityVerlet integrator;
    public SpeciesSpheresMono species;
    public Box box;
    public P2EAM potential;
    public P1ImageHarmonic p1ImageHarmonic;
    public MCMoveAtomSwap mcMoveSwap;

    public SimFe(Crystal crystal, int numAtoms, double temperature, double density, double w, int offsetDim, int numInnerSteps, boolean swap, boolean doHarmonic, double timeStep) {
        super(Space3D.getInstance());
        species = new SpeciesSpheresMono(space, Iron.INSTANCE);
        species.setIsDynamic(true);
        addSpecies(species);
        box = this.makeBox();
        box.setNMolecules(species, numAtoms);
        Primitive primitive = (crystal == HCP) ? new PrimitiveHCP4(space) : new PrimitiveCubic(space);
        Vector l = space.makeVector();
        double[] primitiveSize = primitive.getSize();
        int[] f = new int[]{10, 10, 10};
        if (crystal == HCP) f = new int[]{6, 4, 4};
        for (int i = 0; i < 3; i++) {
            double x = f[i] * primitiveSize[i];
            if (i <= offsetDim) x *= 2;
            l.setX(i, x);
        }
        box.getBoundary().setBoxSize(l);

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();

        double n = 8.7932;
        double m = 8.14475;
        double eps = ElectronVolt.UNIT.toSim(0.0220225);
        double a = 3.48501;
        double C = 28.8474;
        double rc = 6;
        potential = new P2EAM(space, n, m, eps, a, C, rc, rc);

        potentialMaster = new PotentialMasterList(this, 1.2 * rc, space);
        potentialMaster.getNbrCellManager(box).setSuppressBoxLengthWarning(true);
        potentialMaster.setCellRange(2);

        AtomType leafType = species.getLeafType();

        potentialMaster.addPotential(potential, new AtomType[]{leafType, leafType});

        Vector offset = space.makeVector();
        offset.setX(offsetDim, box.getBoundary().getBoxSize().getX(offsetDim) * 0.5);
        p1ImageHarmonic = new P1ImageHarmonic(space, offset, w, !doHarmonic);
        potentialMaster.addPotential(p1ImageHarmonic, new AtomType[]{leafType});

        if (numInnerSteps > 0 && w > 0) {
            if (doHarmonic) {
                if (w > 3e6) {
                    integrator = new IntegratorMDHarmonicMC(potentialMaster, random, timeStep, temperature, box);
                    ((IntegratorMDHarmonicMC) integrator).setP1Harmonic(p1ImageHarmonic);
                    swap = false;
                } else {
                    integrator = new IntegratorImageHarmonicMD(potentialMaster, random, timeStep, temperature, box);
                    ((IntegratorImageHarmonicMD) integrator).setP1Harmonic(p1ImageHarmonic);
                }
            } else {
                integrator = new IntegratorImageMultistepMD(potentialMaster, random, timeStep, temperature, box);
                ((IntegratorImageMultistepMD) integrator).setP1Harmonic(p1ImageHarmonic);
                ((IntegratorImageMultistepMD) integrator).setNumInnerSteps(numInnerSteps);
            }
        } else {
            integrator = new IntegratorVelocityVerlet(potentialMaster, random, timeStep, temperature, box);
        }
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster, box);
        meterPE.setPotentialCalculation(new PotentialCalculationEnergySumEAM(potential));
        integrator.setMeterPotentialEnergy(meterPE);
        integrator.setIsothermal(true);
        integrator.setTemperature(temperature);
        integrator.setThermostatNoDrift(true);
        integrator.getEventManager().addListener(potential.makeIntegratorListener(potentialMaster, box));
        integrator.setForceSum(new PotentialCalculationForceSum());

        this.getController().addActivity(new ActivityIntegrate(integrator));

        p1ImageHarmonic.setZeroForce(true);

        SpaceLattice lat = null;
        if (crystal == Crystal.FCC) {
            lat = new LatticeCubicFcc(space);
        } else if (crystal == Crystal.BCC) {
            lat = new LatticeCubicBcc(space);
        } else if (crystal == HCP) {
            lat = new LatticeHcp4(space);
        } else {
            throw new RuntimeException("Don't know how to do " + crystal);
        }
        ConfigurationLattice config = new ConfigurationLattice(lat, space);
        config.initializeCoordinates(box);

        p1ImageHarmonic.findNOffset(box);

        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));
        potentialMaster.getNeighborManager(box).reset();

        Vector boxLength = box.getBoundary().getBoxSize();
        double lMin = boxLength.getX(0);
        if (boxLength.getX(1) < lMin) lMin = boxLength.getX(1);
        if (boxLength.getX(2) < lMin) lMin = boxLength.getX(2);
        double ww = w / lMin;
        double swapDistance = 1.5 * Math.sqrt(1.5 * temperature / ww);
        if (swapDistance > lMin / 4) swapDistance = lMin / 4;
        if (swapDistance > rc) swapDistance = rc;
        if (swapDistance < 2) swapDistance = 2;
        PotentialMasterCell potentialMasterCell = new PotentialMasterCell(this, swapDistance, space);
        potentialMasterCell.setCellRange(2);
        potentialMasterCell.getNbrCellManager(box).assignCellAll();
        if (swap) {
            mcMoveSwap = new MCMoveAtomSwap(random, potentialMasterCell, space, p1ImageHarmonic);
            mcMoveSwap.setNbrDistance(swapDistance);
            IntegratorMC integratorMC = new IntegratorMC(potentialMaster, random, temperature, box);
            integratorMC.getMoveManager().addMCMove(mcMoveSwap);
            integratorMC.reset();

            integrator.getEventManager().addListener(new IntegratorListener() {
                int countdown = 10, interval = 10;

                public void integratorInitialized(IntegratorEvent e) {
                }

                public void integratorStepStarted(IntegratorEvent e) {
                }

                @Override
                public void integratorStepFinished(IntegratorEvent e) {
                    if (--countdown > 0) {
                        return;
                    }
                    countdown = interval;
                    potentialMasterCell.getNbrCellManager(box).assignCellAll();
                    for (int i = 0; i < 10 * numAtoms; i++) {
                        integratorMC.doStep();
                    }
                    integrator.reset();
                }
            });
        }
    }
    
    public static void main(String[] args) {

        LjMC3DParams params = new LjMC3DParams();
        ParseArgs.doParseArgs(params, args);
        if (args.length==0) {
            params.graphics = false;
            params.numAtoms = 1024;
            params.steps = 10000;
            params.density = 0.150375939849624;
            params.T = 7947.1916;
            params.w = 8048.485777;
            params.crystal = Crystal.BCC;
            params.offsetDim = 2;
            params.numInnerSteps = 100;
            params.swap = true;
            params.nve = false;
            params.rmsdFile = "rmsd.dat";
            params.thermostatInterval = 500;
        }

        final int numAtoms = params.numAtoms;
        final double temperatureK = params.T;
        final double temperature = Kelvin.UNIT.toSim(temperatureK);
        final double density = params.density;
        long steps = params.steps;
        boolean graphics = params.graphics;
        double w = params.w;
        int offsetDim = params.offsetDim;
        Crystal crystal = params.crystal;
        boolean swap = params.swap;
        int numInnerSteps = w > 0 ? params.numInnerSteps : 0;
        boolean nve = params.nve;
        int thermostatInterval = params.thermostatInterval;
        boolean doHarmonic = params.doHarmonic;
        String rmsdFile = params.rmsdFile;
        double timeStep = params.timeStep;
        int interval = 10;
        if (Math.abs(timeStep - 0.001) > 1e-10) {
            int fac = (int) Math.round(0.001 / timeStep);
            thermostatInterval *= fac;
            interval *= fac;
            steps *= fac;
        }

        if (!graphics) {
            System.out.println("Running Iron MC with N="+numAtoms+" at rho="+density+" T="+temperatureK);
            System.out.println(steps+" steps");
            System.out.println("w: "+w);
            if (!doHarmonic && numInnerSteps > 0) System.out.println(numInnerSteps + " inner steps");
            System.out.println("thermostat interval: " + thermostatInterval);
            System.out.println("timeStep: " + timeStep);
        }

        double L = Math.pow(numAtoms/density, 1.0/3.0);
        final SimFe sim = new SimFe(crystal, numAtoms, temperature, density, w, offsetDim, numInnerSteps, swap, doHarmonic, timeStep);
        int[] seeds = sim.getRandomSeeds();
        if (seeds != null) System.out.println("random seeds: " + Arrays.toString(seeds));

        DataSourceEnergies dsEnergies = new DataSourceEnergies(sim.potentialMaster);
        dsEnergies.setPotentialCalculation(new DataSourceEnergies.PotentialCalculationEnergiesEAM(sim.potential));
        dsEnergies.setBox(sim.box);
        IData u = dsEnergies.getData();
        if (!graphics) System.out.println("Fe lattice energy (eV/atom): "+ElectronVolt.UNIT.fromSim(u.getValue(1)/numAtoms));

        MeterStructureFactor meterSfac = new MeterStructureFactor(sim.box, 8);
        if (graphics) {
            sim.integrator.setThermostatInterval(thermostatInterval);
            final String APP_NAME = "SimFe";
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 3);
            ColorScheme colorScheme = new ColorScheme() {
                @Override
                public Color getAtomColor(IAtom a) {
                    int idx0 = a.getLeafIndex();
                    int idx1 = sim.p1ImageHarmonic.getPartner(idx0);
                    return idx0 < idx1 ? Color.RED : Color.BLUE;
                }
            };
            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);
            DiameterHash diameter = simGraphic.getDisplayBox(sim.box).getDiameterHash();
            ((DiameterHashByType)diameter).setDiameter(sim.species.getLeafType(),2);
            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

            DataSourceCountTime tSource = new DataSourceCountTime(sim.integrator);
            AccumulatorHistory energyHist = new AccumulatorHistory(new HistoryCollapsingAverage());
            energyHist.setTimeDataSource(tSource);
            AccumulatorHistory springHist = new AccumulatorHistory(new HistoryCollapsingAverage());
            springHist.setTimeDataSource(tSource);
            DataSplitter splitter = new DataSplitter();
            DataPumpListener energyPump = new DataPumpListener(dsEnergies, splitter, interval);
            sim.integrator.getEventManager().addListener(energyPump);
            splitter.setDataSink(0, springHist);
            splitter.setDataSink(1, energyHist);
            DisplayPlot energyPlot = new DisplayPlot();
            energyPlot.setLabel("Fe");
            energyPlot.setUnit(new CompoundUnit(new Unit[]{new SimpleUnit(Null.DIMENSION,1.0/numAtoms,"why do you want a name.  just use it.","per atom", false)},new double[]{-1}));
//            energyPlot.setUnit(new CompoundUnit(new Unit[]{ElectronVolt.UNIT, new SimpleUnit(Null.DIMENSION,1.0/numAtoms,"why do you want a name.  just use it.","per atom", false)},new double[]{1,-1}));
//            energyPlot.setUnit(ElectronVolt.UNIT);
            energyHist.addDataSink(energyPlot.getDataSet().makeDataSink());
            simGraphic.add(energyPlot);
            energyPlot.setDoLegend(false);
            DisplayPlot springPlot = new DisplayPlot();
            springPlot.setLabel("spring");
            springPlot.setUnit(new CompoundUnit(new Unit[]{new SimpleUnit(Null.DIMENSION,1.0/numAtoms,"why do you want a name.  just use it.","per atom", false)},new double[]{-1}));
            springHist.addDataSink(springPlot.getDataSet().makeDataSink());
            simGraphic.add(springPlot);
            springPlot.setDoLegend(false);

            MeterTemperature meterT = new MeterTemperature(sim.box, 3);
            AccumulatorHistory tHist = new AccumulatorHistory(new HistoryCollapsingAverage());
            tHist.setTimeDataSource(tSource);
            DataPumpListener tPump = new DataPumpListener(meterT, tHist, interval);
            sim.integrator.getEventManager().addListener(tPump);
            DisplayPlot tPlot = new DisplayPlot();
            tPlot.setLabel("T");
            tPlot.setUnit(Kelvin.UNIT);
            simGraphic.add(tPlot);
            tHist.addDataSink(tPlot.getDataSet().makeDataSink());
            tPlot.setDoLegend(false);

            MeterKineticEnergy meterKE = new MeterKineticEnergy(sim.box);
            AccumulatorHistory keHist = new AccumulatorHistory(new HistoryCollapsingAverage());
            keHist.setTimeDataSource(tSource);
            DataPumpListener kePump = new DataPumpListener(meterKE, keHist, interval);
            sim.integrator.getEventManager().addListener(kePump);
            DisplayPlot kePlot = new DisplayPlot();
            kePlot.setLabel("KE");
            kePlot.setUnit(new CompoundUnit(new Unit[]{new SimpleUnit(Null.DIMENSION, 1.0 / numAtoms, "why do you want a name.  just use it.", "per atom", false)}, new double[]{-1}));
            simGraphic.add(kePlot);
            keHist.addDataSink(kePlot.getDataSet().makeDataSink());
            kePlot.setDoLegend(false);

            MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster, sim.box);
            meterPE.setPotentialCalculation(new DataSourceEnergies.PotentialCalculationEnergiesEAM(sim.potential));
            MeterEnergy meterE = new MeterEnergy(sim.potentialMaster, sim.box);
            meterE.setPotential(meterPE);
            AccumulatorHistory eHist = new AccumulatorHistory(new HistoryScrolling(1000));
            eHist.setTimeDataSource(tSource);
            DataPumpListener ePump = new DataPumpListener(meterE, eHist, interval);
            sim.integrator.getEventManager().addListener(ePump);
            DisplayPlot tePlot = new DisplayPlot();
            tePlot.setLabel("total E");
            tePlot.setUnit(new CompoundUnit(new Unit[]{new SimpleUnit(Null.DIMENSION, 1.0 / numAtoms, "why do you want a name.  just use it.", "per atom", false)}, new double[]{-1}));
            simGraphic.add(tePlot);
            eHist.addDataSink(tePlot.getDataSet().makeDataSink());
            tePlot.setDoLegend(false);

            AccumulatorHistory dudwHist = new AccumulatorHistory(new HistoryCollapsingAverage());
            splitter.setDataSink(2, dudwHist);
            dudwHist.setTimeDataSource(tSource);
            DisplayPlot dudwPlot = new DisplayPlot();
            dudwHist.addDataSink(dudwPlot.getDataSet().makeDataSink());
            dudwPlot.setLabel("dudw");
            dudwPlot.setUnit(new CompoundUnit(new Unit[]{new SimpleUnit(Null.DIMENSION,1.0/numAtoms,"why do you want a name.  just use it.","per atom", false)},new double[]{-1}));
            simGraphic.add(dudwPlot);

            AccumulatorAverageFixed avgSfac = new AccumulatorAverageFixed(1);
            avgSfac.setPushInterval(1);
            DataPumpListener pumpSfac = new DataPumpListener(meterSfac, avgSfac, 100 * interval);
            sim.integrator.getEventManager().addListener(pumpSfac);
            DisplayPlot plotSfac = new DisplayPlot();
            avgSfac.addDataSink(plotSfac.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{avgSfac.AVERAGE});
            plotSfac.setLabel("Structure Factor");
            plotSfac.setDoDrawLines(new DataTag[]{meterSfac.getTag()},false);
            plotSfac.getPlot().setYLog(true);
            simGraphic.add(plotSfac);

            if (sim.integrator instanceof IntegratorMDHarmonicMC) {
                DataSourceScalar meterAcc = new DataSourceScalar("acc", Null.DIMENSION) {
                    @Override
                    public double getDataAsScalar() {
                        return ((IntegratorMDHarmonicMC) sim.integrator).getAcceptanceProbability();
                    }
                };
                AccumulatorHistory aHist = new AccumulatorHistory(new HistoryCollapsingDiscard());
                aHist.setTimeDataSource(tSource);
                DataPumpListener aPump = new DataPumpListener(meterAcc, aHist, interval);
                sim.integrator.getEventManager().addListener(aPump);
                DisplayPlot aPlot = new DisplayPlot();
                aPlot.setLabel("Acc");
                simGraphic.add(aPlot);
                aHist.addDataSink(aPlot.getDataSet().makeDataSink());
                aPlot.setDoLegend(false);
            }

            MeterRMSD meterRMSD = new MeterRMSD(sim.box, sim.space);
            AccumulatorHistory rmsdHist = new AccumulatorHistory(new HistoryCollapsingAverage());
            rmsdHist.setTimeDataSource(tSource);
            DataPumpListener rmsdPump = new DataPumpListener(meterRMSD, rmsdHist, interval);
            sim.integrator.getEventManager().addListener(rmsdPump);
            DisplayPlot rmsdPlot = new DisplayPlot();
            rmsdPlot.setLabel("RMSD");
            simGraphic.add(rmsdPlot);
            rmsdHist.addDataSink(rmsdPlot.getDataSet().makeDataSink());
            rmsdPlot.setDoLegend(false);

            simGraphic.makeAndDisplayFrame(APP_NAME);

            return;
        }

        MeterRMSD meterRMSD = new MeterRMSD(sim.box, sim.space);

        sim.integrator.setThermostatInterval(10);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps / 10));

        sim.integrator.resetStepCount();
        sim.integrator.setThermostatInterval(thermostatInterval);
if (nve) {
            sim.integrator.setIsothermal(false);
        }
        if (sim.integrator instanceof IntegratorMDHarmonicMC) {
            ((IntegratorMDHarmonicMC) sim.integrator).resetAcceptance();
        }

        System.out.println("equilibration finished (" + steps / 10 + " steps)");

        long t1 = System.currentTimeMillis();

        if (thermostatInterval > interval * 10) interval *= 2;
        int numBlocks = 100;
        if (steps / numBlocks < thermostatInterval) numBlocks = (int) (steps / thermostatInterval);
        long blockSize = steps / numBlocks / interval;
        if (blockSize==0) blockSize = 1;
        System.out.println("numBlocks: " + numBlocks);

        AccumulatorAverageCovariance accEnergies = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener pumpEnergies = new DataPumpListener(dsEnergies, accEnergies, interval);
        sim.integrator.getEventManager().addListener(pumpEnergies);

        DataLogger dataLogger = null;
        if (rmsdFile != null) {
            dataLogger = new DataLogger();
            DataPumpListener pumpRMSD = new DataPumpListener(meterRMSD, dataLogger, 5 * interval);
            sim.integrator.getEventManager().addListener(pumpRMSD);
            dataLogger.setFileName(rmsdFile);
            DataArrayWriter writer = new DataArrayWriter();
            writer.setIncludeHeader(false);
            dataLogger.setDataSink(writer);
            dataLogger.setAppending(false);
        }
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));

        if (dataLogger != null) {
            dataLogger.cleanUp();
        }


        if (sim.mcMoveSwap != null) {
            System.out.println("swap acceptance: " + sim.mcMoveSwap.getTracker().acceptanceProbability());
        }
        if (sim.integrator instanceof IntegratorMDHarmonicMC) {
            System.out.println("MC accepted fraction " + ((IntegratorMDHarmonicMC) sim.integrator).getAcceptanceProbability());
        }

        IData avgEnergies = accEnergies.getData(accEnergies.AVERAGE);
        IData errEnergies = accEnergies.getData(accEnergies.ERROR);
        IData corEnergies = accEnergies.getData(accEnergies.BLOCK_CORRELATION);
        IData covEnergies = accEnergies.getData(accEnergies.BLOCK_COVARIANCE);

        System.out.println("spring energy: "+avgEnergies.getValue(0)/numAtoms+"   error: "+errEnergies.getValue(0)/numAtoms+"  cor: "+corEnergies.getValue(0));
        System.out.println("Fe energy: "+avgEnergies.getValue(1)/numAtoms+"   error: "+errEnergies.getValue(1)/numAtoms+"  cor: "+corEnergies.getValue(1));
        System.out.println("du/dw: "+avgEnergies.getValue(2)/numAtoms+"   error: "+errEnergies.getValue(2)/numAtoms+"  cor: "+corEnergies.getValue(2));
        double var0 = covEnergies.getValue(0*3+0);
        double var1 = covEnergies.getValue(1*3+1);
        double var2 = covEnergies.getValue(2*3+2);
        double cor01 = 0;
        if (var0*var1>0) cor01 = covEnergies.getValue(0*3+1)/Math.sqrt(covEnergies.getValue(0*3+0)*covEnergies.getValue(1*3+1));
        double cor02 = 0;
        if (var0*var2>0) cor02 = covEnergies.getValue(0*3+2)/Math.sqrt(covEnergies.getValue(0*3+0)*covEnergies.getValue(2*3+2));
        System.out.println("spring correlation: 1 "+cor01+" "+cor02);
        double cor12 = covEnergies.getValue(1*3+2)/Math.sqrt(covEnergies.getValue(1*3+1)*covEnergies.getValue(2*3+2));
        System.out.println("Fe correlation: "+cor01+" 1 "+cor12);
        System.out.println("du/dw correlation: "+cor02+" "+cor12+" 1");

        long t2 = System.currentTimeMillis();
        System.out.println("time: "+(t2-t1)/1000.0+" seconds");
    }

    enum Crystal {FCC, BCC, HCP}

    public static class LjMC3DParams extends ParameterBase {
        public int numAtoms = 500;
        public double T = 2.0;
        public double density = 0.3;
        public long steps = 100000;
        public boolean graphics = false;
        public double w = 1;
        public int offsetDim = 0;
        public int numInnerSteps = 10;
        public Crystal crystal = Crystal.FCC;
        public boolean swap = true;
        public boolean nve = false;
        public int thermostatInterval = 10;
        public boolean doHarmonic = true;
        public String rmsdFile = null;
        public double timeStep = 0.001;
    }

}
