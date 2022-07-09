/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedDensity.positionOrientation;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.config.ConformationChainLinear;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.dcvgcmd.P1WCAWall;
import etomica.graphics.ColorSchemeRandomByMolecule;
import etomica.graphics.DisplayPlotXChart;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveMoleculeRotate;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.BondingInfo;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedForceShifted;
import etomica.potential.PotentialMaster;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeField;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.units.Degree;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.IRandom;

import java.util.Arrays;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
public class LJMCDimer extends Simulation {

    public IntegratorMC integrator;
    public SpeciesGeneral species;
    public Box box;
    public P2LennardJones potential;
    public final PotentialMaster potentialMaster;
    public final double Lz;

    public LJMCDimer(int D, int numMolecules, double density, double Lz, double temperature) {
        super(Space.getInstance(D));

        this.Lz = Lz;
        double a = Degree.UNIT.toSim(45);
        double[] a2 = new double[D - 1];
        Arrays.fill(a2, a);
        species = new SpeciesBuilder(space)
                .addCount(AtomType.simpleFromSim(this), 2)
                .setDynamic(true)
                .withConformation(new ConformationChainLinear(space, 1, a2))
                .build();
        addSpecies(species);

        double sigma = 1.0;

        box = this.makeBox(new BoundaryRectangularSlit(2, getSpace()));
        potentialMaster = new PotentialMasterCell(this.getSpeciesManager(), box, 2, BondingInfo.noBonding());
        P1WCAWall wall = new P1WCAWall(box, 1, 1);
        PotentialComputeField pcField = new PotentialComputeField(getSpeciesManager(), box);
        AtomType leafType = species.getLeafType();
        pcField.setFieldPotential(leafType, wall);
        PotentialComputeAggregate pcAgg = new PotentialComputeAggregate(potentialMaster, pcField);

        integrator = new IntegratorMC(pcAgg, this.getRandom(), 1.0, box);
        integrator.setTemperature(temperature);

        integrator.getMoveManager().addMCMove(new MyMCMoveMolecule(this.getRandom(), pcAgg, box));
        integrator.getMoveManager().addMCMove(new MCMoveMoleculeRotate(this.getRandom(), pcAgg, box));

        box.setNMolecules(species, numMolecules);
        if (D == 2) {
            double Lxy = numMolecules / density / Lz;
            box.getBoundary().setBoxSize(Vector.of(Lxy, Lz - 1));
        } else {
            double Lxy = Math.sqrt(numMolecules / density / Lz);
            box.getBoundary().setBoxSize(Vector.of(Lxy, Lxy, Lz - 1));
        }

        potential = new P2LennardJones(sigma, 1.0);
        P2SoftSphericalTruncatedForceShifted p2 = new P2SoftSphericalTruncatedForceShifted(potential, 3.0);
        potentialMaster.setPairPotential(leafType, leafType, p2);

        ConfigurationLattice configuration = new ConfigurationLattice(D == 3 ? new LatticeCubicFcc(space) : new LatticeOrthorhombicHexagonal(space), space);
        configuration.initializeCoordinates(box);
        Vector l = box.getBoundary().getBoxSize();
        l.setX(2, Lz);
        box.getBoundary().setBoxSize(l);
    }

    public static void main(String[] args) {
        DimerParams params = new DimerParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            // override defaults here
            params.D = 2;
            params.Lz = 4;
            params.density = 0.01;
            params.numMolecules = 10;
            params.temperature = 1.0;
            params.nZdata = 20;
            params.nAngleData = 40;
//            params.D = 2;
//            params.Lz = 20;
//            params.density = 0.35;
//            params.numMolecules = 100;
//            params.temperature = 1.0;
//            params.nZdata = 10;
//            params.nAngleData = 20;

        }
        if(params.D == 2) {
            System.out.println("Estimated accessible density: "+params.density*params.Lz/(params.Lz-2)+", "+params.density*params.Lz/(params.Lz-2)/(2*Math.PI));
        }

        int numMolecules = params.numMolecules;
        double temperature = params.temperature;
        double density = params.density;
        int D = params.D;
        double Lz = params.Lz;
        int nZdata = params.nZdata;
        int nAngleData = params.nAngleData;

        final LJMCDimer sim = new LJMCDimer(D, numMolecules, density, Lz, temperature);

        boolean graphics = true;

        MeterPotentialEnergyFromIntegrator energy = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        MeterHistogramOrientation meterOrientation = new MeterHistogramOrientation(sim.box, 2, nZdata, nAngleData);
        MeterHistogramOrientationMA meterOrientationMA = new MeterHistogramOrientationMA(sim.potentialMaster, sim.box, params.temperature, 2, nZdata, nAngleData);
        MeterDensityProfile meterDensity = new MeterDensityProfile(sim.box, 2, 20 * nZdata);
        MeterDensityProfileMA meterDensityMA = new MeterDensityProfileMA(sim.potentialMaster, sim.box, params.temperature, 2, 20 * nZdata);
        MeterDensityProfileForceSimple meterDensityForceSimple = new MeterDensityProfileForceSimple(sim.potentialMaster, sim.box, params.temperature, 2,  20*nZdata);

        if (graphics) {
            final String APP_NAME = "LJMDDimer";
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 3);

            AccumulatorAverageCollapsing avgEnergy = new AccumulatorAverageCollapsing();
            avgEnergy.setPushInterval(10);
            DataPumpListener pump = new DataPumpListener(energy, avgEnergy, numMolecules);
            sim.integrator.getEventManager().addListener(pump);

            // Orientation-density mapped averaging
            AccumulatorAverageFixed avgHistogramMA = new AccumulatorAverageFixed(1000);
            DataPumpListener pumpHistogramMA = new DataPumpListener(meterOrientationMA, avgHistogramMA, numMolecules);
            sim.integrator.getEventManager().addListener(pumpHistogramMA);
            DataHistogramSplitter splitterMA = new DataHistogramSplitter();
            avgHistogramMA.addDataSink(splitterMA, new AccumulatorAverageFixed.StatType[]{avgHistogramMA.AVERAGE});
            DisplayPlotXChart plotHistogramMA = new DisplayPlotXChart();
            DataDoubleArray zDataMA = meterOrientationMA.getIndependentData(0);
            for (int i = nZdata - 1; i >= 0; i--) {
                splitterMA.setDataSink(i, plotHistogramMA.getDataSet().makeDataSink());
                plotHistogramMA.setLegend(new DataTag[]{splitterMA.getTag(i)}, "z=" + (0.5 * Lz - zDataMA.getValue(i)));
            }
            plotHistogramMA.setLabel("HistogramsMA");
            simGraphic.add(plotHistogramMA);

            // Orientation-density histogramming
            AccumulatorAverageFixed avgHistogram = new AccumulatorAverageFixed(1000);
            DataPumpListener pumpHistogram = new DataPumpListener(meterOrientation, avgHistogram, numMolecules);
            sim.integrator.getEventManager().addListener(pumpHistogram);
            DataHistogramSplitter splitter = new DataHistogramSplitter();
            avgHistogram.addDataSink(splitter, new AccumulatorAverageFixed.StatType[]{avgHistogram.AVERAGE});
            DisplayPlotXChart plotHistogram = new DisplayPlotXChart();
            DataDoubleArray zData = meterOrientation.getIndependentData(0);
            for (int i = nZdata - 1; i >= 0; i--) {
                splitter.setDataSink(i, plotHistogram.getDataSet().makeDataSink());
                plotHistogram.setLegend(new DataTag[]{splitter.getTag(i)}, "z=" + (0.5 * Lz - zData.getValue(i)));
            }
            plotHistogram.setLabel("Histograms");
            simGraphic.add(plotHistogram);

            // Density plot -- conventional
            AccumulatorAverageFixed avgDensityProfile = new AccumulatorAverageFixed(1000);
            DataPumpListener pumpDensityProfile = new DataPumpListener(meterDensity, avgDensityProfile, numMolecules);
            sim.integrator.getEventManager().addListener(pumpDensityProfile);
            DisplayPlotXChart plotDensityProfile = new DisplayPlotXChart();
            avgDensityProfile.addDataSink(plotDensityProfile.getDataSet().makeDataSink(), new AccumulatorAverageFixed.StatType[]{avgHistogram.AVERAGE});
            plotDensityProfile.setLabel("Density profile");
            simGraphic.add(plotDensityProfile);

            // Density plot -- mapped average
            AccumulatorAverageFixed avgDensityProfileMA = new AccumulatorAverageFixed(1000);
            DataPumpListener pumpDensityProfileMA = new DataPumpListener(meterDensityMA, avgDensityProfileMA, numMolecules);
            sim.integrator.getEventManager().addListener(pumpDensityProfileMA);
            DisplayPlotXChart plotDensityProfileMA = new DisplayPlotXChart();
            avgDensityProfileMA.addDataSink(plotDensityProfileMA.getDataSet().makeDataSink(), new AccumulatorAverageFixed.StatType[]{avgHistogram.AVERAGE});
            plotDensityProfileMA.setLabel("Density profile MA");
            simGraphic.add(plotDensityProfileMA);

            // Density plot -- simple force
            AccumulatorAverageFixed avgDensityProfileF = new AccumulatorAverageFixed(1000);
            DataPumpListener pumpDensityProfileF = new DataPumpListener(meterDensityForceSimple, avgDensityProfileF, numMolecules);
            sim.integrator.getEventManager().addListener(pumpDensityProfileF);
            DisplayPlotXChart plotDensityProfileF = new DisplayPlotXChart();
            avgDensityProfileF.addDataSink(plotDensityProfileF.getDataSet().makeDataSink(), new AccumulatorAverageFixed.StatType[]{avgHistogram.AVERAGE});
            plotDensityProfileF.setLabel("Density profile Force");
            simGraphic.add(plotDensityProfileF);

            DataProcessor dpDifficulty = new DataProcessor() {
                DataFunction data;

                @Override
                protected IData processData(IData inputData) {
                    data.E(inputData);
                    data.TE(Math.sqrt(sim.integrator.getStepCount()));
                    return data;
                }

                @Override
                protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
                    data = (DataFunction) inputDataInfo.makeData();
                    return inputDataInfo;
                }
            };
            avgHistogram.addDataSink(dpDifficulty, new AccumulatorAverageFixed.StatType[]{avgHistogram.ERROR});
            DataHistogramSplitter splitterErr = new DataHistogramSplitter();
            dpDifficulty.setDataSink(splitterErr);
            DisplayPlotXChart plotHistogramErr = new DisplayPlotXChart();
            for (int i = nZdata - 1; i >= 0; i--) {
                splitterErr.setDataSink(i, plotHistogramErr.getDataSet().makeDataSink());
                plotHistogramErr.setLegend(new DataTag[]{splitterErr.getTag(i)}, "z=" + (0.5 * Lz - zData.getValue(i)));
            }
            plotHistogramErr.setLabel("Difficulty");
            simGraphic.add(plotHistogramErr);


            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
            simGraphic.getController().getDataStreamPumps().add(pump);
            simGraphic.getController().getDataStreamPumps().add(pumpHistogram);
            simGraphic.getDisplayBox(sim.box).setColorScheme(new ColorSchemeRandomByMolecule(sim.getSpeciesManager(), sim.box, sim.getRandom()));

            simGraphic.makeAndDisplayFrame(APP_NAME);

            DisplayTextBoxesCAE display = new DisplayTextBoxesCAE();
            display.setAccumulator(avgEnergy);
            simGraphic.add(display);
            return;
        }
    }

    public static class DimerParams extends ParameterBase {
        public int numMolecules = 256;
        public double temperature = 2;
        public double density = 0.4;
        public double Lz = 10;
        public int D = 3;
        public int nZdata = 10;
        public int nAngleData = 10;
    }

    // checks that molecule is in box before attempting to compute energy
    class MyMCMoveMolecule extends MCMoveMolecule {
        public MyMCMoveMolecule(IRandom random, PotentialCompute potentialCompute, Box box) {
            super(random, potentialCompute, box);
        }

        public double getChi(double temperature) {
            int D = box.getSpace().D();
            double r1z = Lz/2 + molecule.getChildList().get(0).getPosition().getX(D);
            double r2z = Lz/2 + molecule.getChildList().get(1).getPosition().getX(D);
            if (r1z < 0 || r1z > Lz || r2z < 0 || r2z > Lz) return 0.0;
            return super.getChi(temperature);
        }
    }
}
