/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedDensity.positionOrientation;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.config.ConformationChainLinear;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.DataTag;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDoubleArray;
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

import java.util.Arrays;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
public class LJMCDimer extends Simulation {

    public IntegratorMC integrator;
    public SpeciesGeneral species;
    public Box box;
    public P2LennardJones potential;

    public LJMCDimer(int D, int numMolecules, double density, double Lz, double temperature) {
        super(Space.getInstance(D));

        double a = Degree.UNIT.toSim(45);
        double[] a2 = new double[D-1];
        Arrays.fill(a2, a);
        species = new SpeciesBuilder(space)
                .addCount(AtomType.simpleFromSim(this), 2)
                .setDynamic(true)
                .withConformation(new ConformationChainLinear(space, 1, a2))
                .build();
        addSpecies(species);

        double sigma = 1.0;
        box = this.makeBox(new BoundaryRectangularSlit(2, getSpace()));
        PotentialMaster potentialMaster = new PotentialMasterCell(this.getSpeciesManager(), box, 2, BondingInfo.noBonding());
        P1WCAWall wall = new P1WCAWall(box, 1, 1);
        PotentialComputeField pcField = new PotentialComputeField(getSpeciesManager(), box);
        AtomType leafType = species.getLeafType();
        pcField.setFieldPotential(leafType, wall);
        PotentialComputeAggregate pcAgg = new PotentialComputeAggregate(potentialMaster, pcField);

        integrator = new IntegratorMC(pcAgg, this.getRandom(), 1.0, box);
        integrator.setTemperature(temperature);

        integrator.getMoveManager().addMCMove(new MCMoveMolecule(this.getRandom(), pcAgg, box));
        integrator.getMoveManager().addMCMove(new MCMoveMoleculeRotate(this.getRandom(), pcAgg, box));

        box.setNMolecules(species, numMolecules);
        if (D == 2) {
            double Lxy = numMolecules / density / Lz;
            box.getBoundary().setBoxSize(Vector.of(Lxy, Lz - 1));
        }
        else {
            double Lxy = Math.sqrt(numMolecules / density / Lz);
            box.getBoundary().setBoxSize(Vector.of(Lxy, Lxy, Lz -1));
        }

        potential = new P2LennardJones(sigma, 1.0);
        P2SoftSphericalTruncatedForceShifted p2 = new P2SoftSphericalTruncatedForceShifted(potential, 3.0);
        potentialMaster.setPairPotential(leafType, leafType, p2);

        ConfigurationLattice configuration = new ConfigurationLattice(D==3 ? new LatticeCubicFcc(space) : new LatticeOrthorhombicHexagonal(space), space);
        configuration.initializeCoordinates(box);
        Vector l = box.getBoundary().getBoxSize();
        l.setX(2, Lz);
        box.getBoundary().setBoxSize(l);
    }

    public static void main(String[] args) {
        DimerParams params = new DimerParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            // override defaults here
            params.D = 2;
            params.density = 0.4;
            params.numMolecules = 100;
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


        if (graphics) {
            final String APP_NAME = "LJMDDimer";
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 3);

            AccumulatorAverageCollapsing avgEnergy = new AccumulatorAverageCollapsing();
            avgEnergy.setPushInterval(10);
            DataPumpListener pump = new DataPumpListener(energy, avgEnergy, numMolecules);
            sim.integrator.getEventManager().addListener(pump);

            AccumulatorAverageFixed avgHistogram = new AccumulatorAverageFixed(1);
            avgEnergy.setPushInterval(10);
            DataPumpListener pumpHistogram = new DataPumpListener(meterOrientation, avgHistogram, numMolecules);
            sim.integrator.getEventManager().addListener(pumpHistogram);
            DataHistogramSplitter splitter = new DataHistogramSplitter();
            avgHistogram.addDataSink(splitter, new AccumulatorAverageFixed.StatType[]{avgHistogram.AVERAGE});
            DisplayPlotXChart plotHistogram = new DisplayPlotXChart();
            DataDoubleArray zData = meterOrientation.getIndependentData(0);
            for (int i=nZdata-1; i>=0; i--) {
                splitter.setDataSink(i, plotHistogram.getDataSet().makeDataSink());
                plotHistogram.setLegend(new DataTag[]{splitter.getTag(i)}, "z="+(0.5*Lz-zData.getValue(i)));
            }
            plotHistogram.setLabel("Plot");
            simGraphic.add(plotHistogram);

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
}
