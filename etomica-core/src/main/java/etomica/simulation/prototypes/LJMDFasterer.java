/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountTimeFasterer;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.meter.MeterPotentialEnergyFromIntegratorFasterer;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.DisplayPlotXChart;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorListenerNHC;
import etomica.integrator.IntegratorVelocityVerletFasterer;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.list.PotentialMasterListFasterer;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.SimpleUnit;
import etomica.units.dimensions.Null;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
public class LJMDFasterer extends Simulation {

    public IntegratorVelocityVerletFasterer integrator;
    public SpeciesGeneral species;
    public Box box;
    public IntegratorListenerNHC nhc;


    public LJMDFasterer(Space space, double density, int numAtoms, boolean useNbrLists, boolean doNHC) {
        super(space);

        species = SpeciesGeneral.monatomic(space, AtomType.simple("A", 20), true);
        addSpecies(species);

        box = this.makeBox();
        box.setNMolecules(species, numAtoms);
        new BoxInflate(box, space, density).actionPerformed();
        double rc = 3;
        PotentialMasterFasterer potentialMaster = useNbrLists ? new PotentialMasterListFasterer(getSpeciesManager(), box, 2, rc + 1, BondingInfo.noBonding())
                : new PotentialMasterFasterer(getSpeciesManager(), box, BondingInfo.noBonding());
        integrator = new IntegratorVelocityVerletFasterer(potentialMaster, random, 0.02, 1.0, box);
        if (doNHC) {
            integrator.setIsothermal(false);
            nhc = new IntegratorListenerNHC(integrator, random, 3, 2);
            integrator.getEventManager().addListener(nhc);
        }
        getController().setSleepPeriod(1);
        getController().addActivity(new ActivityIntegrate(integrator));

        double sigma = 1.0;
        Potential2Soft potential = P2LennardJones.makeTruncated(space, sigma, 1.0, new TruncationFactoryForceShift(space, rc));

        AtomType leafType = species.getLeafType();

        potentialMaster.setPairPotential(leafType, leafType, potential);
        if (!useNbrLists) {
            BoxImposePbc imposepbc = new BoxImposePbc(space);
            imposepbc.setBox(box);
            integrator.getEventManager().addListener(new IntegratorListenerAction(imposepbc));
        }

        ConfigurationLattice configuration = new ConfigurationLattice(space.D() == 3 ? new LatticeCubicFcc(space) : new LatticeOrthorhombicHexagonal(space), space);
        configuration.initializeCoordinates(box);
    }

    public static void main(String[] args) {
        final String APP_NAME = "LJMD3D";
        final LJMDFasterer sim = new LJMDFasterer(Space3D.getInstance(), 0.8, 864, true, true);
        final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 3);

        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

        MeterPotentialEnergyFromIntegratorFasterer energy = new MeterPotentialEnergyFromIntegratorFasterer(sim.integrator);
        AccumulatorAverageCollapsing avgEnergy = new AccumulatorAverageCollapsing();
        avgEnergy.setPushInterval(10);
        DataPumpListener pump = new DataPumpListener(energy, avgEnergy, 10);
        sim.integrator.getEventManager().addListener(pump);

        simGraphic.getController().getDataStreamPumps().add(pump);

        DisplayTextBoxesCAE display = new DisplayTextBoxesCAE();
        display.setAccumulator(avgEnergy);
        display.setUnit(new SimpleUnit(Null.DIMENSION, 864, "1/N", "1/N", false));
        simGraphic.add(display);
        DataSourceCountTimeFasterer timer = new DataSourceCountTimeFasterer(sim.integrator);

        if (sim.nhc != null) {
            IntegratorListenerNHC.DataSourceTotalEnergy meterTotalEnergy = new IntegratorListenerNHC.DataSourceTotalEnergy(sim.integrator, sim.nhc);
            AccumulatorHistory historyTotalEnergy = new AccumulatorHistory(new HistoryCollapsingDiscard());
            historyTotalEnergy.setTimeDataSource(timer);
            DataPumpListener pumpTotalEnergy = new DataPumpListener(meterTotalEnergy, historyTotalEnergy);
            sim.integrator.getEventManager().addListener(pumpTotalEnergy);
            DisplayPlotXChart plotTotalEnergy = new DisplayPlotXChart();
            plotTotalEnergy.setLabel("energy");
            historyTotalEnergy.addDataSink(plotTotalEnergy.makeSink("total"));
            simGraphic.add(plotTotalEnergy);
        }

        MeterTemperature meterT = new MeterTemperature(sim.box, 3);
        AccumulatorHistory historyT = new AccumulatorHistory(new HistoryCollapsingDiscard());
        historyT.setTimeDataSource(timer);
        DataPumpListener pumpT = new DataPumpListener(meterT, historyT);
        sim.integrator.getEventManager().addListener(pumpT);
        DisplayPlotXChart plotT = new DisplayPlotXChart();
        plotT.setLabel("T");
        historyT.addDataSink(plotT.makeSink("T"));
        simGraphic.add(plotT);

        simGraphic.makeAndDisplayFrame(APP_NAME);
    }
}
