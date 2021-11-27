/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.reproducibility;


import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayPlotXChart;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveMoleculeRotate;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.LatticeCubicFcc;
import etomica.models.water.P2WaterSPCE;
import etomica.models.water.SpeciesWater3P;
import etomica.potential.*;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeEwaldFourier;
import etomica.potential.ewald.P2Ewald1Real;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.units.*;
import etomica.units.dimensions.Null;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;

public class WaterNPTMC extends Simulation {

    public PotentialCompute potentialMaster;
    public IntegratorMC integrator;
    public ISpecies species;
    public MCMoveMolecule translateMove;
    public MCMoveMoleculeRotate rotateMove;
    public MCMoveVolume volumeMove;

    public WaterNPTMC(int numMolecules, double temperature, double density, double pressure) {
        super(Space3D.getInstance());
        Box box = new Box(space);
        species = SpeciesWater3P.create(true);
        addSpecies(species);
        addBox(box);
        box.setNMolecules(species, numMolecules);
        box.setDensity(density);
        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);


        AtomType oType = species.getTypeByName("O");
        AtomType hType = species.getTypeByName("H");

        double sigmaLJ = P2WaterSPCE.SIGMA;
        double epsilonLJ = P2WaterSPCE.EPSILON;
        double chargeO = P2WaterSPCE.QO;
        double chargeH = P2WaterSPCE.QH;

        PotentialComputeEwaldFourier ewaldFourier = new PotentialComputeEwaldFourier(getSpeciesManager(), box);
        PotentialComputeEwaldFourier.EwaldParams params = ewaldFourier.getOptimalParams(3, 0);

        ewaldFourier.setkCut(params.kCut);
        ewaldFourier.setCharge(oType, chargeO);
        ewaldFourier.setCharge(hType, chargeH);
        ewaldFourier.setAlpha(params.alpha);

        PotentialMaster pm = new PotentialMaster(getSpeciesManager(), box, BondingInfo.noBonding());

        TruncationFactory tf = new TruncationFactorySimple(params.rCut);
        P2LennardJones p2OOLJ = new P2LennardJones(sigmaLJ, epsilonLJ);
        P2Ewald1Real p2OOqq = new P2Ewald1Real(chargeO*chargeO, params.alpha);
        Potential2Soft p2OO = tf.make(p2OOLJ, p2OOqq);
        Potential2Soft p2HH = tf.make(new P2Ewald1Real(chargeH*chargeH, params.alpha));
        P2HardGeneric p2MHC = P2HardSphere.makePotential(0.1);
        Potential2Soft p2OH = tf.make(p2MHC, new P2Ewald1Real(chargeH*chargeO, params.alpha));
        pm.setPairPotential(oType, oType, p2OO);
        pm.setPairPotential(hType, hType, p2HH);
        pm.setPairPotential(oType, hType, p2OH);

        PotentialMasterBonding pmBonding = ewaldFourier.makeIntramolecularCorrection();

        potentialMaster = new PotentialComputeAggregate(pm, ewaldFourier, pmBonding);

        integrator = new IntegratorMC(potentialMaster, random, temperature, box);

        translateMove = new MCMoveMolecule(random, potentialMaster, box);
        integrator.getMoveManager().addMCMove(translateMove);

        rotateMove = new MCMoveMoleculeRotate(random, potentialMaster, box);
        integrator.getMoveManager().addMCMove(rotateMove);

        if (pressure >= 0) {
            volumeMove = new MCMoveVolume(integrator, random, pressure);
        }

    }

    public static void main(String[] args) {

        WaterNPTMC.SimParams params = new WaterNPTMC.SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            // modify parameters here for interactive testing
            params.steps = 120 * 1000 * params.numMolecules;
        }

        double temperatureK = params.temperatureK;
        double temperature = Kelvin.UNIT.toSim(temperatureK);
        Unit pUnit = new PrefixedUnit(Prefix.KILO, Pascal.UNIT);
        double pressureKPa = params.pressureKPa;
        double pressure = pUnit.toSim(pressureKPa);

        Unit dUnit = new SimpleUnit(Null.DIMENSION, 1/((Oxygen.INSTANCE.getMass()+2* Hydrogen.INSTANCE.getMass()) /Constants.AVOGADRO*1e24), "Density", "g/cm^3", false);
        double density = dUnit.toSim(params.density);

        WaterNPTMC sim = new WaterNPTMC(params.numMolecules, temperature, density, pressure);

        if (true) {
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator, 5000000));
            sim.getController().addActionSequential(new IAction() {
                @Override
                public void actionPerformed() {
                    sim.integrator.getMoveManager().addMCMove(sim.volumeMove);
                }
            });
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
            SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "Water NPT");
            ((ColorSchemeByType) graphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.species.getTypeByName("H"), Color.WHITE);
            ((ColorSchemeByType) graphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.species.getTypeByName("O"), Color.RED);

            DataSourceCountSteps timeSource = new DataSourceCountSteps(sim.integrator);

            MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
            AccumulatorHistory accPE = new AccumulatorHistory(new HistoryCollapsingAverage());
            accPE.setTimeDataSource(timeSource);
            DataPumpListener pumpPE = new DataPumpListener(meterPE, accPE);
            sim.integrator.getEventManager().addListener(pumpPE);
            DisplayPlotXChart plotPE = new DisplayPlotXChart();
            plotPE.setLabel("energy");
            accPE.setDataSink(plotPE.getDataSet().makeDataSink());
            graphic.add(plotPE);

            MeterDensity meterDensity = new MeterDensity(sim.box());
            AccumulatorHistory accDensity = new AccumulatorHistory(new HistoryCollapsingAverage());
            accDensity.setTimeDataSource(timeSource);
            DataPumpListener pumpDensity = new DataPumpListener(meterDensity, accDensity, 10);
            sim.integrator.getEventManager().addListener(pumpDensity);

            DisplayPlot historyDensity = new DisplayPlot();
            accDensity.setDataSink(historyDensity.getDataSet().makeDataSink());
            historyDensity.setLabel("Density");
            historyDensity.setUnit(dUnit);
            graphic.add(historyDensity);

            MeterPressure meterPressure = new MeterPressure(sim.box(), sim.potentialMaster);
            meterPressure.setTemperature(temperature);
            AccumulatorHistory accPressure = new AccumulatorHistory(new HistoryCollapsingAverage());
            accPressure.setTimeDataSource(timeSource);
            DataPumpListener pumpPressure = new DataPumpListener(meterPressure, accPressure, 2* params.numMolecules);
            sim.integrator.getEventManager().addListener(pumpPressure);

            DisplayPlot historyPressure = new DisplayPlot();
            accPressure.setDataSink(historyPressure.getDataSet().makeDataSink());
            historyPressure.setLabel("Pressure");
            historyPressure.setUnit(pUnit);
            graphic.add(historyPressure);

            graphic.makeAndDisplayFrame();

            return;
        }
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator));
    }

    public static class SimParams extends ParameterBase {
        public long steps = 10000000;
        public double density = 1;
        public double temperatureK = 300;
        public int numMolecules = 1100;
        public double pressureKPa = 101;
        public double s = 3;
    }

}
