/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.rotation;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPumpListener;
import etomica.data.history.HistoryScrolling;
import etomica.data.meter.MeterEnergyFromIntegrator;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayPlotXChart;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.ShakeListener;
import etomica.models.water.ConformationWater3P;
import etomica.models.water.P2WaterSPC;
import etomica.models.water.SpeciesWater3P;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesAgentManager;
import etomica.species.SpeciesGeneral;
import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.util.Constants;

import java.awt.*;

public class WaterDropletShake {

    public static SimulationGraphic makeWaterDroplet() {
        Space space = Space3D.getInstance();
        Simulation sim = new Simulation(space);
        Box box = new Box(new BoundaryRectangularNonperiodic(sim.getSpace()), space);
        SpeciesGeneral species = SpeciesWater3P.create(true);
        sim.addSpecies(species);
        sim.addBox(box);
        box.setNMolecules(species, 108);
        box.setDensity(1 / 18.0 * Constants.AVOGADRO / 1E24);
        ConfigurationWater108 configFile = new ConfigurationWater108();
        configFile.initializeCoordinates(box);
        PotentialMaster potentialMaster = new PotentialMaster(sim.getSpeciesManager(), box, BondingInfo.noBonding());
        double timeInterval = 0.002;
        int maxIterations = 100;
        final IntegratorVelocityVerlet integrator = new IntegratorVelocityVerlet(potentialMaster, sim.getRandom(), timeInterval, Kelvin.UNIT.toSim(271.2654804973), box);
        double lOH = ConformationWater3P.bondLengthOH;
        double lHH = Math.sqrt(2 * lOH * lOH * (1 - Math.cos(ConformationWater3P.angleHOH)));
        SpeciesAgentManager<ShakeListener.BondConstraints> shakeAgents = new SpeciesAgentManager<>(new SpeciesAgentManager.AgentSource<ShakeListener.BondConstraints>() {
            @Override
            public ShakeListener.BondConstraints makeAgent(ISpecies type) {
                return new ShakeListener.BondConstraints(new int[][]{{0, 2}, {1, 2}, {0, 1}}, new double[]{lOH, lOH, lHH});
            }

            @Override
            public void releaseAgent(ShakeListener.BondConstraints agent, ISpecies type) {

            }
        }, sim.getSpeciesManager());
        ShakeListener shake = new ShakeListener(sim.getSpeciesManager(), shakeAgents, integrator);
        shake.setMaxIterations(maxIterations);
        integrator.getEventManager().addListener(shake);
//        integrator.setOrientAtom((IAtomPositioned)((IMolecule)box.getMoleculeList(speciesOrient).getAtom(0)).getChildList().getAtom(0));
        integrator.setIsothermal(false);
        integrator.setTemperature(Kelvin.UNIT.toSim(273));
//        MeterTemperature meterTemperature = new MeterTemperature(box, 3);
//        System.out.println("T="+Kelvin.UNIT.fromSim(meterTemperature.getDataAsScalar()));
//        integrator.setThermostatInterval(100);
//        System.out.println("using rigid with dt="+dt);
//        System.out.println("h1 at "+((IAtomPositioned)box.getLeafList().getAtom(0)).getPosition());
//        System.out.println("o at "+((IAtomPositioned)box.getLeafList().getAtom(2)).getPosition());


        double chargeOxygen = Electron.UNIT.toSim(-0.82);
        double chargeHydrogen = Electron.UNIT.toSim(0.41);

        AtomType oType = species.getTypeByName("O");
        AtomType hType = species.getTypeByName("H");
        double epsOxygen = new P2WaterSPC(space).getEpsilon();
        double sigOxygen = new P2WaterSPC(space).getSigma();
        P2LennardJones potentialLJOO = new P2LennardJones(space, sigOxygen, epsOxygen);
        P2Electrostatic potentialQOO = new P2Electrostatic(space);
        potentialQOO.setCharge1(chargeOxygen);
        potentialQOO.setCharge2(chargeOxygen);
        potentialMaster.setPairPotential(oType, oType, new P2SoftSphericalSum(space, potentialLJOO, potentialQOO));

        P2Electrostatic potentialQHH = new P2Electrostatic(space);
        potentialQHH.setCharge1(chargeHydrogen);
        potentialQHH.setCharge2(chargeHydrogen);
        potentialMaster.setPairPotential(hType, hType, potentialQHH);


        P2Electrostatic potentialQOH = new P2Electrostatic(space);
        potentialQOH.setCharge1(chargeOxygen);
        potentialQOH.setCharge2(chargeHydrogen);
        potentialMaster.setPairPotential(oType, hType, potentialQOH);

        if (true) {
            sim.getController().setSleepPeriod(2);
            sim.getController().addActivity(new ActivityIntegrate(integrator));
            SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "SHAKE", 1);
            ((ColorSchemeByType) graphic.getDisplayBox(box).getColorScheme()).setColor(species.getTypeByName("H"), Color.WHITE);
            ((ColorSchemeByType) graphic.getDisplayBox(box).getColorScheme()).setColor(species.getTypeByName("O"), Color.RED);

            MeterEnergyFromIntegrator meterEnergy = new MeterEnergyFromIntegrator(integrator);
            AccumulatorHistory historyEnergy = new AccumulatorHistory(new HistoryScrolling(100));
            DataPumpListener pumpEnergy = new DataPumpListener(meterEnergy, historyEnergy);
            integrator.getEventManager().addListener(pumpEnergy);
            MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(integrator);
            AccumulatorHistory historyPE = new AccumulatorHistory(new HistoryScrolling(100));
            DataPumpListener pumpPE = new DataPumpListener(meterPE, historyPE);
            integrator.getEventManager().addListener(pumpPE);
            DisplayPlotXChart plotEnergy = new DisplayPlotXChart();
            plotEnergy.setLabel("energy");
            historyEnergy.setDataSink(plotEnergy.makeSink("energy"));
            historyPE.setDataSink(plotEnergy.makeSink("PE"));
            graphic.add(plotEnergy);

            return graphic;
        }
        sim.getController().runActivityBlocking(new ActivityIntegrate(integrator, Long.MAX_VALUE));
        return null;
    }

    public static void main(String[] args) {
        SimulationGraphic graphic = makeWaterDroplet();
        if (graphic != null) {
            graphic.makeAndDisplayFrame();
        }
    }
}
