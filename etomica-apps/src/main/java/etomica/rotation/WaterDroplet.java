/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.rotation;

import etomica.action.PDBWriter;
import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPump;
import etomica.data.DataSourceCountTime;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterEnergyFromIntegrator;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorRigidIterative;
import etomica.models.water.OrientationCalcWater3P;
import etomica.models.water.P2WaterSPC;
import etomica.models.water.SpeciesWater3P;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.util.Constants;

import java.awt.*;

public class WaterDroplet {

    public static SimulationGraphic makeWaterDroplet() {
        Space space = Space3D.getInstance();
        Simulation sim = new Simulation(space);
        SpeciesGeneral species = SpeciesWater3P.create(true);
        sim.addSpecies(species);
        Box box = new Box(new BoundaryRectangularNonperiodic(sim.getSpace()), space);
        sim.addBox(box);
        box.setNMolecules(species, 108);
        box.setDensity(1 / 18.0 * Constants.AVOGADRO / 1E24);
        ConfigurationWater108 configFile = new ConfigurationWater108();
        configFile.initializeCoordinates(box);
        PotentialMaster potentialMaster = new PotentialMaster(sim.getSpeciesManager(), box, BondingInfo.noBonding());
        double timeInterval = 0.002;
        int maxIterations = 20;
        IntegratorRigidIterative integrator = new IntegratorRigidIterative(sim.getSpeciesManager(), sim.getRandom(), potentialMaster, timeInterval, 1, box);
        integrator.printInterval = 100;
        integrator.setMaxIterations(maxIterations);
        integrator.setIsothermal(true);
        OrientationCalcWater3P calcer = new OrientationCalcWater3P(sim.getSpace());
        integrator.setOrientationCalc(species, calcer);
        integrator.setTemperature(Kelvin.UNIT.toSim(220));
        integrator.setThermostatInterval(1000);
//        System.out.println("h1 at "+((IAtomPositioned)box.getLeafList().getAtom(0)).getPosition());
//        System.out.println("o at "+((IAtomPositioned)box.getLeafList().getAtom(2)).getPosition());

        double chargeOxygen = Electron.UNIT.toSim(-0.82);
        double chargeHydrogen = Electron.UNIT.toSim(0.41);

        AtomType oType = species.getTypeByName("O");
        AtomType hType = species.getTypeByName("H");
        double epsOxygen = P2WaterSPC.epsilonOO;
        double sigOxygen = P2WaterSPC.sigmaOO;
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
            SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "Rigid Fasterer", 1);
            ((ColorSchemeByType) graphic.getDisplayBox(box).getColorScheme()).setColor(species.getTypeByName("H"), Color.WHITE);
            ((ColorSchemeByType) graphic.getDisplayBox(box).getColorScheme()).setColor(species.getTypeByName("O"), Color.RED);

            MeterEnergyFromIntegrator meterE = new MeterEnergyFromIntegrator(integrator);
            AccumulatorHistory history = new AccumulatorHistory(new HistoryCollapsingAverage());
            history.setTimeDataSource(new DataSourceCountTime(integrator));
            DataPump pump = new DataPump(meterE, history);
            DisplayPlot ePlot = new DisplayPlot();
            history.setDataSink(ePlot.getDataSet().makeDataSink());
            IntegratorListenerAction pumpListener = new IntegratorListenerAction(pump);
            pumpListener.setInterval(10);
            integrator.getEventManager().addListener(pumpListener);
            ePlot.setLabel("Energy");
            graphic.add(ePlot);
            return graphic;
        }
        WriteConfiguration writeConfig = new WriteConfiguration(space);
        writeConfig.setBox(box);
        writeConfig.setConfName("water108Eq");
        writeConfig.setDoApplyPBC(false);
        IntegratorListenerAction writeConfigListener = new IntegratorListenerAction(writeConfig);
        writeConfigListener.setInterval(10000);
        integrator.getEventManager().addListener(writeConfigListener);
        PDBWriter writePDB = new PDBWriter();
        writePDB.setBox(box);
        writePDB.setFileName("water108Eq");
        IntegratorListenerAction writePDBListener = new IntegratorListenerAction(writePDB);
        writePDBListener.setInterval(10000);
        integrator.getEventManager().addListener(writePDBListener);
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
