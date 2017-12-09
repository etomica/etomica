/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.rotation;

import etomica.action.PDBWriter;
import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPump;
import etomica.data.DataSourceCountTime;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergyFromIntegrator;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorRigidIterative;
import etomica.listener.IntegratorListenerAction;
import etomica.models.water.OrientationCalcWater3P;
import etomica.models.water.P2WaterSPCSoft;
import etomica.models.water.SpeciesWater3POriented;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.units.Kelvin;
import etomica.util.Constants;

import java.awt.*;

public class WaterDroplet {

    public static SimulationGraphic makeWaterDroplet() {
        Space space = Space3D.getInstance();
        Simulation sim = new Simulation(space);
        Box box = new Box(new BoundaryRectangularNonperiodic(sim.getSpace()), space);
        sim.addBox(box);
        SpeciesWater3POriented species = new SpeciesWater3POriented(sim.getSpace(), true);
        sim.addSpecies(species);
        box.setNMolecules(species, 108);
        box.setDensity(1/18.0*Constants.AVOGADRO/1E24);
        ConfigurationWater108 configFile = new ConfigurationWater108();
        configFile.initializeCoordinates(box);
        PotentialMaster potentialMaster = new PotentialMaster();
        double timeInterval = 0.002;
        int maxIterations = 20;
        IntegratorRigidIterative integrator = new IntegratorRigidIterative(sim, potentialMaster, timeInterval, 1, space);
        integrator.printInterval = 100;
        integrator.setMaxIterations(maxIterations);
        integrator.setBox(box);
        integrator.setIsothermal(true);
        OrientationCalcWater3P calcer = new OrientationCalcWater3P(sim.getSpace());
        integrator.setOrientationCalc(species, calcer);
        integrator.setTemperature(Kelvin.UNIT.toSim(220));
        integrator.setThermostatInterval(1000);
        ActivityIntegrate ai = new ActivityIntegrate(integrator);
        sim.getController().addAction(ai);
//        System.out.println("h1 at "+((IAtomPositioned)box.getLeafList().getAtom(0)).getPosition());
//        System.out.println("o at "+((IAtomPositioned)box.getLeafList().getAtom(2)).getPosition());

        P2WaterSPCSoft p2Water = new P2WaterSPCSoft(sim.getSpace());

        potentialMaster.addPotential(p2Water, new ISpecies[]{species,species});
        if (false) {
            ai.setSleepPeriod(2);
            SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "Rigid", 1);
            ((ColorSchemeByType)graphic.getDisplayBox(box).getColorScheme()).setColor(species.getHydrogenType(), Color.WHITE);
            ((ColorSchemeByType)graphic.getDisplayBox(box).getColorScheme()).setColor(species.getOxygenType(), Color.RED);
    
            MeterEnergy meterE = new MeterEnergy(potentialMaster, box);
            meterE.setKinetic(new MeterKineticEnergyFromIntegrator(integrator));
            meterE.setPotential(new MeterPotentialEnergyFromIntegrator(integrator));
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
        sim.getController().actionPerformed();
        return null;
    }

    public static void main(String[] args) {
        SimulationGraphic graphic = makeWaterDroplet();
        if (graphic != null) {
            graphic.makeAndDisplayFrame();
        }
    }
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
            SimulationGraphic graphic = makeWaterDroplet();

            getContentPane().add(graphic.getPanel());
        }

        private static final long serialVersionUID = 1L;
    }
}
