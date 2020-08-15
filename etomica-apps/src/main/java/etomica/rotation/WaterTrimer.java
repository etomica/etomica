/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.rotation;


import etomica.action.activity.ActivityIntegrate2;
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
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorRigidIterative;
import etomica.models.water.OrientationCalcWater3P;
import etomica.models.water.P2WaterSPCSoft;
import etomica.models.water.SpeciesWater3POriented;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.units.Kelvin;
import etomica.util.Constants;

import java.awt.*;

public class WaterTrimer {

    public static SimulationGraphic makeWaterDroplet() {
        Space space = Space3D.getInstance();
        Simulation sim = new Simulation(space);
        SpeciesWater3POriented species = new SpeciesWater3POriented(sim.getSpace(), true);
        sim.addSpecies(species);
        Box box = new Box(new BoundaryRectangularNonperiodic(sim.getSpace()), space);
        sim.addBox(box);
        box.setNMolecules(species, 3);
        box.setDensity(0.9 / 18.0 * Constants.AVOGADRO / 1E24);
        ConfigurationWater3_3P config = new ConfigurationWater3_3P();
//        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(), space);
        config.initializeCoordinates(box);
        box.getBoundary().setBoxSize(Vector.of(new double[]{15, 15, 15}));
        PotentialMaster potentialMaster = new PotentialMaster();
        double timeInterval = 0.001;
        int maxIterations = 20;
        IntegratorRigidIterative integrator = new IntegratorRigidIterative(sim, potentialMaster, timeInterval, 1, box);
        integrator.printInterval = 1000;
        integrator.setMaxIterations(maxIterations);
//        integrator.setIsothermal(true);
        OrientationCalcWater3P calcer = new OrientationCalcWater3P(sim.getSpace());
        integrator.setOrientationCalc(species, calcer);
        integrator.setTemperature(Kelvin.UNIT.toSim(0));
        integrator.setThermostatInterval(100);
//        System.out.println("h1 at "+((IAtomPositioned)box.getLeafList().getAtom(0)).getPosition());
//        System.out.println("o at "+((IAtomPositioned)box.getLeafList().getAtom(2)).getPosition());

        P2WaterSPCSoft p2Water = new P2WaterSPCSoft(sim.getSpace());

        potentialMaster.addPotential(p2Water, new ISpecies[]{species, species});
//        WriteConfiguration writeConfig = new WriteConfiguration(space);
//        writeConfig.setBox(box);
//        writeConfig.setConfName("water3_3P");
//        writeConfig.setDoApplyPBC(false);
//        integrator.addIntervalAction(writeConfig);
//        integrator.setActionInterval(writeConfig, 1000);
        if (false) {
            sim.getController2().addActivity(new ActivityIntegrate2(integrator)).setSleepPeriod(2);
            SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "Rigid", 1);
            ((ColorSchemeByType) graphic.getDisplayBox(box).getColorScheme()).setColor(species.getHydrogenType(), Color.WHITE);
            ((ColorSchemeByType) graphic.getDisplayBox(box).getColorScheme()).setColor(species.getOxygenType(), Color.RED);

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
//        PDBWriter writePDB = new PDBWriter();
//        writePDB.setBox(box);
//        writePDB.setFileName("water108Eq");
//        integrator.addIntervalAction(writePDB);
//        integrator.setActionInterval(writePDB, 10000);
        sim.getController2().runActivityBlocking(new ActivityIntegrate2(integrator), Long.MAX_VALUE);
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
