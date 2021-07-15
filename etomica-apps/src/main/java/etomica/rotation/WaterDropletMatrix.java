/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.rotation;


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
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorRigidMatrixIterative;
import etomica.models.water.ConformationWaterTIP4P;
import etomica.models.water.OrientationCalcWater4P;
import etomica.models.water.P2WaterTIP4PSoft;
import etomica.models.water.SpeciesWater4P;
import etomica.molecule.MoleculePositionCOM;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesGeneral;
import etomica.units.Kelvin;
import etomica.util.Constants;

import java.awt.*;

public class WaterDropletMatrix {

    public static SimulationGraphic makeWaterDroplet() {
        Space space = Space3D.getInstance();
        Simulation sim = new Simulation(space);
        OrientationCalcWater4P calcer = new OrientationCalcWater4P(sim.getSpace());
        ConformationWaterTIP4P conformationWater = new ConformationWaterTIP4P(space);
        SpeciesGeneral species = SpeciesWater4P.create(true, conformationWater, true);
        sim.addSpecies(species);
        Box box = new Box(new BoundaryRectangularNonperiodic(sim.getSpace()), space);
        sim.addBox(box);
        box.setNMolecules(species, 108);
        box.setDensity(0.7 / 18.0 * Constants.AVOGADRO / 1E24);
        ConfigurationWater108TIP4P config = new ConfigurationWater108TIP4P();
//        Configuration config = new ConfigurationLattice(new LatticeCubicFcc(), space);
        PotentialMaster potentialMaster = new PotentialMaster();
        double timeInterval = 0.002;
        int maxIterations = 20;
        IntegratorRigidMatrixIterative integrator = new IntegratorRigidMatrixIterative(sim.getSpeciesManager(), sim.getRandom(), potentialMaster, timeInterval, 1, box);
        integrator.printInterval = 100;
        integrator.setMaxIterations(maxIterations);
        config.initializeCoordinates(box);
        integrator.setOrientationCalc(species, calcer);
        integrator.setTemperature(Kelvin.UNIT.toSim(298));
//        integrator.setIsothermal(true);
        integrator.setThermostatInterval(100);

        P2WaterTIP4PSoft p2Water = new P2WaterTIP4PSoft(sim.getSpace(), Double.POSITIVE_INFINITY, new MoleculePositionCOM(space));
        potentialMaster.addPotential(p2Water, new ISpecies[]{species, species});

        if (false) {
            sim.getController().setSleepPeriod(2);
            sim.getController().addActivity(new ActivityIntegrate(integrator));
            SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "Matrix", 1);
            ((ColorSchemeByType) graphic.getDisplayBox(box).getColorScheme()).setColor(species.getTypeByName("H"), Color.WHITE);
            ((ColorSchemeByType) graphic.getDisplayBox(box).getColorScheme()).setColor(species.getTypeByName("O"), Color.RED);

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
        sim.getController().runActivityBlocking(new ActivityIntegrate(integrator, Long.MAX_VALUE));
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
