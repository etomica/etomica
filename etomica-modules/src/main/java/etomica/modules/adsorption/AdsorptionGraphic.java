/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.adsorption;

 import java.awt.GridBagLayout;
import java.util.ArrayList;

import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.border.TitledBorder;

import etomica.action.IAction;
import etomica.space.Vector;
import etomica.atom.DiameterHashByType;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPump;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountTime;
import etomica.data.DataTag;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.graphics.DeviceButton;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.math.geometry.Plane;
import etomica.modifier.ModifierGeneral;
import etomica.space.Space;
 import etomica.space3d.Space3D;
import etomica.units.dimensions.Null;
import etomica.units.Pixel;
import etomica.units.Unit;
import etomica.data.history.HistoryCollapsingAverage;

/**
 * Catalysis graphical app.
 * Design by Lev Gelb
 * 
 * @author Andrew Schultz
 */
public class AdsorptionGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Adsorption";
    private final static int REPAINT_INTERVAL = 1;
    protected final DeviceThermoSlider tempSlider;
    protected final MeterProfileByVolumeAdsorption densityProfileMeterA, densityProfileMeterB;
    protected final Adsorption sim;

    public AdsorptionGraphic(final Adsorption simulation, Space _space) {

    	super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

        final IAction resetDataAction = new IAction() {
            public void actionPerformed() {
                getController().getResetAveragesButton().press();
            }
        };
        getController().getResetAveragesButton().setPostAction(new IAction() {
            public void actionPerformed() {
                sim.integratorHybrid.resetStepCount();
            }
        });

    	this.sim = simulation;

        Unit tUnit = Null.UNIT;

        getDisplayBox(sim.box).setPixelUnit(new Pixel(40/sim.box.getBoundary().getBoxSize().getX(1)));
        ((DisplayBoxCanvasG3DSys)getDisplayBox(sim.box).canvas).addPlane(new Plane(space, 0, 1, 0, sim.box.getBoundary().getBoxSize().getX(1)/2-0.00001));
        ((DiameterHashByType)getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.speciesA.getLeafType(), sim.p2AA.getCoreDiameter());
        Vector rMin = space.makeVector();
        Vector rMax = space.makeVector();
        double Ly = sim.box.getBoundary().getBoxSize().getX(1);
        double yMin = -0.5*Ly+0.5*sim.p1WallA.getSigma();
        double yMax = (0.5-0*sim.mcMoveIDA.getZFraction())*Ly - sim.p1WallA.getSigma();
        rMin.Ea1Tv1(-0.5, sim.box.getBoundary().getBoxSize());
        rMax.Ea1Tv1(0.5, sim.box.getBoundary().getBoxSize());
        rMax.setX(1, yMax);
        ((DisplayBoxCanvasG3DSys)getDisplayBox(sim.box).canvas).setBoundingBox(rMin, rMax);
        
        // Simulation Time
        final DisplayTextBox displayCycles = new DisplayTextBox();

        final DataSourceCountTime timer = new DataSourceCountTime(sim.integratorMD);
        displayCycles.setPrecision(7);
        displayCycles.setLabel("Simulation time");
        DataPumpListener pump= new DataPumpListener(timer,displayCycles, 100);
        sim.integratorMD.getEventManager().addListener(pump);

        //temperature selector
        tempSlider = new DeviceThermoSlider(sim.getController(), sim.integratorMC);
        tempSlider.setMinimum(0.0);
        tempSlider.setMaximum(3);
        tempSlider.setPrecision(1);
        tempSlider.setSliderMajorValues(3);
        tempSlider.setSliderMinorValues(2);
        tempSlider.setUnit(tUnit);
        tempSlider.setIsothermalButtonsVisibility(false);

        add(tempSlider);
        
        final EOSSW eos = new EOSSWPatel();
        final EOSSW eosB = new EOSSWPatel();
        MeterExcessAdsorbed meterCountAtoms = new MeterExcessAdsorbed(sim.speciesA, eos);
        meterCountAtoms.setBox(sim.box);
        meterCountAtoms.setRange(1, yMin, yMax);
        AccumulatorAverageCollapsing countAvg = new AccumulatorAverageCollapsing();
        DataPumpListener countPump = new DataPumpListener(meterCountAtoms, countAvg);
        countAvg.setPushInterval(100);
        sim.integratorHybrid.getEventManager().addListener(countPump);
        DisplayTextBoxesCAE displayCount = new DisplayTextBoxesCAE();
        displayCount.setAccumulator(countAvg);
        displayCount.setLabel("Excess adsorbed A");
        dataStreamPumps.add(countPump);

        MeterExcessAdsorbed meterCountAtomsB = new MeterExcessAdsorbed(sim.speciesB, eosB);
        meterCountAtomsB.setBox(sim.box);
        meterCountAtomsB.setRange(1, yMin, yMax);
        final AccumulatorAverageCollapsing countAvgB = new AccumulatorAverageCollapsing();
        final DataPumpListener countPumpB = new DataPumpListener(meterCountAtomsB, countAvgB);
        countAvgB.setPushInterval(100);
        final DisplayTextBoxesCAE displayCountB = new DisplayTextBoxesCAE();
        displayCountB.setLabel("Excess adsorbed B");
        dataStreamPumps.add(countPumpB);

        
        final AccumulatorHistory adsorbedHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        countAvg.addDataSink(adsorbedHistory, new StatType[]{countAvg.MOST_RECENT});
        final DisplayPlot adsorbedHistoryPlot = new DisplayPlot();
        adsorbedHistory.setDataSink(adsorbedHistoryPlot.getDataSet().makeDataSink());
        adsorbedHistoryPlot.setLabel("Adsorption");
        adsorbedHistoryPlot.getPlot().setTitle("Excess adsorbed per area");
        adsorbedHistoryPlot.setLegend(new DataTag[]{adsorbedHistory.getTag()}, "A");
        adsorbedHistory.setTimeDataSource(timer);

        final AccumulatorHistory adsorbedHistoryB = new AccumulatorHistory(new HistoryCollapsingAverage());
        countAvgB.addDataSink(adsorbedHistoryB, new StatType[]{countAvgB.MOST_RECENT});
        adsorbedHistoryB.setTimeDataSource(timer);
        adsorbedHistoryB.setDataSink(adsorbedHistoryPlot.getDataSet().makeDataSink());
        adsorbedHistoryPlot.setLegend(new DataTag[]{adsorbedHistoryB.getTag()}, "B");

        
        densityProfileMeterA = new MeterProfileByVolumeAdsorption(space);
        densityProfileMeterA.setRange(yMin, yMax);
        densityProfileMeterA.setBox(sim.box);
        MeterNMolecules meterNMolecules = new MeterNMolecules();
        meterNMolecules.setSpecies(sim.speciesA);
        densityProfileMeterA.setDataSource(meterNMolecules);
        AccumulatorAverageFixed densityProfileAvg = new AccumulatorAverageFixed(10);
        densityProfileAvg.setPushInterval(100);
        DataPumpListener profilePump = new DataPumpListener(densityProfileMeterA, densityProfileAvg, 10);
        sim.integratorHybrid.getEventManager().addListener(profilePump);
        dataStreamPumps.add(profilePump);

        densityProfileMeterB = new MeterProfileByVolumeAdsorption(space);
        densityProfileMeterB.setRange(yMin, yMax);
        densityProfileMeterB.setBox(sim.box);
        MeterNMolecules meterNMoleculesB = new MeterNMolecules();
        meterNMoleculesB.setSpecies(sim.speciesB);
        densityProfileMeterB.setDataSource(meterNMoleculesB);
        final AccumulatorAverageFixed densityProfileAvgB = new AccumulatorAverageFixed(10);
        densityProfileAvgB.setPushInterval(100);
        final DataPumpListener profilePumpB = new DataPumpListener(densityProfileMeterB, densityProfileAvgB, 10);
        dataStreamPumps.add(profilePumpB);

        final DisplayPlot profilePlot = new DisplayPlot();
        densityProfileAvg.addDataSink(profilePlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{densityProfileAvg.AVERAGE});
        densityProfileAvgB.addDataSink(profilePlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{densityProfileAvgB.AVERAGE});
        profilePlot.setLegend(new DataTag[]{densityProfileAvg.getTag()}, "A");
        profilePlot.setLegend(new DataTag[]{densityProfileAvgB.getTag()}, "B");
        profilePlot.setLabel("Profile");
        
        final JTabbedPane tabs = new JTabbedPane();
        getPanel().controlPanel.add(tabs,SimulationPanel.getVertGBC());
        JPanel tabA = new JPanel(new GridBagLayout());
        tabs.add(tabA, "A");

        final DeviceSlider pSlider = new DeviceSlider(sim.getController());
        pSlider.setShowBorder(true);
        pSlider.setBorderAlignment(TitledBorder.CENTER);
        pSlider.setShowValues(true);
        pSlider.setMinimum(-6);
        pSlider.setPrecision(2);
        pSlider.setMaximum(-1);
        pSlider.setNMajor(5);
        pSlider.setEditValues(true);
        final ModifierPMu pmu = new ModifierPMu(sim.p2AA, sim.integratorMC, eos, sim.mcMoveIDA, meterCountAtoms);
        pmu.setValue(-4);
        pSlider.setModifier(pmu);
        pSlider.setLabel("log10(P)");
        tabA.add(pSlider.graphic(), SimulationPanel.getVertGBC());
        
        final DisplayTextBox pSatDisplay = new DisplayTextBox();
        DataDouble.DataInfoDouble pSatInfo = new DataInfoDouble("P/Psat", Null.DIMENSION);
        pSatDisplay.putDataInfo(pSatInfo);

        final IAction updatePRatio = new IAction() {
            public void actionPerformed() {
                double temperature = sim.integratorMC.getTemperature();
                eos.setTemperature(temperature);
                double sigma = sim.p2AA.getCoreDiameter();
                double pSat = temperature < eos.Tc ? eos.pSat()/(sigma*sigma*sigma) : Double.NaN;
                double p = Math.pow(10, pSlider.getValue());
                data.x = p/pSat;
                pSatDisplay.putData(data);
            }
            DataDouble data = new DataDouble();
        };
        updatePRatio.actionPerformed();

        tabA.add(pSatDisplay.graphic(),SimulationPanel.getVertGBC());
      
        final DataSourceDensityFunction nominalDensity = new DataSourceDensityFunction(eos, yMin, yMax);
        nominalDensity.setPressure(Math.pow(10, pSlider.getValue()));
        profilePlot.setLegend(new DataTag[]{nominalDensity.getTag()}, "A (bulk)");
        DataPumpListener pumpNominalDensity = new DataPumpListener(nominalDensity, profilePlot.getDataSet().makeDataSink(), 100);
        sim.integratorHybrid.getEventManager().addListener(pumpNominalDensity);

        final DataSourceDensityFunction nominalLiquidDensity = new DataSourceDensityFunction(eos, yMin, yMax, true);
        nominalLiquidDensity.setPressure(Math.pow(10, pSlider.getValue()));
        profilePlot.setLegend(new DataTag[]{nominalLiquidDensity.getTag()}, "Aliq (bulk)");
        DataPumpListener pumpNominalLiquidDensity = new DataPumpListener(nominalLiquidDensity, profilePlot.getDataSet().makeDataSink(), 100);
        if (false) sim.integratorHybrid.getEventManager().addListener(pumpNominalLiquidDensity);
        pSlider.setPostAction(new IAction() {
            public void actionPerformed() {
                updatePRatio.actionPerformed();
                double p = Math.pow(10, pSlider.getValue());
                nominalDensity.setPressure(p);
                nominalLiquidDensity.setPressure(p);
            }
        });

        final DeviceSlider epsSlider = new DeviceSlider(sim.getController());
        epsSlider.setShowBorder(true);
        epsSlider.setBorderAlignment(TitledBorder.CENTER);
        epsSlider.setShowValues(true);
        epsSlider.setMinimum(0);
        epsSlider.setPrecision(1);
        epsSlider.setMaximum(10);
        epsSlider.setNMajor(5);
        final ModifierGeneral epsModifier = new ModifierGeneral(sim.p1WallA, "epsilon");
        epsSlider.setModifier(epsModifier);
        epsSlider.setLabel("epsilon");
        tabA.add(epsSlider.graphic(), SimulationPanel.getVertGBC());

        final JPanel tabB = new JPanel(new GridBagLayout());
        tabs.add(tabB, "B");
        final JPanel tabBB = new JPanel(new GridBagLayout());

        final DeviceSlider pSliderB = new DeviceSlider(sim.getController());
        pSliderB.setEnabled(false);
        pSliderB.setShowBorder(true);
        pSliderB.setBorderAlignment(TitledBorder.CENTER);
        pSliderB.setShowValues(true);
        pSliderB.setMinimum(-6);
        pSliderB.setPrecision(2);
        pSliderB.setMaximum(-1);
        pSliderB.setNMajor(5);
        pSliderB.setEditValues(true);
        final ModifierPMu pmuB = new ModifierPMu(sim.p2BB, sim.integratorMC, eosB, sim.mcMoveIDB, meterCountAtomsB);
        pmuB.setValue(-4);
        pSliderB.setModifier(pmuB);
        pSliderB.setLabel("log10(P)");
        tabBB.add(pSliderB.graphic(), SimulationPanel.getVertGBC());

        final DisplayTextBox pSatBDisplay = new DisplayTextBox();
        DataDouble.DataInfoDouble pSatBInfo = new DataInfoDouble("P/Psat", Null.DIMENSION);
        pSatBDisplay.putDataInfo(pSatBInfo);

        final IAction updatePRatioB = new IAction() {
            public void actionPerformed() {
                double temperature = sim.integratorMC.getTemperature();
                eosB.setTemperature(temperature);
                double sigma = sim.p2BB.getCoreDiameter();
                double pSat = temperature < eosB.Tc ? eosB.pSat()/(sigma*sigma*sigma) : Double.NaN;
                double p = Math.pow(10, pSliderB.getValue());
                data.x = p/pSat;
                pSatBDisplay.putData(data);
            }
            DataDouble data = new DataDouble();
        };

        tabBB.add(pSatBDisplay.graphic(),SimulationPanel.getVertGBC());

        final DataSourceDensityFunction nominalDensityB = new DataSourceDensityFunction(eosB, yMin, yMax);
        nominalDensityB.setPressure(Math.pow(10, pSliderB.getValue()));
        pSliderB.setPostAction(new IAction() {
            public void actionPerformed() {
                updatePRatioB.actionPerformed();
                double p = Math.pow(10, pSliderB.getValue());
                nominalDensityB.setPressure(p);
            }
        });
        final DataPumpListener pumpNominalDensityB = new DataPumpListener(nominalDensityB, profilePlot.getDataSet().makeDataSink(), 100);
        profilePlot.setLegend(new DataTag[]{nominalDensityB.getTag()}, "B (bulk)");

        final DeviceSlider epsSliderB = new DeviceSlider(sim.getController());
        epsSliderB.setEnabled(false);
        epsSliderB.setShowBorder(true);
        epsSliderB.setBorderAlignment(TitledBorder.CENTER);
        epsSliderB.setShowValues(true);
        epsSliderB.setMinimum(0);
        epsSliderB.setPrecision(1);
        epsSliderB.setMaximum(10);
        epsSliderB.setNMajor(5);
        final ModifierGeneral epsModifierB = new ModifierGeneral(sim.p1WallB, "epsilon");
        epsSliderB.setModifier(epsModifierB);
        epsSliderB.setLabel("epsilon");
        tabBB.add(epsSliderB.graphic(), SimulationPanel.getVertGBC());

        
        final DeviceButton mixButton = new DeviceButton(sim.getController());
        mixButton.setLabel("Add Component B");
        mixButton.setAction(new IAction() {
            public void actionPerformed() {
                if (!pSliderB.isEnabled()) {
                    pSliderB.setEnabled(true);
                    epsSliderB.setEnabled(true);
                    sim.integratorHybrid.getEventManager().addListener(pumpNominalDensityB);
                    sim.integratorHybrid.getEventManager().addListener(profilePumpB);
                    sim.integratorHybrid.getEventManager().addListener(countPumpB);
                    sim.integratorMC.getMoveManager().addMCMove(sim.mcMoveIDB);
                    displayCountB.setAccumulator(countAvgB);
                    updatePRatioB.actionPerformed();

                    mixButton.setLabel("Remove Component B");
                }
                else {
                    pSliderB.setEnabled(false);
                    epsSliderB.setEnabled(false);
                    sim.integratorHybrid.getEventManager().removeListener(pumpNominalDensityB);
                    sim.integratorHybrid.getEventManager().removeListener(profilePumpB);
                    sim.integratorHybrid.getEventManager().removeListener(countPumpB);
                    sim.integratorMC.getMoveManager().removeMCMove(sim.mcMoveIDB);
                    displayCountB.setAccumulator(null);
                    profilePlot.getDataSet().reset();
                    adsorbedHistoryPlot.getDataSet().reset();
                    sim.box.setNMolecules(sim.speciesB, 0);
                    sim.integratorHybrid.reset();
                    mixButton.setLabel("Add Component B");
                }
            }
        });
        tabB.add(mixButton.graphic(), SimulationPanel.getVertGBC());


        tempSlider.setSliderPostAction(new IAction() {
            public void actionPerformed() {
                updatePRatio.actionPerformed();
                // force pmu to update.  the slider value hasn't changed,
                // but the effect of that value has.
                pmu.setValue(pmu.getValue());

                updatePRatioB.actionPerformed();
                // force pmu to update.  the slider value hasn't changed,
                // but the effect of that value has.
                pmuB.setValue(pmuB.getValue());
            }
        });

        tabB.add(tabBB, SimulationPanel.getVertGBC());
        
        MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integratorMD);
        final AccumulatorAverageCollapsing peAvg = new AccumulatorAverageCollapsing();
        final DataPumpListener pePump = new DataPumpListener(meterPE, peAvg);
        sim.integratorHybrid.getEventManager().addListener(pePump);
        dataStreamPumps.add(pePump);
        DisplayTextBoxesCAE peDisplay = new DisplayTextBoxesCAE();
        peDisplay.setAccumulator(peAvg);

//		MeterTemperature thermometer = new MeterTemperature(sim, sim.box, space.D());
//        final AccumulatorHistory temperatureHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
//        final DataPumpListener temperaturePump = new DataPumpListener(thermometer,temperatureHistory, 1000);
//        sim.integratorMD.getEventManager().addListener(temperaturePump);
//		dataStreamPumps.add(temperaturePump);
//        DisplayPlot temperatureHistoryPlot = new DisplayPlot();
//        temperatureHistory.setDataSink(temperatureHistoryPlot.getDataSet().makeDataSink());
//        temperatureHistoryPlot.setLabel("Temperature");
//        temperatureHistoryPlot.setDoLegend(false);
//        temperatureHistoryPlot.setUnit(Kelvin.UNIT);
//        temperatureHistoryPlot.getPlot().setYLabel("Temperature (K)");
        
        add(displayCycles);
        add(displayCount);
        add(displayCountB);
        add(peDisplay);
        add(adsorbedHistoryPlot);
        add(profilePlot);
//        add(temperatureHistoryPlot);
    }

    public static void main(String[] args) {
        Space space = Space3D.getInstance();


        AdsorptionGraphic adsGraphic = new AdsorptionGraphic(new Adsorption(space), space);
		SimulationGraphic.makeAndDisplayFrame
		        (adsGraphic.getPanel(), APP_NAME);
    }
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
	        getRootPane().putClientProperty(
	                        "defeatSystemEventQueueCheck", Boolean.TRUE);
            Space sp = Space3D.getInstance();
            AdsorptionGraphic swmdGraphic = new AdsorptionGraphic(new Adsorption(sp), sp);

		    getContentPane().add(swmdGraphic.getPanel());
	    }

        private static final long serialVersionUID = 1L;
    }

}
