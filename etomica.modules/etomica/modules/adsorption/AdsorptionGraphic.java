package etomica.modules.adsorption;

 import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.util.ArrayList;

import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import etomica.action.IAction;
import etomica.api.IVectorMutable;
import etomica.atom.DiameterHashByType;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPump;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountTime;
import etomica.data.meter.MeterNMolecules;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.math.geometry.Plane;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.units.Null;
import etomica.units.Pixel;
import etomica.units.Unit;
import etomica.util.HistoryCollapsingAverage;

/**
 * Catalysis graphical app.
 * Design by Lev Gelb
 * 
 * @author Andrew Schultz
 */
public class AdsorptionGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Adsorption";
    private final static int REPAINT_INTERVAL = 1;
    protected DeviceThermoSlider tempSlider;
    protected final MeterProfileByVolumeAdsorption densityProfileMeter;
    protected Adsorption sim;

    public AdsorptionGraphic(final Adsorption simulation, ISpace _space) {

    	super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL, _space, simulation.getController());

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

        final IAction resetDataAction = new IAction() {
            public void actionPerformed() {
                getController().getResetAveragesButton().press();
            }
        };

    	this.sim = simulation;

        Unit tUnit = Null.UNIT;

        getDisplayBox(sim.box).setPixelUnit(new Pixel(40/sim.box.getBoundary().getBoxSize().getX(1)));
        ((DisplayBoxCanvasG3DSys)getDisplayBox(sim.box).canvas).addPlane(new Plane(space, 0, 1, 0, sim.box.getBoundary().getBoxSize().getX(1)/2-0.00001));
        ((DiameterHashByType)getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species.getLeafType(), sim.p2.getCoreDiameter());
        IVectorMutable rMin = space.makeVector();
        IVectorMutable rMax = space.makeVector();
        double Ly = sim.box.getBoundary().getBoxSize().getX(1);
        double yMin = -0.5*Ly+0.5*sim.p1Wall.getSigma();
        double yMax = (0.5-sim.mcMoveID.getZFraction())*Ly - sim.p1Wall.getSigma();
        rMin.Ea1Tv1(-0.5, sim.box.getBoundary().getBoxSize());
        rMax.Ea1Tv1(0.5, sim.box.getBoundary().getBoxSize());
        rMax.setX(1, yMax);
        ((DisplayBoxCanvasG3DSys)getDisplayBox(sim.box).canvas).setBoundingBox(rMin, rMax);
        
        // Simulation Time
        final DisplayTextBox displayCycles = new DisplayTextBox();

        final DataSourceCountTime meterCycles = new DataSourceCountTime(sim.integratorMD);
        displayCycles.setPrecision(7);
        DataPumpListener pump= new DataPumpListener(meterCycles,displayCycles, 100);
        sim.integratorMD.getEventManager().addListener(pump);
        displayCycles.setLabel("Simulation time");

        //temperature selector
        tempSlider = new DeviceThermoSlider(sim.getController(), sim.integratorMC);
//        tempSlider.setUnit(Kelvin.UNIT);
        tempSlider.setMinimum(0.0);
        tempSlider.setMaximum(5);
        tempSlider.setPrecision(1);
        tempSlider.setSliderMajorValues(5);
        tempSlider.setUnit(tUnit);
        tempSlider.setIsothermalButtonsVisibility(false);

        JPanel statePanel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc2 = new GridBagConstraints();
        gbc2.gridx = 0;  gbc2.gridy = 0;
//        statePanel.add(tempSlider.graphic(), gbc2);
        
        add(tempSlider);
        
        final EOSSW eos = new EOSSWPatel();
        MeterExcessAdsorbed meterCountAtoms = new MeterExcessAdsorbed(eos);
        meterCountAtoms.setBox(sim.box);
        meterCountAtoms.setRange(1, yMin, yMax);
        AccumulatorAverageCollapsing countAvg = new AccumulatorAverageCollapsing();
        DataPumpListener countPump = new DataPumpListener(meterCountAtoms, countAvg);
        countAvg.setPushInterval(100);
        sim.integratorHybrid.getEventManager().addListener(countPump);
        DisplayTextBoxesCAE displayCount = new DisplayTextBoxesCAE();
        displayCount.setAccumulator(countAvg);
        dataStreamPumps.add(countPump);
        
        final AccumulatorHistory adsorbedHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        countAvg.addDataSink(adsorbedHistory, new StatType[]{countAvg.MOST_RECENT});
        DisplayPlot adsorbedHistoryPlot = new DisplayPlot();
        adsorbedHistory.setDataSink(adsorbedHistoryPlot.getDataSet().makeDataSink());
        adsorbedHistoryPlot.setLabel("Adsorption");
        adsorbedHistoryPlot.getPlot().setTitle("Excess adsorbed per area");
        adsorbedHistoryPlot.setDoLegend(false);
        DataSourceCountTime timer = new DataSourceCountTime(sim.integratorMD);
        adsorbedHistory.setTimeDataSource(timer);

        densityProfileMeter = new MeterProfileByVolumeAdsorption(space);
        densityProfileMeter.setRange(yMin, yMax);
        densityProfileMeter.setBox(sim.box);
        MeterNMolecules meterNMolecules = new MeterNMolecules();
        meterNMolecules.setSpecies(sim.species);
        densityProfileMeter.setDataSource(meterNMolecules);
        AccumulatorAverageFixed densityProfileAvg = new AccumulatorAverageFixed(10);
        densityProfileAvg.setPushInterval(100);
        DataPumpListener profilePump = new DataPumpListener(densityProfileMeter, densityProfileAvg, 10);
        sim.integratorHybrid.getEventManager().addListener(profilePump);
        dataStreamPumps.add(profilePump);

        DisplayPlot profilePlot = new DisplayPlot();
        densityProfileAvg.addDataSink(profilePlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{densityProfileAvg.AVERAGE});
        profilePlot.setDoLegend(false);
        profilePlot.setLabel("Profile");

        final DeviceSlider pSlider = new DeviceSlider(sim.getController());
        pSlider.setShowBorder(true);
        pSlider.setBorderAlignment(TitledBorder.CENTER);
        pSlider.setShowValues(true);
        pSlider.setMinimum(-4);
        pSlider.setPrecision(1);
        pSlider.setMaximum(0);
        pSlider.setNMajor(5);
        final ModifierPMu pmu = new ModifierPMu(sim.p2, sim.integratorMC, eos, sim.mcMoveID, meterCountAtoms);
        pmu.setValue(-4);
        pSlider.setModifier(pmu);
        pSlider.setLabel("log10(P)");
        add(pSlider);
        
        final DisplayTextBox pSatDisplay = new DisplayTextBox();
        DataDouble.DataInfoDouble pSatInfo = new DataInfoDouble("P/Psat", Null.DIMENSION);
        pSatDisplay.putDataInfo(pSatInfo);

        final IAction updatePRatio = new IAction() {
            public void actionPerformed() {
                double temperature = sim.integratorMC.getTemperature();
                eos.setTemperature(temperature);
                double sigma = sim.p2.getCoreDiameter();
                double pSat = temperature < eos.Tc ? eos.pSat()/(sigma*sigma*sigma) : Double.NaN;
                double p = Math.pow(10, pSlider.getValue());
                data.x = p/pSat;
                pSatDisplay.putData(data);
            }
            DataDouble data = new DataDouble();
        };
        tempSlider.setSliderPostAction(new IAction() {
            public void actionPerformed() {
                updatePRatio.actionPerformed();
                pmu.setValue(pmu.getValue());
            }
        });
        updatePRatio.actionPerformed();
        
        getPanel().controlPanel.add(pSatDisplay.graphic(),SimulationPanel.getVertGBC());
      
        final DataSourceDensityFunction nominalDensity = new DataSourceDensityFunction(eos, yMin, yMax);
        nominalDensity.setPressure(Math.pow(10, pSlider.getValue()));
        pSlider.setPostAction(new IAction() {
            public void actionPerformed() {
                updatePRatio.actionPerformed();
                double p = Math.pow(10, pSlider.getValue());
                nominalDensity.setPressure(p);
            }
        });
        DataPumpListener pumpNominalDensity = new DataPumpListener(nominalDensity, profilePlot.getDataSet().makeDataSink(), 100);
        sim.integratorHybrid.getEventManager().addListener(pumpNominalDensity);
        
//
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
        add(adsorbedHistoryPlot);
        add(profilePlot);
//        add(temperatureHistoryPlot);
    }

    public static void main(String[] args) {
        ISpace space = Space3D.getInstance();


        AdsorptionGraphic adsGraphic = new AdsorptionGraphic(new Adsorption(space), space);
		SimulationGraphic.makeAndDisplayFrame
		        (adsGraphic.getPanel(), APP_NAME);
    }
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
	        getRootPane().putClientProperty(
	                        "defeatSystemEventQueueCheck", Boolean.TRUE);
            ISpace sp = Space3D.getInstance();
            AdsorptionGraphic swmdGraphic = new AdsorptionGraphic(new Adsorption(sp), sp);

		    getContentPane().add(swmdGraphic.getPanel());
	    }

        private static final long serialVersionUID = 1L;
    }

}
