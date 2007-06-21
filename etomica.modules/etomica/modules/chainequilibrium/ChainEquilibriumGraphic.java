package etomica.modules.chainequilibrium;

import java.awt.GridBagConstraints;
import java.util.ArrayList;

import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.SwingConstants;
import javax.swing.border.TitledBorder;

import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.DataSinkTable;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceThermoSelector;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTable;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.units.PrefixedUnit;
import etomica.util.Constants.CompassDirection;

/**
 * @author William Scharmach
 * 
 * TODO To change the template for this generated type comment go to Window -
 * Preferences - Java - Code Style - Code Templates
 */
public class ChainEquilibriumGraphic extends SimulationGraphic {

	private static final String APP_NAME = "Chain Reaction Equilibrium";
	private static final int REPAINT_INTERVAL = 35;

    protected ChainEquilibriumSim sim;

    public ChainEquilibriumGraphic(ChainEquilibriumSim simulation) {

		super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);
        this.sim = simulation;
        
        ArrayList dataStreamPumps = getController().getDataStreamPumps();
        
        getDisplayPhase(sim.phase).setPixelUnit(new Pixel(10));

        GridBagConstraints horizGBC = SimulationPanel.getHorizGBC();
        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        // ********* Data Declaration Section *******	
        int eMin = 0, eMax = 1000,  majorSpacing = 15, minorSpacing = 5;

        // **** Stuff that Modifies the Simulation

        // Sliders on Well depth page
        final DeviceSlider AASlider = sliders(sim, eMin, eMax, "RR epsilon", sim.AAbonded);
        final DeviceSlider ABSlider = sliders(sim, eMin, eMax, "RB epsilon", sim.ABbonded);
        final DeviceSlider BBSlider = sliders(sim, eMin, eMax, "BB epsilon", sim.BBbonded);
        
        // Sliders on Core size page
        DeviceSlider AAWellSlider = sliders(sim, 5, 95, "RR core", majorSpacing, minorSpacing, sim.AAbonded);
        DeviceSlider ABWellSlider = sliders(sim, 5, 95, "RB core", majorSpacing, minorSpacing, sim.ABbonded);
        DeviceSlider BBWellSlider = sliders(sim, 5, 95, "BB core", majorSpacing, minorSpacing, sim.BBbonded);

        // The Species Editors
        MySpeciesEditor AEditor = new MySpeciesEditor(this, sim.speciesA.getAgent(sim.phase), "Red");
        MySpeciesEditor BEditor = new MySpeciesEditor(this, sim.speciesB.getAgent(sim.phase), "Black");
		
        // the Atom Diameter Modifer
        DiameterModifier sizeModifier = new DiameterModifier(sim.AAbonded,sim.ABbonded, sim.BBbonded, sim.speciesA, sim.speciesB);

        DeviceSlider sizeSlider = new DeviceSlider(sim.getController(), sizeModifier);
		
        DisplayBox tBox = new DisplayBox();

        JPanel speciesEditors = new JPanel(new java.awt.GridLayout(0, 1));
        JPanel epsilonSliders = new JPanel(new java.awt.GridLayout(0, 1));
        JPanel temperaturePanel = new JPanel(new java.awt.GridBagLayout());
        JPanel sizeSliders = new JPanel(new java.awt.GridLayout(0, 1));
        JPanel lambdaSliders = new JPanel(new java.awt.GridLayout(0, 1));
        temperaturePanel.setBorder(new javax.swing.border.TitledBorder("Temperature (K)"));

        DataPump tPump = new DataPump (sim.thermometer, tBox);
        dataStreamPumps.add(tPump);

        AccumulatorAverage accumulator = new AccumulatorAverage();
        DataPump pump = new DataPump(sim.molecularCount,accumulator);
        dataStreamPumps.add(pump);
        sim.integratorHard1.addIntervalAction(pump);

        DataSinkTable dataTable = new DataSinkTable();
        accumulator.addDataSink(dataTable.makeDataSink(),new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        DisplayTable THING = new DisplayTable(dataTable);
        THING.setTransposed(false);
        THING.setShowingRowLabels(false);
        THING.setPrecision(7);

        DisplayPlot compositionPlot = new DisplayPlot();
        accumulator.addDataSink(compositionPlot.getDataSet().makeDataSink(),new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        compositionPlot.setDoLegend(false);

        getController().getReinitButton().setPostAction(getDisplayPhasePaintAction(sim.phase));

        // Things to Do while simulation is on (It is adding the DataPumps, which run off the meters)
        sim.integratorHard1.addIntervalAction(tPump);
        // Setting up how often it operates. 
        sim.integratorHard1.setActionInterval(tPump, 100);

		DeviceThermoSelector tSelect = setup(sim);
        tSelect.setIntegrator(sim.integratorHard1);
        ColorSchemeByType colorScheme = (ColorSchemeByType)getDisplayPhase(sim.phase).getColorScheme();
        colorScheme.setColor(sim.speciesA.getMoleculeType(), java.awt.Color.red);
        colorScheme.setColor(sim.speciesB.getMoleculeType(), java.awt.Color.black);

        int ms = 20;
        AEditor.nSlider.getSlider().setMajorTickSpacing(ms);
        BEditor.nSlider.getSlider().setMajorTickSpacing(ms);
        AEditor.nSlider.getSlider().setMinorTickSpacing(2);
        BEditor.nSlider.getSlider().setMinorTickSpacing(2);
        AEditor.nSlider.getSlider().setLabelTable(AEditor.nSlider.getSlider().createStandardLabels(ms));
        BEditor.nSlider.getSlider().setLabelTable(AEditor.nSlider.getSlider().createStandardLabels(ms));
        sizeSlider.setLabel("Atom size");
        sizeSlider.setPrecision(2);
        sizeSlider.setMinimum(0.0);
        sizeSlider.setMaximum(3.0);
        sizeSlider.setNMajor(3);
        sizeSlider.setValue(3.0);
        sizeSlider.setShowValues(true);
        sizeSlider.setEditValues(true);
        sizeModifier.setDisplay(getDisplayPhase(sim.phase));

        tBox.setUnit(Kelvin.UNIT);
        tBox.setLabel("Measured value");
        tBox.setLabelPosition(CompassDirection.NORTH);

        compositionPlot.setLabel("Composition");

        getPanel().tabbedPane.add(getDisplayPhase(sim.phase).getLabel(), getDisplayPhase(sim.phase).graphic());
        getPanel().tabbedPane.add(compositionPlot.getLabel(), compositionPlot.graphic());
        getPanel().tabbedPane.add("Averages", THING.graphic());

        //panel for the temperature control/display 
        temperaturePanel.add(tSelect.graphic(null), horizGBC);
        temperaturePanel.add(tBox.graphic(null), horizGBC);

        speciesEditors.add(AEditor);
        speciesEditors.add(BEditor);
        speciesEditors.setBorder(new TitledBorder("Species Adjustment"));
        epsilonSliders.add(AASlider.graphic(null));
        epsilonSliders.add(ABSlider.graphic(null));
        epsilonSliders.add(BBSlider.graphic(null));
        lambdaSliders.add(AAWellSlider.graphic(null));
        lambdaSliders.add(ABWellSlider.graphic(null));
        lambdaSliders.add(BBWellSlider.graphic(null));
        sizeSliders.add(sizeSlider.graphic(null));

        final JTabbedPane sliderPanel = new JTabbedPane();
        sliderPanel.setBorder(new TitledBorder("Potential Adjustment"));
        sliderPanel.add("Well depth (K)", epsilonSliders);
        sliderPanel.add("Core size (%)", lambdaSliders);
        sliderPanel.add("Monomer size (\u00C5)", sizeSliders);
        sliderPanel.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent event) {
                sliderPanel.invalidate();
                sliderPanel.validate();
            }
        });
        sliderPanel.addKeyListener(new java.awt.event.KeyListener() {
            public void keyPressed(java.awt.event.KeyEvent e) {}

            public void keyTyped(java.awt.event.KeyEvent e) {
                etomica.units.Prefix prefix = etomica.units.Prefix.keySelect(e.getKeyChar());
                if (prefix == null || prefix.value() > 1000.0
                        || prefix.value() < 1.0)
                    return;
                etomica.units.Unit unit = new PrefixedUnit (prefix, etomica.units.Kelvin.UNIT);
                AASlider.setUnit(unit);
                ABSlider.setUnit(unit);
                BBSlider.setUnit(unit);
                sliderPanel.setTitleAt(0, "Well depth (" + unit.symbol() + ")");
                AASlider.setLabel("RR epsilon");
                ABSlider.setLabel("RB epsilon");
                BBSlider.setLabel("BB epsilon");
                AASlider.getModifier().setValue(
                        unit.toSim(AASlider.getSlider().getValue()));
                ABSlider.getModifier().setValue(
                        unit.toSim(ABSlider.getSlider().getValue()));
                BBSlider.getModifier().setValue(
                        unit.toSim(BBSlider.getSlider().getValue()));
                sliderPanel.repaint();
            }

            public void keyReleased(java.awt.event.KeyEvent e) {}
		});

        //panel for all the controls
        getPanel().controlPanel.add(temperaturePanel, vertGBC);
        getPanel().controlPanel.add(sliderPanel, vertGBC);
        getPanel().controlPanel.add(speciesEditors, vertGBC);

        //set the number of significant figures displayed on the table.
        javax.swing.table.DefaultTableCellRenderer numberRenderer = new javax.swing.table.DefaultTableCellRenderer() {
            java.text.NumberFormat formatter;
            {
                formatter = java.text.NumberFormat.getInstance();
                formatter.setMaximumFractionDigits(6);
            }

            public void setValue(Object value) {
                setText((value == null) ? "" : formatter.format(value));
            }
        };

        numberRenderer.setHorizontalAlignment(SwingConstants.RIGHT);
    }

    public DeviceThermoSelector setup(ChainEquilibriumSim sim){
        DeviceThermoSelector tSelect = new DeviceThermoSelector(sim.controller1, Kelvin.UNIT, true);
        tSelect.setTemperatures(new double[] { 50., 100., 150., 200., 300.,500., 700., 1000., 1200. });
        tSelect.setUnit(Kelvin.UNIT);
        tSelect.setSelected(3); 	//sets 300K as selected temperature
        tSelect.getLabel().setText("Set value");
        return tSelect;
    }
    public DeviceSlider sliders(ChainEquilibriumSim sim, int eMin, int eMax, String s, P2SquareWellBonded p){

        DeviceSlider AASlider = new DeviceSlider(sim.getController(), new WellDepthModifier(p));
        AASlider.setUnit((Kelvin.UNIT));
        AASlider.setShowBorder(true);
        AASlider.setLabel(s);
        AASlider.setMinimum(eMin);
        AASlider.setMaximum(eMax);
        AASlider.setNMajor(4);
        AASlider.getSlider().setSnapToTicks(true);

        return AASlider;
    }
    public DeviceSlider sliders(ChainEquilibriumSim sim, int eMin, int eMax, String s, int majorSpacing, int minorSpacing,
            P2SquareWellBonded p){
	
        DeviceSlider AASlider = new DeviceSlider(sim.getController(), new WellModifier(p));
        AASlider.setUnit((Kelvin.UNIT));
        AASlider.setShowBorder(true);
        AASlider.setLabel(s);
        AASlider.setMinimum(eMin);
        AASlider.setMaximum(eMax);
        AASlider.setNMajor(5);
        AASlider.getSlider().setSnapToTicks(true);
        AASlider.getSlider().setValue(50);
        AASlider.getSlider().setMajorTickSpacing(majorSpacing);
        AASlider.getSlider().setMinorTickSpacing(minorSpacing);
        AASlider.getSlider().setLabelTable(AASlider.getSlider().createStandardLabels(majorSpacing));
        AASlider.getSlider().setLabelTable(AASlider.getSlider().createStandardLabels(majorSpacing));

        return 	AASlider;
    }

    public static void main(String[] args) {
        ChainEquilibriumSim sim = new ChainEquilibriumSim();
        ChainEquilibriumGraphic graphic = new ChainEquilibriumGraphic(sim);
        SimulationGraphic.makeAndDisplayFrame(graphic.getPanel(), APP_NAME);
    }

    public static class Applet extends javax.swing.JApplet {

        public void init() {
			getRootPane().putClientProperty("defeatSystemEventQueueCheck", Boolean.TRUE);
	        ChainEquilibriumSim sim = new ChainEquilibriumSim();
			getContentPane().add(new ChainEquilibriumGraphic(sim).getPanel());
        }
    }
}