package etomica.modules.chainequilibrium;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
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
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTable;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.modifier.Modifier;
import etomica.units.Dimension;
import etomica.units.Kelvin;
import etomica.units.Null;
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
	private static final int REPAINT_INTERVAL = 1;

    protected ChainEquilibriumSim sim;

    public ChainEquilibriumGraphic(ChainEquilibriumSim simulation) {

		super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);
        this.sim = simulation;
        
        ArrayList dataStreamPumps = getController().getDataStreamPumps();
        
        getDisplayBox(sim.box).setPixelUnit(new Pixel(10));

        GridBagConstraints horizGBC = SimulationPanel.getHorizGBC();
        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        // ********* Data Declaration Section *******	
        int eMin = 0, eMax = 1000,  majorSpacing = 15, minorSpacing = 5;

        // **** Stuff that Modifies the Simulation
        
        final DeviceSlider delaySlider = new DeviceSlider(sim.controller1, new Modifier() {
            public double getValue() {return Math.pow(sim.activityIntegrate.getSleepPeriod(),1./2.);}
            public void setValue(double d) {sim.activityIntegrate.setSleepPeriod((int)(d*d));}
            public Dimension getDimension() {return Null.DIMENSION;}
            public String getLabel() {return "Sleep period";}
        });
        delaySlider.setMinimum(0);
        delaySlider.setMaximum(10);
        delaySlider.setNMajor(5);
        delaySlider.setValue(0);
        delaySlider.setShowMinorTicks(true);
        JPanel delaySliderPanel = new JPanel();
        delaySliderPanel.add(delaySlider.graphic());
        delaySliderPanel.setBorder(new javax.swing.border.TitledBorder("Delay"));
        

        // Sliders on Well depth page
        final DeviceSlider AASlider = sliders(sim, eMin, eMax, "RR epsilon", sim.AAbonded);
        final DeviceSlider ABSlider = sliders(sim, eMin, eMax, "RB epsilon", sim.ABbonded);
        final DeviceSlider BBSlider = sliders(sim, eMin, eMax, "BB epsilon", sim.BBbonded);
        
        // Sliders on Core size page
        DeviceSlider AAWellSlider = sliders(sim, 5, 95, "RR core", majorSpacing, minorSpacing, sim.AAbonded);
        DeviceSlider ABWellSlider = sliders(sim, 5, 95, "RB core", majorSpacing, minorSpacing, sim.ABbonded);
        DeviceSlider BBWellSlider = sliders(sim, 5, 95, "BB core", majorSpacing, minorSpacing, sim.BBbonded);

        // The Species Editors
        MySpeciesEditor AEditor = new MySpeciesEditor(this, sim.box, sim.speciesA, "Red");
        MySpeciesEditor BEditor = new MySpeciesEditor(this, sim.box, sim.speciesB, "Black");
		
        // the Atom Diameter Modifer
        DiameterModifier sizeModifier = new DiameterModifier(sim.AAbonded,sim.ABbonded, sim.BBbonded, sim.speciesA, sim.speciesB);

        DeviceSlider sizeSlider = new DeviceSlider(sim.getController(), sizeModifier);
		
        DisplayTextBox tBox = new DisplayTextBox();

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
        
        sim.integratorHard1.setTimeStep(0.01);

        DataSinkTable dataTable = new DataSinkTable();
        accumulator.addDataSink(dataTable.makeDataSink(),new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        DisplayTable THING = new DisplayTable(dataTable);
        THING.setTransposed(false);
        THING.setShowingRowLabels(false);
        THING.setPrecision(7);

        DisplayPlot compositionPlot = new DisplayPlot();
        accumulator.addDataSink(compositionPlot.getDataSet().makeDataSink(),new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        compositionPlot.setDoLegend(false);

        getController().getReinitButton().setPostAction(getPaintAction(sim.box));

        // Things to Do while simulation is on (It is adding the DataPumps, which run off the meters)
        sim.integratorHard1.addIntervalAction(tPump);
        // Setting up how often it operates. 
        sim.integratorHard1.setActionInterval(tPump, 100);

        DeviceThermoSelector tSelect = new DeviceThermoSelector(sim.controller1, Kelvin.UNIT, true);
        tSelect.setTemperatures(new double[] { 50., 100., 150., 200., 300.,500., 700., 1000., 1200. });
        tSelect.setUnit(Kelvin.UNIT);
        tSelect.setIntegrator(sim.integratorHard1);
        tSelect.setSelected(3);     //sets 150K as selected temperature (0th selection is "Adiabatic")
        tSelect.getLabel().setText("Set value");
        
        ColorSchemeByType colorScheme = (ColorSchemeByType)getDisplayBox(sim.box).getColorScheme();
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
        sizeModifier.setDisplay(getDisplayBox(sim.box));

        tBox.setUnit(Kelvin.UNIT);
        tBox.setLabel("Measured value");
        tBox.setLabelPosition(CompassDirection.NORTH);

        compositionPlot.setLabel("Composition");

        getPanel().tabbedPane.add(getDisplayBox(sim.box).getLabel(), getDisplayBox(sim.box).graphic());
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
//        ((JPanel)getPanel().tabbedPane.getComponentAt(0)).add(delaySliderPanel);
//        JPanel panel1 = new JPanel(new GridBagLayout());
//        panel1.add(temperaturePanel,vertGBC);
//        panel1.add(sliderPanel,vertGBC);
//        getPanel().controlPanel.add(panel1, horizGBC);
        getPanel().controlPanel.add(delaySliderPanel, vertGBC);
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