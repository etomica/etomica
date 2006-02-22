package etomica.modules.chainequilibrium;

import javax.swing.JPanel;
import javax.swing.SwingConstants;

import etomica.action.Action;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.DataSinkTable;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceThermoSelector;
import etomica.graphics.DeviceTrioControllerButton;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTable;
import etomica.integrator.IntervalActionAdapter;
import etomica.units.Kelvin;
import etomica.units.PrefixedUnit;
import etomica.util.Constants.CompassDirection;

/**
 * @author William Scharmach
 * 
 * TODO To change the template for this generated type comment go to Window -
 * Preferences - Java - Code Style - Code Templates
 */
public class ReactionEquilibriumGraphic {
    public JPanel panel = new JPanel();
    boolean initializing;
    public DisplayPhase displayPhase1;
    public ReactionEquilibrium simulation;

    public ReactionEquilibriumGraphic(ReactionEquilibrium sim) {
		
        simulation = sim;
        sim.register(sim.integratorHard1);
        // ********* Data Declaration Section *******	
        initializing = true;
        int eMin = 0, eMax = 1000,  majorSpacing = 15, minorSpacing = 5;

        // **** Stuff that Modifies the Simulation

        // All the SLIDERS
        final DeviceSlider AASlider = sliders(sim, eMin, eMax, "RR epsilon", sim.AAbonded);
        final DeviceSlider ABSlider = sliders(sim, eMin, eMax, "RB epsilon", sim.ABbonded);
        final DeviceSlider BBSlider = sliders(sim, eMin, eMax, "BB epsilon", sim.BBbonded);
        DeviceSlider AAWellSlider = sliders(sim, 5, 95, "RR core", majorSpacing, minorSpacing, sim.AAbonded);
        DeviceSlider ABWellSlider = sliders(sim, 5, 95, "RB core", majorSpacing, minorSpacing, sim.ABbonded);
        DeviceSlider BBWellSlider = sliders(sim, 5, 95, "BB core", majorSpacing, minorSpacing, sim.BBbonded);

        // The Species Editors
        MySpeciesEditor AEditor = new MySpeciesEditor(this, sim.speciesA.getAgent(sim.phase1), "Red");
        MySpeciesEditor BEditor = new MySpeciesEditor(this, sim.speciesB.getAgent(sim.phase1), "Black");
		
        // the Atom Diameter Modifer
        DiameterModifier sizeModifier = new DiameterModifier(sim.AAbonded,sim.ABbonded, sim.BBbonded, sim.speciesA, sim.speciesB);
		
        // All the Devices
        DeviceTrioControllerButton control = new DeviceTrioControllerButton(sim);
        DeviceSlider sizeSlider = new DeviceSlider(sim.getController(), sizeModifier);
		
        // *** all the display stuff
		
        // DISPLAYPHASE is all the stuff you see on the Right, the pic, the graphs
        // the FINAL line makes the Pairent and everything below is attached to DisplayPhase1
        displayPhase1 = new DisplayPhase(sim.phase1,sim.getDefaults().pixelUnit);

        DisplayBox tBox = new DisplayBox();


        // DISPLAYPANEL is all the stuff you see to the Left, the controls, the temperature box
        // the FINAL line makes the Pairent and everything below is attached
        final javax.swing.JTabbedPane displayPanel = new javax.swing.JTabbedPane();
        JPanel startPanel = (JPanel) control.graphic();
        JPanel speciesEditors = new JPanel(new java.awt.GridLayout(0, 1));
        JPanel epsilonSliders = new JPanel(new java.awt.GridLayout(0, 1));
        JPanel controlPanel = new JPanel(new java.awt.GridBagLayout());
        JPanel temperaturePanel = new JPanel(new java.awt.GridBagLayout());
        JPanel sizeSliders = new JPanel(new java.awt.GridLayout(0, 1));
        JPanel lambdaSliders = new JPanel(new java.awt.GridLayout(0, 1));
        JPanel topPanel = new JPanel(new java.awt.GridBagLayout());
        temperaturePanel.setBorder(new javax.swing.border.TitledBorder("Temperature (K)"));
        java.awt.GridBagConstraints gbc1 = new java.awt.GridBagConstraints();

        DataPump tPump = new DataPump (sim.thermometer, tBox);

        AccumulatorAverage accumulator = new AccumulatorAverage(sim);
        DataPump pump = new DataPump(sim.molecularCount,accumulator);
        new IntervalActionAdapter(pump,sim.integratorHard1);

        DataSinkTable dataTable = new DataSinkTable();
        accumulator.addDataSink(dataTable,new AccumulatorAverage.StatType[]{AccumulatorAverage.StatType.AVERAGE});
        DisplayTable THING = new DisplayTable(dataTable);
        THING.setRowLabels(new String[] { "monomer", "dimer", "trimer", "4-mer", "5-mer", "6-mer", "7-10-mer", "11-13-mer", "14-25-mer",">25-mer"});
        THING.setTransposed(false);
        THING.setShowingRowLabels(true);
        THING.setPrecision(7);

        DisplayPlot compositionplot = new DisplayPlot(dataTable);
		
        // Stuff that Happens AS SIMULATION RUNS !
        // Action Preform Method only will work on FINAL types of data which never change mid program
        // It works on Display, which is a FINAL type picture
        // There may be a way to set up a FINAL System.out. Something that would print Data to the screen
        Action RefreshPicture =  new Action() {
            public void actionPerformed() {
               	displayPhase1.repaint();
            }
            public String getLabel() {return "";}
        };

        // Things to Do while simulation is on (It is adding the DataPumps, which run off the meters)
        IntervalActionAdapter tAdapter = new IntervalActionAdapter (tPump, sim.integratorHard1);
        sim.integratorHard1.addListener(new IntervalActionAdapter(RefreshPicture));
        
        // Setting up how often it operates. 
		tAdapter.setActionInterval(100);

		DeviceThermoSelector tSelect = setup(sim);
        tSelect.setIntegrator(sim.integratorHard1);
        ColorSchemeByType colorScheme = (ColorSchemeByType)displayPhase1.getColorScheme();
        colorScheme.setColor(sim.speciesA.getMoleculeType(), java.awt.Color.red);
        colorScheme.setColor(sim.speciesB.getMoleculeType(), java.awt.Color.black);

        int ms = 10;
        AEditor.nSlider.getSlider().setMajorTickSpacing(ms);
        BEditor.nSlider.getSlider().setMajorTickSpacing(ms);
        AEditor.nSlider.getSlider().setMinorTickSpacing(1);
        BEditor.nSlider.getSlider().setMinorTickSpacing(1);
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
        sizeModifier.setDisplay(displayPhase1);

        tBox.setUnit(Kelvin.UNIT);
        tBox.setLabel("Measured value");
        tBox.setLabelPosition(CompassDirection.NORTH);
        //dimerfractionaccum.addDataSink(table.getDataTable().makeColumn(Dimension.FRACTION));

        //dimerfractionhistory.addDataSink (compositionplot.getDataTable().makeColumn(Dimension.FRACTION));
        //tAverageDimer.addDataSink(dimerfractionhistory);

        compositionplot.setLabel("Composition");

        displayPanel.add(displayPhase1.getLabel(), displayPhase1.graphic());
        displayPanel.add(compositionplot.getLabel(), compositionplot.graphic(compositionplot.getPlot()));
        displayPanel.add("Averages", THING.graphic());

        //workaround for JTabbedPane bug in JDK 1.2
        displayPanel.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent event) {
                displayPanel.invalidate();
                displayPanel.validate();
            }
        });

        //panel for the temperature control/display
        // ********* Marker 
        gbc1.gridx = 0;
        gbc1.gridy = 0;
        gbc1.gridwidth = 1;
        temperaturePanel.add(tSelect.graphic(null), gbc1);
        gbc1.gridx = 0;
        gbc1.gridy = 1;
        temperaturePanel.add(tBox.graphic(null));

        speciesEditors.add(AEditor);
        speciesEditors.add(BEditor);
        speciesEditors.setBorder(new javax.swing.border.TitledBorder("Species Adjustment"));
        epsilonSliders.add(AASlider.graphic(null));
        epsilonSliders.add(ABSlider.graphic(null));
        epsilonSliders.add(BBSlider.graphic(null));
        lambdaSliders.add(AAWellSlider.graphic(null));
        lambdaSliders.add(ABWellSlider.graphic(null));
        lambdaSliders.add(BBWellSlider.graphic(null));
        sizeSliders.add(sizeSlider.graphic(null));

        final javax.swing.JTabbedPane sliderPanel = new javax.swing.JTabbedPane();
        sliderPanel.setBorder(new javax.swing.border.TitledBorder(
                "Potential Adjustment"));
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

		//top panel for control, temperature, potential adjustment
        java.awt.GridBagConstraints gbc2 = new java.awt.GridBagConstraints();
        gbc2.gridx = 0;
        gbc2.gridy = 0;
        gbc2.gridheight = 1;
        topPanel.add(startPanel, gbc2);
        gbc2.gridx = 0;
        gbc2.gridy = 1;
        gbc2.gridheight = 1;
        topPanel.add(temperaturePanel, gbc2);
        gbc2.gridx = 1;
        gbc2.gridy = 0;
        gbc2.gridheight = 2;
        topPanel.add(sliderPanel, gbc2);

        //panel for all the controls
        gbc2.gridx = 0;
        gbc2.gridy = java.awt.GridBagConstraints.RELATIVE;
        controlPanel.add(topPanel, gbc2);
        controlPanel.add(speciesEditors, gbc2);

        panel.add(controlPanel, java.awt.BorderLayout.WEST);
        panel.add(displayPanel, java.awt.BorderLayout.EAST);

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
        initializing = false;
    }
    
    public DeviceThermoSelector setup(ReactionEquilibrium sim){
        DeviceThermoSelector tSelect = new DeviceThermoSelector(sim.controller1, Kelvin.UNIT, true);
        tSelect.setTemperatures(new double[] { 50., 100., 150., 200., 300.,500., 700., 1000., 1200. });
        tSelect.setUnit(Kelvin.UNIT);
        tSelect.setSelected(3); 	//sets 300K as selected temperature
        tSelect.getLabel().setText("Set value");
        return tSelect;
    }
    public DeviceSlider sliders(ReactionEquilibrium sim, int eMin, int eMax, String s, P2SquareWellBonded p){

        DeviceSlider AASlider = new DeviceSlider(sim.getController(), new WellDepthModifier(p));
        AASlider.setUnit((Kelvin.UNIT));
        AASlider.setShowBorder(true);
        AASlider.setLabel(s);
        AASlider.setMinimum(eMin);
        AASlider.setMaximum(eMax);
        AASlider.setNMajor(5);
        AASlider.getSlider().setSnapToTicks(true);

        return AASlider;
    }
    public DeviceSlider sliders(ReactionEquilibrium sim, int eMin, int eMax, String s, int majorSpacing, int minorSpacing,
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
        javax.swing.JFrame f = new javax.swing.JFrame(); //create a window
        f.setSize(800, 550);
        
        ReactionEquilibrium sim = new ReactionEquilibrium();
        
        ReactionEquilibriumGraphic graphic = new ReactionEquilibriumGraphic(sim);
        f.getContentPane().add(graphic.panel);
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {
            //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {
                System.exit(0);
            }
        });
    }

    public static class Applet extends javax.swing.JApplet {

        public void init() {
			getRootPane().putClientProperty("defeatSystemEventQueueCheck", Boolean.TRUE);
			getContentPane().add(new ReactionEquilibriumGraphic(new ReactionEquilibrium()).panel);
        }
    }
}