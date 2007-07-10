package etomica.modules.reactionequilibrium;
	
import java.awt.GridBagConstraints;

import javax.swing.JPanel;

import etomica.action.Action;
import etomica.action.SimulationRestart;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.DataSourceCountTime;
import etomica.data.DataSplitter;
import etomica.data.DataTag;
import etomica.data.types.DataTable;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceThermoSelector;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTable;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.modifier.Modifier;
import etomica.potential.P2SquareWell;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Dimension;
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
public class ReactionEquilibriumGraphic extends SimulationGraphic {

	private static final String APP_NAME = "Reaction Equilibrium";
	private static final int REPAINT_INTERVAL = 1000;

	public ReactionEquilibriumGraphic(ReactionEquilibrium sim) {

		super(sim, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);

		GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        getDisplayBox(sim.box).setPixelUnit(new Pixel(10));

        Configuration config = new ConfigurationLattice(new LatticeOrthorhombicHexagonal());
        config.initializeCoordinates(sim.box);

		DeviceThermoSelector tSelect = new DeviceThermoSelector(sim, sim.integratorHard1);
        sim.integratorHard1.addIntervalAction(this.getDisplayBoxPaintAction(sim.box));
        
		tSelect.setTemperatures(new double[] { 50., 100., 300., 600., 1000.,
				1200., 1600., 2000., 2500. });
		tSelect.setUnit(Kelvin.UNIT);
		tSelect.setSelected(3); //sets 300K as selected temperature
		tSelect.getLabel().setText("Set value");
        ((ColorSchemeByType)getDisplayBox(sim.box).getColorScheme()).setColor(sim.speciesA.getMoleculeType(), java.awt.Color.red);
        ((ColorSchemeByType)getDisplayBox(sim.box).getColorScheme()).setColor(sim.speciesB.getMoleculeType(), java.awt.Color.black);

		//	adjustment of species properties
		MySpeciesEditor AEditor = new MySpeciesEditor(sim, sim.box, sim.speciesA, "Red");
		MySpeciesEditor BEditor = new MySpeciesEditor(sim, sim.box, sim.speciesB, "Black");
		int ms = 10;
		AEditor.nSlider.getSlider().setMajorTickSpacing(ms);
		BEditor.nSlider.getSlider().setMajorTickSpacing(ms);
		AEditor.nSlider.getSlider().setMinorTickSpacing(1);
		BEditor.nSlider.getSlider().setMinorTickSpacing(1);
		AEditor.nSlider.getSlider().setLabelTable(
				AEditor.nSlider.getSlider().createStandardLabels(ms));
		BEditor.nSlider.getSlider().setLabelTable(
				AEditor.nSlider.getSlider().createStandardLabels(ms));

		//sliders to adjust potentials well depth
		int eMin = 0;
		int eMax = 1000;
		//    final DeviceSlider AASlider = new DeviceSlider(AAbonded, "epsilon");
		final DeviceSlider AASlider = new DeviceSlider(sim.getController(), new WellDepthModifier(
				sim.AAbonded));
		AASlider.setUnit((Kelvin.UNIT));
		AASlider.setShowBorder(true);
		AASlider.setLabel("RR epsilon");
		AASlider.setMinimum(eMin);
		AASlider.setMaximum(eMax);
		AASlider.setNMajor(5);
		AASlider.getSlider().setSnapToTicks(true);
		//    final DeviceSlider ABSlider = new DeviceSlider(ABbonded, "epsilon");
		final DeviceSlider ABSlider = new DeviceSlider(sim.getController(), new WellDepthModifier(
				sim.ABbonded));
		ABSlider.setUnit((Kelvin.UNIT));
		ABSlider.setShowBorder(true);
		ABSlider.setLabel("RB epsilon");
		ABSlider.setMinimum(eMin);
		ABSlider.setMaximum(eMax);
		ABSlider.setNMajor(5);
		ABSlider.getSlider().setSnapToTicks(true);
		//    final DeviceSlider BBSlider = new DeviceSlider(BBbonded, "epsilon");
		final DeviceSlider BBSlider = new DeviceSlider(sim.getController(), new WellDepthModifier(
				sim.BBbonded));
		BBSlider.setUnit((Kelvin.UNIT));
		BBSlider.setShowBorder(true);
		BBSlider.setLabel("BB epsilon");
		BBSlider.setMinimum(eMin);
		BBSlider.setMaximum(eMax);
		BBSlider.setNMajor(5);
		BBSlider.getSlider().setSnapToTicks(true);

		DiameterModifier sizeModifier = new DiameterModifier(sim.AAbonded,
				sim.ABbonded, sim.BBbonded, sim.speciesA, sim.speciesB);
		DeviceSlider sizeSlider = new DeviceSlider(sim.getController(), sizeModifier);
		//       sizeSlider.setShowBorder(true); //doesn't look right with just one
		// slider in
		// pane
		sizeSlider.setLabel("Atom size");
		sizeSlider.setPrecision(2);
		sizeSlider.setMinimum(0.0);
		sizeSlider.setMaximum(3.0);
		sizeSlider.setNMajor(3);
		sizeSlider.setValue(3.0);
		sizeSlider.setShowValues(true);
		sizeSlider.setEditValues(true);

		//sliders to adjust widths of wells
		eMin = 5;
		eMax = 95;
		int majorSpacing = 15;
		int minorSpacing = 5;
		DeviceSlider AAWellSlider = new DeviceSlider(sim.getController(),new WellModifier(sim.AAbonded));
		DeviceSlider ABWellSlider = new DeviceSlider(sim.getController(),new WellModifier(sim.ABbonded));
		DeviceSlider BBWellSlider = new DeviceSlider(sim.getController(),new WellModifier(sim.BBbonded));
		AAWellSlider.setShowBorder(true);
		ABWellSlider.setShowBorder(true);
		BBWellSlider.setShowBorder(true);
		AAWellSlider.setLabel("RR core");
		ABWellSlider.setLabel("RB core");
		BBWellSlider.setLabel("BB core");
		AAWellSlider.setMinimum(eMin);
		AAWellSlider.setMaximum(eMax);
		ABWellSlider.setMinimum(eMin);
		ABWellSlider.setMaximum(eMax);
		BBWellSlider.setMinimum(eMin);
		BBWellSlider.setMaximum(eMax);
		AAWellSlider.getSlider().setSnapToTicks(true);
		ABWellSlider.getSlider().setSnapToTicks(true);
		BBWellSlider.getSlider().setSnapToTicks(true);
		AAWellSlider.getSlider().setValue(50);
		ABWellSlider.getSlider().setValue(50);
		BBWellSlider.getSlider().setValue(50);
		AAWellSlider.getSlider().setMajorTickSpacing(majorSpacing);
		AAWellSlider.getSlider().setMinorTickSpacing(minorSpacing);
		ABWellSlider.getSlider().setMajorTickSpacing(majorSpacing);
		ABWellSlider.getSlider().setMinorTickSpacing(minorSpacing);
		BBWellSlider.getSlider().setMajorTickSpacing(majorSpacing);
		BBWellSlider.getSlider().setMinorTickSpacing(minorSpacing);
		AAWellSlider.getSlider().setLabelTable(
				AAWellSlider.getSlider().createStandardLabels(majorSpacing));
		AAWellSlider.getSlider().setLabelTable(
				AAWellSlider.getSlider().createStandardLabels(majorSpacing));
		ABWellSlider.getSlider().setLabelTable(
				ABWellSlider.getSlider().createStandardLabels(majorSpacing));
		ABWellSlider.getSlider().setLabelTable(
				ABWellSlider.getSlider().createStandardLabels(majorSpacing));
		BBWellSlider.getSlider().setLabelTable(
				BBWellSlider.getSlider().createStandardLabels(majorSpacing));
		BBWellSlider.getSlider().setLabelTable(
				BBWellSlider.getSlider().createStandardLabels(majorSpacing));

		//so that display is updated when slider changes atom sizes
		sizeModifier.setDisplay(getDisplayBox(sim.box));
        
		DisplayTextBox tBox = new DisplayTextBox(sim.thermometer.getDataInfo());
		DataPump tPump = new DataPump (sim.thermometer, tBox);
        sim.integratorHard1.addIntervalAction(tPump);
        sim.integratorHard1.setActionInterval(tPump, 100);
        tBox.setUnit(Kelvin.UNIT);
		tBox.setLabel("Measured value");
		tBox.setLabelPosition(CompassDirection.NORTH);

// NOTE : THE FOLLOWING IMPLEMENTATION IS CAUSING THE GRAPHIC
// TO BE UPDATED.  NEED TO REMOVE THE UPDATE OF THE GRAPHIC
// FROM THIS IMPLEMENTATION OF AVERAGE DISPLAY.
		//display of averages
        DataFork dimerFork = new DataFork();
		DataPump dimerPump = new DataPump (sim.meterDimerFraction, dimerFork);
        sim.integratorHard1.addIntervalAction(dimerPump);
        sim.integratorHard1.setActionInterval(dimerPump, 100);
        AccumulatorAverage dimerfractionaccum = new AccumulatorAverage();
        dimerfractionaccum.setPushInterval(10);
        dimerFork.addDataSink(dimerfractionaccum);
		DisplayTable table = new DisplayTable();
		dimerfractionaccum.setDataSink(table.getDataTable().makeDataSink());

        DataSplitter splitter = new DataSplitter();
        
        dimerFork.addDataSink(splitter);
         
		//display for history of mole fractions
        DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integratorHard1);
        DisplayPlot plot = new DisplayPlot();
        plot.setLabel("Composition");
        plot.setDoLegend(true);
        int nData = sim.meterDimerFraction.getDataInfo().getLength();
        DataTable.DataInfoTable dimerInfo = (DataTable.DataInfoTable)sim.meterDimerFraction.getDataInfo();

         for (int i=0; i<nData; i++) {
            AccumulatorHistory dimerfractionhistory = new AccumulatorHistory();
            dimerfractionhistory.setTimeDataSource(timeCounter);
            
    		splitter.setDataSink(i, dimerfractionhistory);
    		dimerfractionhistory.addDataSink (plot.getDataSet().makeDataSink());
    		plot.setLegend(new DataTag[]{dimerfractionhistory.getTag()}, dimerInfo.getRowHeader(i));
        }

		//************* Lay out components ****************//

		//panel for the temperature control/display
		JPanel temperaturePanel = new JPanel(new java.awt.GridBagLayout());
		temperaturePanel.setBorder(new javax.swing.border.TitledBorder(
				"Temperature (K)"));
		temperaturePanel.add(tSelect.graphic(null), SimulationPanel.getHorizGBC());
		temperaturePanel.add(tBox.graphic(null), SimulationPanel.getHorizGBC());

		//panel for the species editors
		JPanel speciesEditors = new JPanel(new java.awt.GridLayout(0, 1));
		speciesEditors.add(AEditor);
		speciesEditors.add(BEditor);
		speciesEditors.setBorder(new javax.swing.border.TitledBorder(
				"Species Adjustment"));

		//panel of well-depth sliders
		JPanel epsilonSliders = new JPanel(new java.awt.GridLayout(0, 1));
		epsilonSliders.add(AASlider.graphic(null));
		epsilonSliders.add(ABSlider.graphic(null));
		epsilonSliders.add(BBSlider.graphic(null));

		//panel of well-width sliders
		JPanel lambdaSliders = new JPanel(new java.awt.GridLayout(0, 1));
		lambdaSliders.add(AAWellSlider.graphic(null));
		lambdaSliders.add(ABWellSlider.graphic(null));
		lambdaSliders.add(BBWellSlider.graphic(null));

		//panel for size slider
		JPanel sizeSliders = new JPanel(new java.awt.GridLayout(0, 1));
		sizeSliders.add(sizeSlider.graphic(null));

		//tabbed pane for both sets of sliders
		final javax.swing.JTabbedPane sliderPanel = new javax.swing.JTabbedPane();
		sliderPanel.setBorder(new javax.swing.border.TitledBorder(
				"Potential Adjustment"));
		sliderPanel.add("Well depth (K)", epsilonSliders);
		sliderPanel.add("Core size (%)", lambdaSliders);
		sliderPanel.add("Atom size (\u00C5)", sizeSliders);
		sliderPanel.addChangeListener(new javax.swing.event.ChangeListener() {
			public void stateChanged(javax.swing.event.ChangeEvent event) {
				sliderPanel.invalidate();
				sliderPanel.validate();
			}
		});
		sliderPanel.addKeyListener(new java.awt.event.KeyListener() {
			public void keyPressed(java.awt.event.KeyEvent e) {
			}

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

			public void keyReleased(java.awt.event.KeyEvent e) {
			}
		});

		//top panel for control, temperature, potential adjustment
		getPanel().controlPanel.add(temperaturePanel, vertGBC);
		getPanel().controlPanel.add(sliderPanel, vertGBC);
		getPanel().controlPanel.add(speciesEditors, vertGBC);
		add(plot);
		add(table);


        getController().getReinitButton().setPostAction(getDisplayBoxPaintAction(sim.box));

		//***************set all the colors******************
		/*
		 * java.awt.Color background = Constants.KHAKI.brighter().brighter();
		 * java.awt.Color contrast = Constants.DARK_RED; java.awt.Color
		 * sliderColor = Constants.TAN; // java.awt.Color tabColor =
		 * java.awt.Color.black; java.awt.Color tabColor = Constants.DARK_RED;
		 * java.awt.Color tabTextColor = background; java.awt.Color buttonColor =
		 * tabColor; java.awt.Color buttonTextColor = tabTextColor;
		 * java.awt.Color panelColor = Constants.TAN;
		 * panel.setBackground(background);
		 * devicePanel.setBackground(background);
		 * displayPanel.setBackground(tabColor);
		 * displayPanel.setForeground(tabTextColor);
		 * displayBoxPanel.setBackground(background);
		 * 
		 * sliderPanel.setBackground(tabColor);
		 * sliderPanel.setForeground(tabTextColor);
		 * sliderPanel.setOpaque(false); epsilonSliders.setOpaque(false);
		 * lambdaSliders.setOpaque(false);
		 * AAWellSlider.graphic(null).setBackground(sliderColor);
		 * ABWellSlider.graphic(null).setBackground(sliderColor);
		 * BBWellSlider.graphic(null).setBackground(sliderColor);
		 * AASlider.graphic(null).setBackground(sliderColor);
		 * ABSlider.graphic(null).setBackground(sliderColor);
		 * BBSlider.graphic(null).setBackground(sliderColor);
		 * 
		 * controlPanel.setBackground(background);
		 * 
		 * startPanel.setBackground(contrast); //border color
		 * startPanel.setOpaque(false);
		 * 
		 * temperaturePanel.setBackground(contrast); //border color
		 * temperaturePanel.setOpaque(false);
		 * tBox.graphic(null).setBackground(background);
		 * tSelect.graphic(null).setBackground(background);
		 * 
		 * displayBox1.graphic(null).setBackground(panelColor);
		 * plot.graphic(null).setBackground(panelColor);
		 * plot.getPlot().setBackground(panelColor); //doesn't have intended
		 * effect table.graphic(null).setBackground(panelColor);
		 * speciesA.setColor(Constants.BRIGHT_RED);
		 * speciesB.setColor(java.awt.Color.black);
		 * 
		 * controller1.getButton().graphic(null).setBackground(buttonColor);
		 * controller1.getButton().graphic(null).setForeground(buttonTextColor);
		 * restart.graphic(null).setBackground(buttonColor);
		 * restart.graphic(null).setForeground(buttonTextColor);
		 * resetAverages.graphic(null).setBackground(buttonColor);
		 * resetAverages.graphic(null).setForeground(buttonTextColor); //
		 */

		//set the number of significant figures displayed on the table.
		javax.swing.table.DefaultTableCellRenderer numberRenderer = new javax.swing.table.DefaultTableCellRenderer() {
			java.text.NumberFormat formatter;
			//initializer for inner class
			{
				formatter = java.text.NumberFormat.getInstance();
				formatter.setMaximumFractionDigits(6);//this is the only part
				// that
				// differs from the default
			}

			public void setValue(Object value) {
				setText((value == null) ? "" : formatter.format(value));
			}
		};
		numberRenderer.setHorizontalAlignment(javax.swing.JLabel.RIGHT);
		//table.table.setDefaultRenderer(Number.class, numberRenderer);

		initializing = false;
		//panel.removeAll();
		//panel.add(panel);
	}

	public static void main(String[] args) {

		ReactionEquilibrium sim = new ReactionEquilibrium();
		ReactionEquilibriumGraphic graphic = new ReactionEquilibriumGraphic(sim);
		SimulationGraphic.makeAndDisplayFrame(graphic.getPanel(), APP_NAME);
	}

	//=================================================================
	//panel containing species-editing devices

	class MySpeciesEditor extends javax.swing.JPanel {

		//	public DeviceSlider nSlider;
		public DeviceNSelector nSlider;

		public Species species;

		public final javax.swing.JTextField mass = new javax.swing.JTextField(
				"40");

		//    public java.awt.TextField mass = new java.awt.TextField("40");

        public MySpeciesEditor(final ReactionEquilibrium sim, Box box, Species s, String label) {
            super();
            species = s;
            nSlider = new DeviceNSelector(sim.getController());
            nSlider.setResetAction(new SimulationRestart(sim));
            nSlider.setSpecies(species);
            nSlider.setBox(box);
            //nSlider.setDisplayBox(DisplayBox1);
            nSlider.setMinimum(0);
            nSlider.setMaximum(40);
            nSlider.setPostAction(new Action() {
                public void actionPerformed() {
                    AtomAgentManager agentManager = sim.getAgentManager();
                    AtomIteratorLeafAtoms iter = new AtomIteratorLeafAtoms(sim.box);
                    iter.reset();
                    for (IAtom a = iter.nextAtom(); a != null; a = iter.nextAtom()) {
                        //                      System.out.println(iter.peek().toString());
                        agentManager.setAgent(a, null);
                    }
                    try {
                    	sim.integratorHard1.reset();
                    } catch(ConfigurationOverlapException e) {}
                    getDisplayBox(sim.box).repaint();
                }
           });

			//listener for changes to mass textbox
			java.awt.event.ActionListener myListener = new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent event) {
					if (initializing)
						return;
					int value;
					try {
						value = Integer.parseInt(mass.getText());
					} catch (NumberFormatException ex) {
						return;
					}
					if (value < 1)
						value = 1;
					if (value > 1000000)
						value = 1000000;
					final double newMass = value;
					mass.setText(Integer.toString(value));
					((ElementSimple)((AtomTypeLeaf)species.getMoleculeType()).getElement()).setMass(newMass);
                     try {
                         sim.integratorHard1.reset();
                     } catch(ConfigurationOverlapException e) {}
				}
			};
			mass.addActionListener(myListener);
			mass.setBorder(new javax.swing.border.TitledBorder("Mass"));
			mass.setColumns(6);
			mass.setOpaque(false);
			setLayout(new java.awt.FlowLayout());
			add(nSlider.graphic(null));
			add(mass);
			setBorder(new javax.swing.border.TitledBorder(label));
		}
	} //end of MySpeciesEditor

	public static class Applet extends javax.swing.JApplet {

		public void init() {
			getRootPane().putClientProperty("defeatSystemEventQueueCheck",
					Boolean.TRUE);
			ReactionEquilibrium sim = new ReactionEquilibrium();
			ReactionEquilibriumGraphic graphic = new ReactionEquilibriumGraphic(sim);
			SimulationGraphic.makeAndDisplayFrame(graphic.getPanel(), APP_NAME);
		}
	}//end of Applet

	//---------------------------------------------

	class DiameterModifier implements Modifier {
		P2SquareWellBonded potentialRR, potentialRB, potentialBB;

		SpeciesSpheresMono speciesR, speciesB;

		DisplayBox display;

		DiameterModifier(P2SquareWellBonded potentialRR,
				P2SquareWellBonded potentialRB, P2SquareWellBonded potentialBB,
				SpeciesSpheresMono speciesR, SpeciesSpheresMono speciesB) {
			this.potentialRR = potentialRR;
			this.potentialRB = potentialRB;
			this.potentialBB = potentialBB;
			this.speciesR = speciesR;
			this.speciesB = speciesB;
		}

		public Dimension getDimension() {
			return etomica.units.Length.DIMENSION;
		}

		public void setValue(double d) {
			if (d == 0.0)
				d = 0.01;
			double changeFraction = d
					/ (potentialRR.getCoreDiameter() * potentialRR.getLambda());
			double newCoreDiameter = changeFraction
					* potentialRR.getCoreDiameter();
			potentialRR.setCoreDiameter(newCoreDiameter);
			potentialRB.setCoreDiameter(newCoreDiameter);
			potentialBB.setCoreDiameter(newCoreDiameter);
			((AtomTypeSphere)speciesR.getMoleculeType()).setDiameter(d);
			((AtomTypeSphere)speciesB.getMoleculeType()).setDiameter(d);
			if (display != null)
				display.repaint();
		}

		public double getValue() {
			return ((AtomTypeSphere)speciesR.getMoleculeType()).getDiameter();
		}

		public void setDisplay(DisplayBox display) {
			this.display = display;
		}

		public String getLabel() {
			return "Diameter";
		}
	}

	//---------------------------------------------

	class WellDepthModifier implements Modifier {

		P2SquareWell potential;

		WellDepthModifier(P2SquareWell pot) {
			potential = pot;
		}

		public String getLabel() {
			return "WellDepth";
		}

		public Dimension getDimension() {
			return etomica.units.Energy.DIMENSION;
		}

		public void setValue(double d) {
			potential.setEpsilon(d);
		}

		public double getValue() {
			return potential.getEpsilon();
		}
	}

	class WellModifier implements Modifier {

		double currentValue;

		double fullDiameter;

		P2SquareWell potential;

		WellModifier(P2SquareWell pot) {
			potential = pot;
			fullDiameter = potential.getCoreDiameter() * potential.getLambda();
		}

		public String getLabel() {
			return "WellMod";
		}

		public Dimension getDimension() {
			return etomica.units.Null.DIMENSION;
		}

		public void setValue(double d) {
			if (initializing)
				return;
			currentValue = d;
			double x = 0.01 * currentValue;
			fullDiameter = potential.getCoreDiameter() * potential.getLambda();
			potential.setCoreDiameter(x * fullDiameter);
			potential.setLambda(1.0 / x);
		}

		public double getValue() {
			return currentValue;
		}
	}//end of WellModulator
	boolean initializing;

}