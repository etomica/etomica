package etomica.modules.reactionequilibrium;

import javax.swing.JPanel;

import etomica.action.Action;
import etomica.action.SimulationRestart;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.atom.SpeciesAgent;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.chem.elements.ElementSimple;
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
import etomica.graphics.DeviceTrioControllerButton;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTable;
import etomica.integrator.IntervalActionAdapter;
import etomica.modifier.Modifier;
import etomica.potential.P2SquareWell;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Dimension;
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

	public ReactionEquilibriumGraphic(ReactionEquilibrium sim) {
		sim.register(sim.integratorHard1);
		initializing = true;
		DeviceTrioControllerButton control = new DeviceTrioControllerButton(sim);
		DeviceThermoSelector tSelect = new DeviceThermoSelector(sim, sim.integratorHard1);
		displayPhase1 = new DisplayPhase(sim.phase1,sim.getDefaults().pixelUnit);
        sim.integratorHard1.addListener(new IntervalActionAdapter(
                new Action() {
                    public void actionPerformed() {
                        displayPhase1.repaint();
                    }
                }));
        
		tSelect.setTemperatures(new double[] { 50., 100., 300., 600., 1000.,
				1200., 1600., 2000., 2500. });
		tSelect.setUnit(Kelvin.UNIT);
		tSelect.setSelected(3); //sets 300K as selected temperature
		tSelect.getLabel().setText("Set value");
        ((ColorSchemeByType)displayPhase1.getColorScheme()).setColor(sim.speciesA.getMoleculeType(), java.awt.Color.red);
        ((ColorSchemeByType)displayPhase1.getColorScheme()).setColor(sim.speciesB.getMoleculeType(), java.awt.Color.black);

		//	adjustment of species properties
		MySpeciesEditor AEditor = new MySpeciesEditor(sim, 
				sim.speciesA.getAgent(sim.phase1), "Red");
		MySpeciesEditor BEditor = new MySpeciesEditor(sim, 
				sim.speciesB.getAgent(sim.phase1), "Black");
//		AEditor.nSlider.getSlider().setValue(21);
//		BEditor.nSlider.getSlider().setValue(21);
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
		sizeModifier.setDisplay(displayPhase1);
        
		DisplayBox tBox = new DisplayBox(sim.thermometer.getDataInfo());
		DataPump tPump = new DataPump (sim.thermometer, tBox);
		IntervalActionAdapter tAdapter = new IntervalActionAdapter (tPump, sim.integratorHard1);
		tAdapter.setActionInterval(100);
        tBox.setUnit(Kelvin.UNIT);
		tBox.setLabel("Measured value");
		tBox.setLabelPosition(CompassDirection.NORTH);

		//display of averages
        DataFork dimerFork = new DataFork();
		DataPump dimerPump = new DataPump (sim.meterDimerFraction, dimerFork);
        IntervalActionAdapter dAdapter = new IntervalActionAdapter (dimerPump, sim.integratorHard1);
        dAdapter.setActionInterval(100);
        AccumulatorAverage dimerfractionaccum = new AccumulatorAverage(sim);
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

		//tabbed pane for the big displays
		final javax.swing.JTabbedPane displayPanel = new javax.swing.JTabbedPane();
		displayPanel.add(displayPhase1.getLabel(), displayPhase1.graphic(null));
		displayPanel.add(plot.getLabel(), plot.graphic(null));
		displayPanel.add("Averages", table.graphic(null));
		//workaround for JTabbedPane bug in JDK 1.2
		displayPanel.addChangeListener(new javax.swing.event.ChangeListener() {
			public void stateChanged(javax.swing.event.ChangeEvent event) {
				displayPanel.invalidate();
				displayPanel.validate();
			}
		});

		JPanel startPanel = (JPanel) control.graphic();
		//panel for the temperature control/display
		JPanel temperaturePanel = new JPanel(new java.awt.GridBagLayout());
		temperaturePanel.setBorder(new javax.swing.border.TitledBorder(
				"Temperature (K)"));
		java.awt.GridBagConstraints gbc1 = new java.awt.GridBagConstraints();
		gbc1.gridx = 0;
		gbc1.gridy = 0;
		gbc1.gridwidth = 1;
		temperaturePanel.add(tSelect.graphic(null), gbc1);
		gbc1.gridx = 0;
		gbc1.gridy = 1;
		temperaturePanel.add(tBox.graphic(null), gbc1);

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
		JPanel topPanel = new JPanel(new java.awt.GridBagLayout());
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
		//    JPanel controlPanel = new JPanel(new java.awt.GridLayout(0,1));
		JPanel controlPanel = new JPanel(new java.awt.GridBagLayout());
		gbc2 = new java.awt.GridBagConstraints();
		gbc2.gridx = 0;
		gbc2.gridy = java.awt.GridBagConstraints.RELATIVE;
		/// controlPanel.add(startPanel,gbc2);
		/// controlPanel.add(temperaturePanel, gbc2);
		controlPanel.add(topPanel, gbc2);
		controlPanel.add(speciesEditors, gbc2);
		/// controlPanel.add(sliderPanel, gbc2);

		panel.add(controlPanel, java.awt.BorderLayout.WEST);
		panel.add(displayPanel, java.awt.BorderLayout.EAST);

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
		 * displayPhase1.graphic(null).setBackground(panelColor);
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

	//=================================================================
	//panel containing species-editing devices

	class MySpeciesEditor extends javax.swing.JPanel {

		//	public DeviceSlider nSlider;
		public DeviceNSelector nSlider;

		public SpeciesAgent species;

		public final javax.swing.JTextField mass = new javax.swing.JTextField(
				"40");

		//    public java.awt.TextField mass = new java.awt.TextField("40");

        public MySpeciesEditor(final ReactionEquilibrium sim, SpeciesAgent s, String label) {
            super();
            species = s;
            nSlider = new DeviceNSelector(sim.getController());
            nSlider.setResetAction(new SimulationRestart(sim));
            nSlider.setSpeciesAgent(species);
            //nSlider.setDisplayPhase(DisplayPhase1);
            nSlider.setMinimum(0);
            nSlider.setMaximum(40);
            nSlider.setPostAction(new Action() {
                public void actionPerformed() {
                    IAtom[] agents = sim.getAgents(sim.phase1);
                    AtomIteratorLeafAtoms iter = new AtomIteratorLeafAtoms(sim.phase1);
                    iter.reset();
                    while (iter.hasNext()) {
                        //                      System.out.println(iter.peek().toString());
                        agents[iter.nextAtom().getGlobalIndex()] = null;
                    }
                    try {
                    	sim.integratorHard1.reset();
                    } catch(ConfigurationOverlapException e) {}
                    displayPhase1.repaint();
                }
           });
//            nSlider.addChangeListener(new javax.swing.event.ChangeListener() {
//                public void stateChanged(javax.swing.event.ChangeEvent evt) {
//                    AtomIteratorLeafAtoms iter = new AtomIteratorLeafAtoms(sim.phase1);
//                    iter.reset();
//                    while (iter.hasNext()) {
//                        //                      System.out.println(iter.peek().toString());
//                        Atom[] a = (Atom[])iter.nextAtom().allatomAgents[sim.idx];
//                        a[0] = null;
//                        a[1] = null;
//                    }
//                     try {
//                         sim.integratorHard1.reset();
//                     } catch(ConfigurationOverlapException e) {}
//                     displayPhase1.repaint();
//                }
//            });
			//            nSlider = new DeviceSlider(new NMoleculeModulator(s));
			//            nSlider.setShowBorder(true);
			//// nSlider.setLabel(label);
			//			nSlider.setLabel("Atom count");
			//	        nSlider.setMinimum(0);
			//	        nSlider.setMaximum(40);
			//	        nSlider.getSlider().setSnapToTicks(true);
			//	        nSlider.graphic(null).setSize(new java.awt.Dimension(40,30));

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
					((ElementSimple)((AtomTypeLeaf)species.getType().getSpecies().getMoleculeType()).getElement()).setMass(newMass);
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
			getContentPane().add(new ReactionEquilibriumGraphic(new ReactionEquilibrium()).panel);
		}
	}//end of Applet

	//---------------------------------------------

	class DiameterModifier implements Modifier {
		P2SquareWellBonded potentialRR, potentialRB, potentialBB;

		SpeciesSpheresMono speciesR, speciesB;

		DisplayPhase display;

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

		public void setDisplay(DisplayPhase display) {
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
    protected DisplayPhase displayPhase1;
}