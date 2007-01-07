package etomica.modules.osmosis;

import java.awt.Color;
import java.awt.Graphics;

import javax.swing.JPanel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import etomica.action.SimulationRestart;
import etomica.atom.AtomTypeSphere;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.DataSourceCountTime;
import etomica.data.meter.MeterLocalMoleFraction;
import etomica.data.meter.MeterTemperature;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceThermoSelector;
import etomica.graphics.DeviceTrioControllerButton;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxesCAE;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayTimer;
import etomica.graphics.Drawable;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.DisplayBox.LabelType;
import etomica.integrator.Integrator;
import etomica.integrator.IntervalActionAdapter;
import etomica.lattice.LatticeCubicSimple;
import etomica.math.geometry.Rectangle;
import etomica.modifier.Modifier;
import etomica.potential.P1HardBoundary;
import etomica.potential.P2HardSphere;
import etomica.space.Vector;
import etomica.space2d.Vector2D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Dimension;
import etomica.units.Kelvin;
import etomica.units.Length;
import etomica.units.Unit;
import etomica.util.Constants.CompassDirection;

/**
 * Osmosis module.
 * @author Jhumpa Adhikari
 * @author Andrew Schultz
 */

public class Osmosis {

    public DisplayPhase displayPhase;
    public DataSourceCountTime cycles;
    public DisplayBox displayCycles;
    public MeterOsmoticPressure osmosisPMeter;
    public MeterLocalMoleFraction moleFraction;
    public OsmosisSim sim;
    public JPanel panel;

    public Osmosis() {

        sim = new OsmosisSim();
        sim.activityIntegrate.setDoSleep(true);
        sim.activityIntegrate.setSleepPeriod(1);
        sim.register(sim.integrator);

        final int thickness = 4;

        Unit tUnit = Kelvin.UNIT;

        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicSimple(2, 1.0));

        DeviceTrioControllerButton button = new DeviceTrioControllerButton(sim);
        SimulationRestart simRestart = new SimulationRestart(sim);
        
        simRestart.setConfiguration(config);
        button.getReinitButton().setAction(simRestart);
        button.setShape("VERTICAL"); //three choices "HORIZONTAL","AUTOMATIC"           

	    //display of phase
	    displayPhase = new DisplayPhase(sim.phase);
        ColorSchemeByType colorScheme = new ColorSchemeByType();

        colorScheme.setColor(sim.speciesA.getMoleculeType(), Color.blue);
        colorScheme.setColor(sim.speciesB.getMoleculeType(), Color.red);
        displayPhase.setColorScheme(colorScheme);
        displayPhase.setAlign(1,DisplayPhase.BOTTOM);
        displayPhase.getOriginShift()[0] = thickness;
        displayPhase.getOriginShift()[1] = -thickness;
        displayPhase.addDrawable(new MyWall());


        cycles = new DataSourceCountTime();
        displayCycles = new DisplayTimer(sim.integrator);
        displayCycles.setLabelType(LabelType.BORDER);
        displayCycles.setLabel("Time");
	    displayCycles.setPrecision(6);	

        osmosisPMeter = new MeterOsmoticPressure(sim.getSpace(), new P1HardBoundary[]{sim.boundaryHardLeftA}, 
                new P1HardBoundary[]{sim.boundaryHardRightA, sim.boundaryHardB});
        osmosisPMeter.setIntegrator(sim.integrator);
        AccumulatorAverage osmosisPMeterAvg = new AccumulatorAverage(sim);
        DataPump pump = new DataPump(osmosisPMeter, osmosisPMeterAvg);
        sim.register(osmosisPMeter, pump);
        IntervalActionAdapter adapter = new IntervalActionAdapter(pump);
        sim.integrator.addListener(adapter);
        DisplayBoxesCAE dBox = new DisplayBoxesCAE();
        dBox.setAccumulator(osmosisPMeterAvg);
        dBox.setPrecision(6);

        moleFraction = new MeterLocalMoleFraction();
//        moleFraction = new MeterMoleFraction();
        moleFraction.setPhase(sim.phase);
        Vector dimensions = sim.phase.getBoundary().getDimensions();
        moleFraction.setShape(new Rectangle(sim.getSpace(), dimensions.x(0)*0.5, dimensions.x(1)));
        moleFraction.setShapeOrigin(new Vector2D(dimensions.x(0)*0.25, 0));
        moleFraction.setSpecies(sim.speciesB);
        AccumulatorAverage moleFractionAvg = new AccumulatorAverage(sim);
        pump = new DataPump(moleFraction, moleFractionAvg);
        sim.register(moleFraction, pump);
        adapter = new IntervalActionAdapter(pump);
        sim.integrator.addListener(adapter);
        DisplayBoxesCAE mfBox = new DisplayBoxesCAE();
        mfBox.setAccumulator(moleFractionAvg);
        mfBox.setPrecision(8);
	    
	    DeviceThermoSelector tSelect = new DeviceThermoSelector(sim, sim.integrator);
	    tSelect.setTemperatures(new double[] {50.,100.,300.,600.,1000.});
	    tSelect.setUnit(tUnit);
	    tSelect.setSelected(0); //sets adiabatic as selected temperature
		MeterTemperature thermometer = new MeterTemperature();
		thermometer.setPhase(sim.phase);
		DisplayBox tBox = new DisplayBox();
        pump = new DataPump(thermometer, tBox);
        adapter = new IntervalActionAdapter(pump);
        sim.integrator.addListener(adapter);
		tBox.setUnit(tUnit);
		tBox.setLabel("Measured value");
		tBox.setLabelPosition(CompassDirection.NORTH);
	    tSelect.getLabel().setText("Set value");

        DeviceNSelector nASelector = new DeviceNSelector(sim.getController());
        nASelector.setResetAction(simRestart);
        nASelector.setSpeciesAgent(sim.phase.getAgent(sim.speciesA));
        nASelector.getSlider().addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                displayPhase.repaint();
            }
        });
        nASelector.setMaximum(30);
        
        DeviceNSelector nBSelector = new DeviceNSelector(sim.getController());
        nBSelector.setResetAction(simRestart);
        nBSelector.setSpeciesAgent(sim.phase.getAgent(sim.speciesB));
        nBSelector.getSlider().addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent evt) {
                displayPhase.repaint();
            }
        });
        nBSelector.setMaximum(10);
	    
	    //DevicePotentialSelector potentialmeter = new DevicePotentialSelector(this);
        
        DiameterModifier diameterModifier = new DiameterModifier(sim.potentialAA,sim.potentialBB,sim.potentialAB,
                                                            sim.boundarySemiB,
                                                            sim.boundaryHardTopBottomA,sim.boundaryHardLeftA,sim.boundaryHardRightA,
                                                            sim.boundaryHardB,
                                                            sim.speciesA,sim.speciesB);
        diameterModifier.setDisplay(displayPhase);
        DeviceSlider sliderDiameter = new DeviceSlider(sim.getController());
        sliderDiameter.setModifier(diameterModifier);
//        sliderDiameter.getSlider().setMaximum(4);
//        sliderDiameter.getSlider().setMinimum(0);
//        sliderDiameter.getSlider().setValue(3);
//        sliderDiameter.getSlider().setPaintTicks(true);
//        sliderDiameter.setPrecision(2);
//        sliderDiameter.setNMajor(1);
//	    sliderDiameter.getSlider().setMajorTickSpacing(1);
//	    sliderDiameter.getSlider().setLabelTable(sliderDiameter.getSlider().createStandardLabels(1));

		sliderDiameter.setPrecision(2);
		sliderDiameter.setMaximum(4);
		sliderDiameter.setMinimum(0);
		sliderDiameter.setValue(3);
//		sliderDiameter.setPaintTicks(true);
		sliderDiameter.setNMajor(4);
//		sliderDiameter.setMajorTickSpacing(1);
//		sliderDiameter.getSlider().setLabelTable(sliderDiameter.getSlider().createStandardLabels(1));
        
        //************* Lay out components ****************//

        

        //tabbed pane for the big displays

    	final javax.swing.JTabbedPane displayPanel = new javax.swing.JTabbedPane();

    	javax.swing.JPanel displayPhasePanel = new javax.swing.JPanel(new java.awt.BorderLayout());

    	displayPhasePanel.add(displayPhase.graphic(),java.awt.BorderLayout.WEST);

//    	displayPhasePanel.add(scaleSliderPanel,java.awt.BorderLayout.EAST);

    	displayPanel.add(displayPhase.getLabel(), displayPhasePanel);


        displayPanel.addChangeListener(

            new javax.swing.event.ChangeListener() {

                public void stateChanged(javax.swing.event.ChangeEvent event) {

                    displayPanel.invalidate();

                    displayPanel.validate();

                }

        });

        JPanel startPanel = new JPanel(new java.awt.GridBagLayout());

        java.awt.GridBagConstraints gbc0 = new java.awt.GridBagConstraints();

      //  startPanel.setBorder(new javax.swing.border.TitledBorder("Control"));

        gbc0.gridx = 0; gbc0.gridy = 0;

        startPanel.add(button.graphic(), gbc0);

        

        //panel for the temperature control/display

        JPanel cyclesPanel = new JPanel(new java.awt.FlowLayout());

        cyclesPanel.setBorder(new javax.swing.border.TitledBorder("Cycles"));

        cyclesPanel.add(displayCycles.graphic(null));

        JPanel temperaturePanel = new JPanel(new java.awt.GridBagLayout());
        temperaturePanel.setBorder(new javax.swing.border.TitledBorder("Temperature (K)"));
        java.awt.GridBagConstraints gbc1 = new java.awt.GridBagConstraints();
        gbc1.gridx = 0;  gbc1.gridy = 1;
        gbc1.gridwidth = 1;
        temperaturePanel.add(tSelect.graphic(null),gbc1);
        gbc1.gridx = 1;  gbc1.gridy = 1;
        temperaturePanel.add(tBox.graphic(null),gbc1);
        

        //panel for the meter displays

        JPanel osmoticPanel = new JPanel(new java.awt.FlowLayout());
        osmoticPanel.setBorder(new javax.swing.border.TitledBorder("Osmotic Pressure (PV/Nk)"));
        osmoticPanel.add(dBox.graphic(null));

        JPanel moleFractionPanel = new JPanel(new java.awt.FlowLayout());
        javax.swing.border.TitledBorder titleBorder = new javax.swing.border.TitledBorder("Mole Fraction (nSolute/nSolution)");
        moleFractionPanel.setBorder(titleBorder);
        moleFractionPanel.add(mfBox.graphic(null));
        
        javax.swing.JTabbedPane tabPaneMeter = new javax.swing.JTabbedPane();
        tabPaneMeter.addTab("Osmotic Pressure", osmoticPanel);
        tabPaneMeter.addTab("Mole Fraction", moleFractionPanel);
        tabPaneMeter.addTab("Cycles", cyclesPanel);

        //panel for sliders

        JPanel sliderPanelA = new JPanel(new java.awt.GridLayout(0,1));
        nASelector.setShowBorder(false);
        sliderPanelA.add(nASelector.graphic(null));
        sliderPanelA.setBorder(new javax.swing.border.TitledBorder("Set "+nASelector.getLabel()));

        
        JPanel sliderPanelB = new JPanel(new java.awt.GridLayout(0,1));
        nBSelector.setShowBorder(false);
        sliderPanelB.add(nBSelector.graphic(null));
        sliderPanelB.setBorder(new javax.swing.border.TitledBorder("Set "+nBSelector.getLabel()));


        JPanel sliderDiaPanel = new JPanel(new java.awt.GridLayout(0,1));
        sliderDiameter.setShowBorder(false);
        sliderDiaPanel.add(sliderDiameter.graphic(null));
        sliderDiaPanel.setBorder(new javax.swing.border.TitledBorder("Set Diameter"));

        
        javax.swing.JTabbedPane tabPaneSliders = new javax.swing.JTabbedPane();
        tabPaneSliders.addTab(nASelector.getLabel(), sliderPanelA);
        tabPaneSliders.addTab(nBSelector.getLabel(), sliderPanelB);
        tabPaneSliders.addTab("Diameter", sliderDiaPanel);
        
        
        //panel for all the controls

//        JPanel controlPanel = new JPanel(new java.awt.GridLayout(0,1));

        JPanel controlPanel = new JPanel(new java.awt.GridBagLayout());

        java.awt.GridBagConstraints gbc2 = new java.awt.GridBagConstraints();

        gbc2.gridx = 0;

        gbc2.gridy = java.awt.GridBagConstraints.RELATIVE;

        controlPanel.add(startPanel,gbc2);

        /*controlPanel.add(cyclesPanel, gbc2);
        controlPanel.add(osmoticPanel, gbc2);
        controlPanel.add(moleFractionPanel, gbc2);

        controlPanel.add(sliderPanelA, gbc2);
        controlPanel.add(sliderPanelB, gbc2);
        controlPanel.add(sliderDiaPanel, gbc2);*/
        
        javax.swing.JTabbedPane tabPane = new javax.swing.JTabbedPane();
        tabPane.addTab("Meters", tabPaneMeter);
        tabPane.addTab("Sliders", tabPaneSliders);
        
        controlPanel.add(temperaturePanel, gbc2);
        controlPanel.add(tabPane,gbc2);
        
   /*     JPanel potentialPanel = new JPanel();
        potentialPanel.add(potentialmeter.graphic(null));
	    potentialPanel.setBorder(new javax.swing.border.TitledBorder("Potential selection"));
        controlPanel.add(potentialPanel,gbc2);*/

        
        panel = new JPanel();
        panel.setLayout(new java.awt.BorderLayout());
        panel.add(controlPanel, java.awt.BorderLayout.WEST);

        panel.add(displayPanel, java.awt.BorderLayout.EAST);

        Thread repainter = new Thread() {
            public void run() {
                while (true) {
                    panel.repaint();
                    try{Thread.sleep(10);}
                    catch(InterruptedException e){}
                }
            }
        };
        repainter.start();
    
    }

	//drawable that puts a line down the middle of the box, where the
	//semipermeable membrane potential acts
    protected class MyWall implements Drawable {
    	public void draw(Graphics g, int[] origin, double scale) {
    		int x1 = origin[0]+(int)(0.5*scale*sim.phase.getBoundary().getDimensions().x(0));
    		int y1 = origin[1];
			int h = (int)(scale*sim.phase.getBoundary().getDimensions().x(1));
			int w = 4;
			g.setColor(java.awt.Color.green);
    		g.fillRect(x1-w, y1, w, h);
    	}
    }

    

    public static void main(String[] args) {
        Osmosis osmosis = new Osmosis();
        SimulationGraphic.makeAndDisplayFrame(osmosis.panel);
    }

    
    public class DiameterModifier implements Modifier {
        
        P2HardSphere potentialAA,potentialBB,potentialAB;
        P1HardWall membraneB;
        P1HardBoundary boundaryHard1A,boundaryHard2A,boundaryHard3A;
        P1HardBoundary boundaryHardB;
        SpeciesSpheresMono speciesA,speciesB;
        DisplayPhase display;
        Integrator integrator;
        
        public DiameterModifier(P2HardSphere potentialAA ,P2HardSphere potentialBB ,P2HardSphere potentialAB ,
                          P1HardWall membraneB,
                          P1HardBoundary boundaryHard1A, P1HardBoundary boundaryHard2A, P1HardBoundary boundaryHard3A, 
                          P1HardBoundary boundaryHardB,
                          SpeciesSpheresMono speciesA ,SpeciesSpheresMono speciesB) {
            this.potentialAA = potentialAA;
            this.potentialBB = potentialBB;
            this.potentialAB = potentialAB;
            this.membraneB = membraneB;
            this.boundaryHard1A = boundaryHard1A;
            this.boundaryHard2A = boundaryHard2A;
            this.boundaryHard3A = boundaryHard3A;
            this.boundaryHardB = boundaryHardB;
            this.speciesA = speciesA;
            this.speciesB = speciesB;
        }
        
        public String getLabel() {
            return "a label";
        }
        
        public Dimension getDimension() {
            return Length.DIMENSION;
        }
        
        public void setValue(double d) {
            potentialAA.setCollisionDiameter(d);
            potentialBB.setCollisionDiameter(d);
            potentialAB.setCollisionDiameter(d);
            boundaryHard1A.setCollisionRadius(0.5*d);
            boundaryHard2A.setCollisionRadius(0.5*d);
            boundaryHard3A.setCollisionRadius(0.5*d);
            membraneB.setCollisionRadius(0.5*d);
            boundaryHardB.setCollisionRadius(0.5*d);
            ((AtomTypeSphere)speciesA.getMoleculeType()).setDiameter(d);
            ((AtomTypeSphere)speciesB.getMoleculeType()).setDiameter(d);
            if (display != null) {
                display.repaint();
            }
            if (integrator != null) {
                try {
                    integrator.reset();
                }
                catch (ConfigurationOverlapException e) {
                    //overlaps are likely after increasing the diameter
                }
            }
        }
        
        public double getValue() {
            return ((AtomTypeSphere)speciesA.getMoleculeType()).getDiameter();
        }
        
        public void setDisplay(DisplayPhase newDisplay){
            display = newDisplay;
        }
        
        public void setIntegrator(Integrator newIntegrator) {
            integrator = newIntegrator;
        }
    }

    public static class Applet extends javax.swing.JApplet {
	    public void init() {
		    getContentPane().add(new Osmosis().panel);
	    }

        private static final long serialVersionUID = 1L;
    }

}








