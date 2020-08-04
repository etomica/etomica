/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.insertion;

import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.data.*;
import etomica.data.histogram.HistogramNotSoSimple;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.*;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorHard.Agent;
import etomica.integrator.IntegratorMD;
import etomica.math.DoubleRange;
import etomica.math.function.Function;
import etomica.modifier.Modifier;
import etomica.potential.P2HardSphere;
import etomica.potential.P2SquareWell;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.units.Pixel;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.awt.event.ItemListener;
import java.util.ArrayList;

public class InsertionGraphic extends SimulationGraphic {

    private final static String APP_NAME = "Square-Well Molecular Dynamics";
    private final static int REPAINT_INTERVAL = 1;
    protected DeviceThermoSlider tempSlider;
    public ItemListener potentialChooserListener;
    public JComboBox potentialChooser;
    protected P2HardSphere potentialHS;
    protected P2SquareWell potentialSW;
    public DeviceBox lamBox;
    public double lambda;
    protected Insertion sim;
    public MeterWidom meterWidom;
    protected DisplayPlotXChart widomHistPlot;
    public DisplayPlotXChart widom2Plot;
    
    public InsertionGraphic(final Insertion simulation) {

    	super(simulation, TABBED_PANE, APP_NAME, REPAINT_INTERVAL);

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();

        final IAction resetDataAction = new IAction() {
            public void actionPerformed() {
                getController().getResetAveragesButton().press();
            }
        };

    	this.sim = simulation;

        lambda = 1.5;
        double sigma = 1.0;

        if (sim.getSpace().D() == 2) {
            getDisplayBox(sim.box).setPixelUnit(new Pixel(400/sim.box.getBoundary().getBoxSize().getX(1)));
        }
        else {
            getDisplayBox(sim.box).setPixelUnit(new Pixel(40/sim.box.getBoundary().getBoxSize().getX(1)));
        }
        

        sim.activityIntegrate.setSleepPeriod(0);
       
        //combo box to select potentials
        final String repulsionOnly = "Repulsion only";
        final String repulsionAttraction = "Repulsion and attraction";
        potentialChooser = new javax.swing.JComboBox(new String[] {
                repulsionOnly, repulsionAttraction});

        lamBox = new DeviceBox();
        // Unselectable because "Ideal gas" is selected initially
        potentialChooser.setSelectedIndex(0);
        lamBox.setEditable(false);

        // Simulation Time

        //temperature selector
        tempSlider = new DeviceThermoSlider(sim.getController(), sim.integrator);
        tempSlider.setIsothermalButtonsVisibility(false);
        tempSlider.setPrecision(1);
        tempSlider.setMinimum(0.0);
        tempSlider.setMaximum(3.0);
        tempSlider.setSliderMajorValues(3);
        tempSlider.setIsothermal();

        JPanel statePanel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc2 = new GridBagConstraints();
        gbc2.gridx = 0;  gbc2.gridy = 0;
        statePanel.add(tempSlider.graphic(), gbc2);

        GridBagConstraints vertGBC = SimulationPanel.getVertGBC();

        JPanel potentialPanel = new JPanel(new GridBagLayout());
        potentialPanel.add(potentialChooser,vertGBC);
        JPanel parameterPanel = new JPanel(new GridLayout(0,1));
        parameterPanel.add(lamBox.graphic());
        potentialPanel.add(parameterPanel,vertGBC);
        
        //
        // Tabbed pane for state, potential, controls pages
        //
        JTabbedPane setupPanel = new JTabbedPane();
        setupPanel.add(statePanel, "State");
        setupPanel.add(potentialPanel, "Potential");

        potentialSW = new P2SquareWell(sim.getSpace(),sigma,lambda,1.0,true);
        potentialHS = new P2HardSphere(sim.getSpace(),sigma,true);
        
        if(potentialChooserListener != null) potentialChooser.removeItemListener(potentialChooserListener);
        
        potentialChooserListener = new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                if(evt.getStateChange() == java.awt.event.ItemEvent.DESELECTED) return; 
                setPotential((String)evt.getItem());
                if(evt.getItem() == repulsionOnly) {
                    lamBox.setEditable(false);
                }
                else {
                    lamBox.setEditable(true);   
                }
            }
        };
        potentialChooser.addItemListener(potentialChooserListener);
        setPotential((String)potentialChooser.getSelectedItem());

        Modifier lamModifier = new Modifier() {
            public void setValue(double newValue) {
                potentialSW.setLambda(newValue);
                sim.potentialGhost.setLambda(newValue);
            }
            
            public double getValue() {
                return potentialSW.getLambda();
            }
            
            public String getLabel() {
                return "lambda";
            }
            
            public Dimension getDimension() {
                return Length.DIMENSION;
            }
        };
        lamBox.setModifier(lamModifier);
        lamBox.setController(sim.getController());

        //display of box, timer
        ColorSchemeByType colorScheme = new ColorSchemeByType();
        colorScheme.setColor(sim.species.getLeafType(),java.awt.Color.red);
        getDisplayBox(sim.box).setColorScheme(new ColorScheme() {
            MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.integrator.getPotentialMaster());
            public Color getAtomColor(IAtom a) {
                if ((sim.integrator.getEventManager().firingEvent() || !sim.getController().isActive()) && a.getType().getSpecies() == sim.speciesGhost) {
                    sim.potentialGhost.setEpsilonCore(-Double.POSITIVE_INFINITY);
                    if (potentialChooser.getSelectedItem().equals("Repulsion and attraction")) {
                        sim.potentialGhost.setEpsilonWell(1.0);
                    }
                    meterPE.setBox(sim.box);
                    meterPE.setTarget(a);
                    double u = meterPE.getDataAsScalar();
                    if (u==Double.POSITIVE_INFINITY) {
                        sim.potentialGhost.setEpsilonCore(0.0);
                        sim.potentialGhost.setEpsilonWell(0.0);
                        return Color.RED;
                    }
                    if (u < -6) {
                        u = -6;
                    }
                    sim.potentialGhost.setEpsilonCore(0.0);
                    sim.potentialGhost.setEpsilonWell(0.0);
                    return new Color(0.0f, (float)((6+u)/6.0), (float)(-u/6.0));
                }
                return Color.BLACK;
            }
        });


        final DeviceNSelector nSlider = new DeviceNSelector(sim.getController());
        nSlider.setResetAction(new IAction() {
            SimulationRestart simRestart = new SimulationRestart(sim);

            public void actionPerformed() {
                sim.box.setNMolecules(sim.speciesGhost, 0);
                simRestart.actionPerformed();
                sim.box.setNMolecules(sim.speciesGhost, 1);
                sim.integrator.reset();
                meterWidom.reset();
            }
        });
        nSlider.setSpecies(sim.species);
        nSlider.setBox(sim.box);
        nSlider.setMinimum(0);
        nSlider.setMaximum(sim.getSpace().D() == 3 ? 500 : 168);
        nSlider.setLabel("Number of Atoms");
        nSlider.setShowBorder(true);
        nSlider.setShowValues(true);
        // add a listener to adjust the thermostat interval for different
        // system sizes (since we're using ANDERSEN_SINGLE.  Smaller systems 
        // don't need as much thermostating.
        nSlider.setPostAction(new IAction() {
            
            public void actionPerformed() {
                final int n = (int)nSlider.getValue() > 0 ? (int)nSlider.getValue() : 1;
                sim.integrator.setThermostatInterval(n > 40 ? 1 : 40/n);

                getDisplayBox(sim.box).repaint();
            }
        });
        JPanel nSliderPanel = new JPanel(new GridLayout(0,1));
        nSliderPanel.setBorder(new TitledBorder(null, "Number of Molecules", TitledBorder.CENTER, TitledBorder.TOP));
        nSlider.setShowBorder(false);
        nSlider.setNMajor(4);
        nSliderPanel.add(nSlider.graphic());
        gbc2.gridx = 0;  gbc2.gridy = 1;
        statePanel.add(nSliderPanel, gbc2);

        //************* Lay out components ****************//

        getDisplayBox(sim.box).setScale(0.7);


	    final IAction temperatureAction = new IAction() {
		    public void actionPerformed() {
		        resetDataAction.actionPerformed();
		    }
		};
		tempSlider.setSliderPostAction(temperatureAction);
        tempSlider.setRadioGroupPostAction(temperatureAction);

        final IAction resetAction = new IAction() {
        	public void actionPerformed() {
        	    sim.integrator.reset();
                meterWidom.reset();
        		getDisplayBox(sim.box).graphic().repaint();
        	}
        };

        this.getController().getReinitButton().setPostAction(new IAction() {
            public void actionPerformed() {
                sim.integrator.setThermostat(IntegratorMD.ThermostatType.ANDERSEN);
                sim.integrator.doThermostat();
                sim.integrator.setThermostat(IntegratorMD.ThermostatType.ANDERSEN_SCALING);
                resetAction.actionPerformed();
                meterWidom.reset();
            }
        });
        this.getController().getResetAveragesButton().setPostAction(resetAction);

        DeviceDelaySlider delaySlider = new DeviceDelaySlider(sim.getController());
        delaySlider.setMaxSleep(1000);

        getPanel().controlPanel.add(setupPanel, vertGBC);
        getPanel().controlPanel.add(delaySlider.graphic(), vertGBC);

        meterWidom = new MeterWidom();
        sim.integrator.addCollisionListener(meterWidom);
        AccumulatorAverageCollapsing accWidom = new AccumulatorAverageCollapsing();
        DataPumpListener pumpWidom = new DataPumpListener(meterWidom, accWidom);
        sim.integrator.getEventManager().addListener(pumpWidom);
        dataStreamPumps.add(pumpWidom);
        DisplayTextBoxesCAE widomBoxes = new DisplayTextBoxesCAE();
        widomBoxes.setAccumulator(accWidom);
        widomBoxes.setDoShowCurrent(false);
        widomBoxes.setLabel("Widom insertion");
        widomBoxes.setPrecision(6);
        getPanel().controlPanel.add(widomBoxes.graphic(), vertGBC);
        
        DataProcessor dpMu = new DataProcessorFunction(new Function() {
            
            public double f(double x) {
                if (x==0) return Double.NaN;
                return -sim.integrator.getTemperature()*Math.log(x);
            }
        });
        accWidom.addDataSink(dpMu, new AccumulatorAverage.StatType[]{accWidom.AVERAGE});
        
        AccumulatorHistory historyMu = new AccumulatorHistory(new HistoryCollapsingDiscard());
        dpMu.setDataSink(historyMu);
        historyMu.setTimeDataSource(new DataSourceCountTime(sim.integrator));
        DisplayPlotXChart historyWidomPlot = new DisplayPlotXChart();
        historyMu.setDataSink(historyWidomPlot.getDataSet().makeDataSink());
        historyWidomPlot.setLabel("History");
        historyWidomPlot.setDoLegend(false);
        add(historyWidomPlot);
        
        HistogramDataSource widomHist = new HistogramDataSource(meterWidom.hist);
        widomHistPlot = new DisplayPlotXChart();
        DataPumpListener widomHistPump = new DataPumpListener(widomHist, widomHistPlot.getDataSet().makeDataSink(), 10);
        sim.integrator.getEventManager().addListener(widomHistPump);
        widomHistPlot.setLabel("widom");
        widomHistPlot.getPlot().setXRange(0,1);
        widomHistPlot.getPlot().setYLog(true);
        widomHistPlot.setDoLegend(false);
        widomHistPlot.setXLabel("u");
        add(widomHistPlot);
        
        DataProcessor widomWeight = new DataProcessor() {
            DataFunction data;

            protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
                data = (DataFunction)inputDataInfo.makeData();
                dataInfo = inputDataInfo.getFactory().makeDataInfo();
                dataInfo.addTag(tag);
                return dataInfo;
            }
            
            protected IData processData(IData inputData) {
                double[] y1 = ((DataFunction)inputData).getData();
                double[] y2 = data.getData();
                double[] x = ((DataInfoFunction)dataInfo).getXDataSource().getIndependentData(0).getData();
                for (int i=0; i<x.length; i++) {
                    if (x[i] > 0) {
                        y2[i] = 0;
                    }
                    else {
                        y2[i] = y1[i]*Math.exp(-x[i]/sim.integrator.getTemperature());
                    }
                }
                return data;
            }
        };
        DataPumpListener widomHist2Pump = new DataPumpListener(widomHist, widomWeight, 10);
        widom2Plot = new DisplayPlotXChart();
        widomWeight.setDataSink(widom2Plot.getDataSet().makeDataSink());
        sim.integrator.getEventManager().addListener(widomHist2Pump);
        widom2Plot.getPlot().setXRange(0,1);
        widom2Plot.setXLabel("u");
        widom2Plot.setDoLegend(false);
        add(widom2Plot);
        
        JPanel myPlotPanel = new JPanel(new GridLayout(0, 1));
        myPlotPanel.add(widomHistPlot.graphic());
        myPlotPanel.add(widom2Plot.graphic());
        JScrollPane plotsPane = new JScrollPane(myPlotPanel);
        java.awt.Dimension d = widomHistPlot.getPlot().getPreferredSize();
        d.width -= 100;
        d.height = 210;
        d.width += 40;
        d.height = d.height * 2 + 40;
        plotsPane.setPreferredSize(d);
        getPanel().tabbedPane.add("Histograms", plotsPane);

    }

    public void setPotential(String potentialDesc) {
        final boolean HS = potentialDesc.equals("Repulsion only"); 
        final boolean SW = potentialDesc.equals("Repulsion and attraction"); 
        sim.getController().doActionNow( new IAction() {
            public void actionPerformed() {
                if (HS) {
                    potentialHS.setBox(sim.box);
                    sim.potentialWrapper.setWrappedPotential(potentialHS);
                    if (widomHistPlot!=null) {
                        widomHistPlot.getPlot().setXRange(0,1);
                        widom2Plot.getPlot().setXRange(0,1);
                    }
                }
                else if (SW) {
                    potentialSW.setBox(sim.box);
                    sim.potentialWrapper.setWrappedPotential(potentialSW);
                    widomHistPlot.getPlot().setXRange(-7,1);
                    widom2Plot.getPlot().setXRange(-7,1);
                }
                else {
                    throw new RuntimeException("oops");
                }
                try {
                    sim.integrator.reset();
                } catch(ConfigurationOverlapException e) {}
                if (meterWidom!=null) meterWidom.reset();
            }
        });
    }
    
    public class MeterWidom extends DataSourceScalar implements IntegratorHard.CollisionListener {
        protected final MeterPotentialEnergy meterPE;
        protected double lastTime, lastTime2;
        protected double sum;
        protected int currentWells, currentCores;
        public HistogramNotSoSimple hist;
        protected final Vector dr, dv;
        
        public MeterWidom() {
            super("widom", Null.DIMENSION);
            meterPE = new MeterPotentialEnergy(sim.integrator.getPotentialMaster());
            hist = new HistogramNotSoSimple(12, new DoubleRange(-10.5, 1.5));
            hist.setDoAveraging(false);
            dr = sim.getSpace().makeVector();
            dv = sim.getSpace().makeVector();
            reset();
        }
        
        public void reset() {
            meterPE.setBox(sim.box);
            meterPE.setTarget(sim.box.getMoleculeList(sim.getSpecies(1)).get(0).getChildList().get(0));
            sim.potentialGhost.setEpsilonWell(1.0);
            currentWells = -(int)Math.round(meterPE.getDataAsScalar());
            sim.potentialGhost.setEpsilonWell(0);
            sim.potentialGhost.setEpsilonCore(1);
            currentCores = -(int)Math.round(meterPE.getDataAsScalar());
            sim.potentialGhost.setEpsilonCore(0);
            lastTime = sim.integrator.getCurrentTime();
            lastTime2 = lastTime;
            sum = 0;
            hist.reset();
        }

        public void collisionAction(Agent colliderAgent) {
            IAtom atom = colliderAgent.atom;
            if (atom.getParentGroup().getType().getIndex() != 1) {
                atom = colliderAgent.collisionPartner;
                if (atom == null) return;
            }
            if (atom.getParentGroup().getType().getIndex() == 1) {
                dr.Ev1Mv2(colliderAgent.atom.getPosition(), colliderAgent.collisionPartner.getPosition());
                dv.Ev1Mv2(((IAtomKinetic)colliderAgent.atom).getVelocity(), ((IAtomKinetic)colliderAgent.collisionPartner).getVelocity());
                dr.PEa1Tv1(colliderAgent.collisionTime(), dv);
                double r2 = dr.squared();
                boolean core = (Math.abs(r2-1) < 1e-9); 
                if (core) {
                    if (r2>1) {
                        currentCores--;
                        currentWells++;
                    }
                    else {
                        currentCores++;
                        currentWells--;
                    }
                }
                else {
                    double l = sim.potentialGhost.getLambda();
                    if (r2 < l*l) {
                        currentWells++;
                    }
                    else {
                        currentWells--;
                    }
                }
                if (currentWells < 0 || currentCores < 0) {
                    throw new RuntimeException("oops");
                }
                double u = 0;
                if (currentCores>0) u = Double.POSITIVE_INFINITY;
                else if (potentialChooser.getSelectedItem().equals("Repulsion and attraction")) u = -currentWells;
                double t = sim.integrator.getCurrentTime()+colliderAgent.collisionTime();
                hist.addValue(u==Double.POSITIVE_INFINITY?1:u, (t-lastTime));
                sum += (t-lastTime)*Math.exp(-u/sim.integrator.getTemperature());
                lastTime = t;
            }
        }

        public double getDataAsScalar() {
            double t = sim.integrator.getCurrentTime();
            double v = sum/(t-lastTime2);
            lastTime2 = t;
            sum = 0;
            return v;
        }
    }

    public static void main(String[] args) {
        Space space = Space2D.getInstance();
        if(args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    space = Space3D.getInstance();
                }
            } catch(NumberFormatException e) {}
        }

        InsertionGraphic swmdGraphic = new InsertionGraphic(new Insertion(space));
		SimulationGraphic.makeAndDisplayFrame
		        (swmdGraphic.getPanel(), APP_NAME);
    }
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
	        getRootPane().putClientProperty(
	                        "defeatSystemEventQueueCheck", Boolean.TRUE);
            String dimStr = getParameter("dim");
            int dim = 3;
            if (dimStr != null) {
                dim = Integer.valueOf(dimStr).intValue();
            }
            Space sp = Space.getInstance(dim);
            InsertionGraphic swmdGraphic = new InsertionGraphic(new Insertion(sp));

		    getContentPane().add(swmdGraphic.getPanel());
	    }

        private static final long serialVersionUID = 1L;
    }

}
