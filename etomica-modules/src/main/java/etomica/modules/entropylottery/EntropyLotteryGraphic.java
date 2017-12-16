/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.entropylottery;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagLayout;
import java.util.ArrayList;

import javax.swing.JPanel;

import etomica.action.ActionGroupSeries;
import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPump;
import etomica.data.DataSourceCountSteps;
import etomica.data.DataTag;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.integrator.IntegratorListenerAction;
import etomica.space.Space;
import etomica.space1d.Space1D;
import etomica.units.Pixel;
import etomica.data.history.HistoryCollapsingDiscard;

public class EntropyLotteryGraphic extends SimulationGraphic {

    private static final String APP_NAME = "Entropy Lottery";
    private static final int REPAINT_ACTION = 5;

	protected final EntropyLottery sim;

	public EntropyLotteryGraphic(final EntropyLottery simulation, Space _space) {

		super(simulation, GRAPHIC_ONLY, APP_NAME, REPAINT_ACTION);
        this.sim = simulation;

        sim.activityIntegrate.setSleepPeriod(10);

        ArrayList<DataPump> dataStreamPumps = getController().getDataStreamPumps();
        
        this.getController().getSimRestart().setConfiguration(new ConfigurationZero(space));

        this.getController().getControllerButton().setPostAction(getPaintAction(sim.box));

	    //display of box, timer
        getDisplayBox(sim.box).setPixelUnit(new Pixel(300/sim.box.getBoundary().getBoxSize().getX(0)));
        getDisplayBox(sim.box).setDrawingHeight(300);
        DisplayBoxCanvas1DBins canvas = new DisplayBoxCanvas1DBins(getDisplayBox(sim.box), sim.getController());
        Component controllerButtons = getPanel().graphicsPanel.getComponent(0);
        final DisplayBox displayBox = getDisplayBox(sim.box);
        remove(displayBox);
        getPanel().graphicsPanel.remove(controllerButtons);
        displayBox.setBoxCanvas(canvas);

        MeterEntropy meterEntropy = new MeterEntropy();
        meterEntropy.setBox(sim.box);
        AccumulatorHistory entropyHistory = new AccumulatorHistory(new HistoryCollapsingDiscard(100));
        DataSourceCountSteps stepCounter = new DataSourceCountSteps(sim.integrator);
        entropyHistory.setTimeDataSource(stepCounter);
        DataPump pump = new DataPump(meterEntropy, entropyHistory);
        IntegratorListenerAction pumpListener = new IntegratorListenerAction(pump);
        sim.integrator.getEventManager().addListener(pumpListener);
        pumpListener.setInterval(20);
        dataStreamPumps.add(pump);
        
        final DataSourceProbabilityDensity probabilityDensity = new DataSourceProbabilityDensity();
        probabilityDensity.setBox(sim.box);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(probabilityDensity));
//        sim.integrator.setIntervalActionPriority(probabilityDensity, 0);
        EntropyProcessor entropyProcessor = new EntropyProcessor();
        pump = new DataPump(probabilityDensity, entropyProcessor);
        pumpListener = new IntegratorListenerAction(pump);
        sim.integrator.getEventManager().addListener(pumpListener);
        pumpListener.setInterval(20);
        AccumulatorHistory probabilityEntropyHistory = new AccumulatorHistory(new HistoryCollapsingDiscard(100));
        stepCounter = new DataSourceCountSteps(sim.integrator);
        probabilityEntropyHistory.setTimeDataSource(stepCounter);
        entropyProcessor.setDataSink(probabilityEntropyHistory);
        dataStreamPumps.add(pump);
        canvas.setExtraData(probabilityDensity);
        
        IAction resetEntropy = new IAction() {
            public void actionPerformed() {
                probabilityDensity.reset();
            }
        };
        this.getController().getReinitButton().setPostAction(new ActionGroupSeries(new IAction[]{getPaintAction(sim.box), resetEntropy}));
        
        DisplayPlot entropyPlot = new DisplayPlot();
        entropyHistory.setDataSink(entropyPlot.getDataSet().makeDataSink());
        entropyPlot.setLegend(new DataTag[]{meterEntropy.getTag()}, "measured");
        probabilityEntropyHistory.setDataSink(entropyPlot.getDataSet().makeDataSink());
        entropyPlot.setLegend(new DataTag[]{probabilityDensity.getTag()}, "predicted");
        entropyPlot.getPlot().setTitle("Entropy");
        
        DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
        nSelector.setBox(sim.box);
        nSelector.setSpecies(sim.species);
        nSelector.setLabel("Number of Balls");
        nSelector.setShowBorder(true);
        
        final DeviceSlider nUrnSelector = new DeviceSlider(sim.getController(), new ModifierDimensions(space, sim.box));
        nUrnSelector.setMinimum(3);
        nUrnSelector.setMaximum(30);
        nUrnSelector.setLabel("Number of Urns");
        nUrnSelector.setShowBorder(true);
        IAction resetDisplay = new IAction() {
            public void actionPerformed() {
                double nUrn = nUrnSelector.getValue();
                double a2p = 300.0/nUrn;
                displayBox.setPixelUnit(new Pixel(a2p));
                double yScale = nUrn*nUrn/(6*sim.box.getNMolecules(sim.species));
                if (yScale > 6) {
                    yScale = 6;
                }
                ((DisplayBoxCanvas1DBins)displayBox.canvas).setYScale(yScale);
                displayBox.repaint();
            }
        };

        SimulationRestart restartAction =
        	(SimulationRestart)getController().getReinitButton().getAction();

        nUrnSelector.setPostAction(new ActionGroupSeries(new IAction[]
                                   {resetDisplay, restartAction, resetEntropy}));
        nUrnSelector.doUpdate();

        nSelector.setResetAction(new ActionGroupSeries(new IAction[]
                                  {resetDisplay, restartAction, resetEntropy}));
        nSelector.doUpdate();

        resetDisplay.actionPerformed();

        // SimulationGraphic has already added the box graphic.
        // Remove it from the graphicsPanel and add it to this
        // classes' bigPanel.
        //tabbed pane for the big displays
        JPanel bigPanel = new JPanel(new GridBagLayout());
        
        bigPanel.add(entropyPlot.getPlot(), SimulationPanel.getVertGBC());
        displayBox.canvas.setPreferredSize(new Dimension(500, 300));
        bigPanel.add(displayBox.canvas, SimulationPanel.getVertGBC());
        bigPanel.add(controllerButtons, SimulationPanel.getVertGBC());

        add(nSelector);
        add(nUrnSelector);
        getPanel().graphicsPanel.add(bigPanel);

    }    

    public static void main(String[] args) {
    	Space sp = Space1D.getInstance();
    	EntropyLottery sim = new EntropyLottery(sp);
        EntropyLotteryGraphic entropyLotteryGraphic = new EntropyLotteryGraphic(sim, sp);
		SimulationGraphic.makeAndDisplayFrame(entropyLotteryGraphic.getPanel(), APP_NAME);
    }
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
        	Space sp = Space1D.getInstance();
        	EntropyLottery sim = new EntropyLottery(sp);
            EntropyLotteryGraphic entropyLotteryGraphic = new EntropyLotteryGraphic(sim, sp);
		    getContentPane().add(entropyLotteryGraphic.getPanel());
	    }

        private static final long serialVersionUID = 1L;
    }

}


