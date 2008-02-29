package etomica.modules.entropylottery;

import java.awt.GridLayout;
import java.util.ArrayList;

import javax.swing.JPanel;

import etomica.action.Action;
import etomica.action.ActionGroupSeries;
import etomica.action.SimulationRestart;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPump;
import etomica.data.DataSourceCountSteps;
import etomica.data.DataTag;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.space.Space;
import etomica.space1d.Space1D;
import etomica.units.Pixel;
import etomica.util.HistoryCollapsing;

public class EntropyLotteryGraphic extends SimulationGraphic {

    private static final String APP_NAME = "Entropy Lottery";
    private static final int REPAINT_ACTION = 5;

	protected final EntropyLottery sim;

	public EntropyLotteryGraphic(final EntropyLottery simulation, Space _space) {

		super(simulation, GRAPHIC_ONLY, APP_NAME, REPAINT_ACTION, _space);
        this.sim = simulation;

        sim.activityIntegrate.setSleepPeriod(10);

        ArrayList dataStreamPumps = getController().getDataStreamPumps();
        
        this.getController().getSimRestart().setConfiguration(new ConfigurationZero(space));

        this.getController().getReinitButton().setPostAction(getPaintAction(sim.box));

        this.getController().getControllerButton().setPostAction(getPaintAction(sim.box));

	    //display of box, timer
        getDisplayBox(sim.box).setPixelUnit(new Pixel(300/sim.box.getBoundary().getDimensions().x(0)));
        getDisplayBox(sim.box).setDrawingHeight(300);
        DisplayBoxCanvas1DBins canvas = new DisplayBoxCanvas1DBins(getDisplayBox(sim.box));
        getDisplayBox(sim.box).setBoxCanvas(canvas);

        //tabbed pane for the big displays
        JPanel bigPanel = new JPanel(new GridLayout(2,0));
        
        MeterEntropy meterEntropy = new MeterEntropy();
        meterEntropy.setBox(sim.box);
        AccumulatorHistory entropyHistory = new AccumulatorHistory(new HistoryCollapsing(100));
        DataSourceCountSteps stepCounter = new DataSourceCountSteps(sim.integrator);
        entropyHistory.setTimeDataSource(stepCounter);
        DataPump pump = new DataPump(meterEntropy, entropyHistory);
        sim.integrator.addIntervalAction(pump);
        sim.integrator.setActionInterval(pump, 20);
        dataStreamPumps.add(pump);
        
        DataSourceProbabilityDensity probabilityDensity = new DataSourceProbabilityDensity();
        sim.integrator.addNonintervalListener(probabilityDensity);
        sim.integrator.addIntervalAction(probabilityDensity);
        sim.integrator.setIntervalActionPriority(probabilityDensity, 0);
        EntropyProcessor entropyProcessor = new EntropyProcessor();
        pump = new DataPump(probabilityDensity, entropyProcessor);
        sim.integrator.addIntervalAction(pump);
        sim.integrator.setActionInterval(pump, 20);
        AccumulatorHistory probabilityEntropyHistory = new AccumulatorHistory(new HistoryCollapsing(100));
        stepCounter = new DataSourceCountSteps(sim.integrator);
        probabilityEntropyHistory.setTimeDataSource(stepCounter);
        entropyProcessor.setDataSink(probabilityEntropyHistory);
        dataStreamPumps.add(pump);
        canvas.setExtraData(probabilityDensity);
        
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
        
        final DeviceSlider nUrnSelector = new DeviceSlider(sim.getController(), new ModifierDimensions(sim.box));
        nUrnSelector.setMinimum(3);
        nUrnSelector.setMaximum(30);
        nUrnSelector.setLabel("Number of Urns");
        nUrnSelector.setShowBorder(true);
        Action resetDisplay = new Action() {
            public void actionPerformed() {
                double nUrn = nUrnSelector.getValue();
                double a2p = 300.0/nUrn;
                getDisplayBox(sim.box).setPixelUnit(new Pixel(a2p));
                double yScale = nUrn*nUrn/(6*sim.box.getNMolecules(sim.species));
                if (yScale > 6) {
                    yScale = 6;
                }
                ((DisplayBoxCanvas1DBins)getDisplayBox(sim.box).canvas).setYScale(yScale);
                getDisplayBox(sim.box).repaint();
            }
        };

        SimulationRestart restartAction =
        	(SimulationRestart)getController().getReinitButton().getAction();

        nUrnSelector.setPostAction(new ActionGroupSeries(new Action[]
                                   {resetDisplay, restartAction}));
        nUrnSelector.doUpdate();

        nSelector.setResetAction(new ActionGroupSeries(new Action[]
                                  {resetDisplay, restartAction}));
        nSelector.doUpdate();

        resetDisplay.actionPerformed();

        // SimulationGraphic has already added the box graphic.
        // Remove it from the graphicsPanel and add it to this
        // classes' bigPanel.
        getPanel().graphicsPanel.removeAll();
        bigPanel.add(entropyPlot.getPlot());
        bigPanel.add(getDisplayBox(sim.box).graphic());

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


