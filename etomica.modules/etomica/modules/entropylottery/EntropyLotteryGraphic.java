package etomica.modules.entropylottery;
import java.awt.BorderLayout;

import javax.swing.JPanel;

import etomica.action.Action;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPump;
import etomica.data.DataSourceCountSteps;
import etomica.data.DataTag;
import etomica.graphics.DeviceTrioControllerButton;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntervalActionAdapter;
import etomica.space1d.Space1D;
import etomica.util.HistoryCollapsing;

public class EntropyLotteryGraphic {
    
    public void init(final EntropyLottery sim) {

        sim.getDefaults().blockSize = 100;
        sim.getDefaults().doSleep = true;
        sim.activityIntegrate.setDoSleep(true);
        
        sim.register(sim.integrator);
        
        panel = new JPanel();
        panel.setLayout(new BorderLayout());
        
        DeviceTrioControllerButton control = new DeviceTrioControllerButton(sim);
        
	    //display of phase, timer
	    final DisplayPhase displayPhase = new DisplayPhase(sim.phase,sim.getDefaults().pixelUnit);
        displayPhase.setDrawingHeight(300);
        DisplayPhaseCanvas1DBins canvas = new DisplayPhaseCanvas1DBins(displayPhase);
        displayPhase.setPhaseCanvas(canvas);

        sim.integrator.addListener(new IntervalActionAdapter(new Action() {
            public void actionPerformed() {displayPhase.repaint();}
            public String getLabel() {return "Phase";}
        }));
        
   /*     DisplayTimer timer = new DisplayTimer(integrator);
        timer.setLabelType(DisplayBox.BORDER);
        timer.setUnit(new Unit(LennardJones.Time.UNIT));
	*/    
        DataSourceCountSteps timeCounter = new DataSourceCountSteps();
        sim.integrator.addListener(timeCounter);
		
        //tabbed pane for the big displays
        javax.swing.JPanel bigPanel = new javax.swing.JPanel();
        
    	JPanel displayPanel = new JPanel();
        
        MeterEntropy meterEntropy = new MeterEntropy();
        meterEntropy.setPhase(sim.phase);
        AccumulatorHistory entropyHistory = new AccumulatorHistory(HistoryCollapsing.FACTORY, 100);
        DataSourceCountSteps stepCounter = new DataSourceCountSteps();
        sim.integrator.addListener(stepCounter);
        entropyHistory.setTimeDataSource(stepCounter);
        DataPump pump = new DataPump(meterEntropy, entropyHistory);
        IntervalActionAdapter iaa = new IntervalActionAdapter(pump);
        iaa.setActionInterval(20);
        sim.integrator.addListener(iaa);
        
        DataSourceProbabilityDensity probabilityDensity = new DataSourceProbabilityDensity();
        sim.integrator.addListener(probabilityDensity);
        EntropyProcessor entropyProcessor = new EntropyProcessor();
        pump = new DataPump(probabilityDensity, entropyProcessor);
        iaa = new IntervalActionAdapter(pump);
        iaa.setActionInterval(20);
        sim.integrator.addListener(iaa);
        AccumulatorHistory probabilityEntropyHistory = new AccumulatorHistory(HistoryCollapsing.FACTORY, 100);
        stepCounter = new DataSourceCountSteps();
        sim.integrator.addListener(stepCounter);
        probabilityEntropyHistory.setTimeDataSource(stepCounter);
        entropyProcessor.setDataSink(probabilityEntropyHistory);
        
        DisplayPlot entropyPlot = new DisplayPlot();
        entropyHistory.setDataSink(entropyPlot.getDataSet().makeDataSink());
        entropyPlot.setLegend(new DataTag[]{meterEntropy.getTag()}, "measured");
        probabilityEntropyHistory.setDataSink(entropyPlot.getDataSet().makeDataSink());
        entropyPlot.setLegend(new DataTag[]{probabilityDensity.getTag()}, "predicted");
        displayPanel.add(entropyPlot.getPlot());
        
        
        bigPanel.add(displayPhase.graphic());
        bigPanel.add(displayPanel);
        
        //panel for the start buttons
        JPanel startPanel = (JPanel)control.graphic();
        
        
        JPanel controlPanel = new JPanel(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gbc2 = new java.awt.GridBagConstraints();
        gbc2.gridy = 0;
        gbc2.gridx = java.awt.GridBagConstraints.RELATIVE;
        controlPanel.add(startPanel,gbc2);
        panel.add(controlPanel, java.awt.BorderLayout.NORTH);
        panel.add(bigPanel, java.awt.BorderLayout.EAST);

    }    
    
    public static void main(String[] args) {
            
        EntropyLotteryGraphic entropyLotteryGraphic = new EntropyLotteryGraphic();
        entropyLotteryGraphic.init(new EntropyLottery(Space1D.getInstance()));
		SimulationGraphic.makeAndDisplayFrame(entropyLotteryGraphic.panel);
    }
    
    public static class Applet extends javax.swing.JApplet {

        public void init() {
            EntropyLotteryGraphic entropyLotteryGraphic = new EntropyLotteryGraphic();
            entropyLotteryGraphic.init(new EntropyLottery(Space1D.getInstance()));
		    getContentPane().add(entropyLotteryGraphic.panel);
	    }

        private static final long serialVersionUID = 1L;
    }
    
    protected JPanel panel;
}


