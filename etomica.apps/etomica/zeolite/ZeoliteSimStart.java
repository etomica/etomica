package etomica.zeolite;

import etomica.action.IntegratorActionAdapter;
import etomica.atom.iterator.AtomIteratorMolecule;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPump;
import etomica.data.meter.MeterEnergy;
import etomica.graphics.DisplayPlot;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorPhase;
import etomica.species.Species;
import etomica.util.HistoryCollapsing;

public class ZeoliteSimStart extends IntegratorActionAdapter{

    public ZeoliteSimStart(ZeoliteSimulation sim,zeoliteSimGraphic graphic) {
        this(sim.integrator,sim.getDefaults().ignoreOverlap,sim.getInterval());
        this.sim = sim;
        this.graphic = graphic;
	}
    
    public ZeoliteSimStart(IntegratorMD integrator, boolean ignoreOverlap,int interval) {
        super(integrator);
        this.ignoreOverlap = ignoreOverlap;
        this.interval = interval;
        
    }

    public void actionPerformed() {
        if(integrator != null) {
        	//Turn off PhaseDisplay
        	//graphic.remove(graphic.display);
        	
        	
        	//Set Max Steps Based on SimulationParameters.xls
        	sim.activityIntegrate.setMaxSteps(3270000/2);
        	sim.activityIntegrate.setDoSleep(false);
        	((IntegratorMD)integrator).setThermostatInterval(327000);
        	//Keeping another graphic of the total energy drift
        	MeterEnergy eMeter = new MeterEnergy(((IntegratorPhase)integrator).getPotential());
        	eMeter.setPhase(sim.phase);
        	AccumulatorHistory energyHistory = new AccumulatorHistory(new HistoryCollapsing());
        	energyHistory.getHistory().setHistoryLength(sim.getInterval()*500);
        	//AccumulatorAverage enAcc = new AccumulatorAverage(sim);
        	//enAcc.setPushInterval(20);
        	//DataFork enFork = new DataFork(new DataSink[]{energyHistory, enAcc});
        	DataPump energyPump = new DataPump(eMeter, energyHistory);
            sim.integrator.addIntervalAction(energyPump);
            sim.integrator.setActionInterval(energyPump, 10);
        	energyHistory.setPushInterval(50);
        	sim.register(eMeter,energyPump);
        	
        	DisplayPlot hisPlot = new DisplayPlot();
        	hisPlot.setLabel("Energy");
        	energyHistory.setDataSink(hisPlot.getDataSet().makeDataSink());
        	hisPlot.setDoLegend(true);
        	graphic.add(hisPlot);
        	
        	interval = 327;
        		
        	
        	filename = sim.getFileName();
        	sp = sim.getSpeciesRMS();
        	System.out.println(filename);
        	//sim.integrator.setTimeStep(0.00);
        	MSDCoordWriter coordWriter = new MSDCoordWriter(sim.getSpace(), filename);
        	coordWriter.setPhase(sim.phase);
            coordWriter.setIterator(new AtomIteratorMolecule(sp));
            coordWriter.setIntegrator(sim.integrator);
            coordWriter.setWriteInterval(interval);
            coordWriter.openFile();
            System.out.println("created MSDCoordWriter");
        	
        }
    }

    /**
     * @return Returns the ignoreOverlap.
     */
    public boolean isIgnoreOverlap() {
        return ignoreOverlap;
    }

    /**
     * @param ignoreOverlap The ignoreOverlap to set.
     */
    public void setIgnoreOverlap(boolean ignoreOverlap) {
        this.ignoreOverlap = ignoreOverlap;
    }

    private static final long serialVersionUID = 1L;
    private boolean ignoreOverlap;
    private int interval;
    private ZeoliteSimulation sim;
    private String filename;
    private Species[] sp;
    private zeoliteSimGraphic graphic;
}
