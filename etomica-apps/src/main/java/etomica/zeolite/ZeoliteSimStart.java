/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.zeolite;

import etomica.action.IntegratorActionAdapter;
import etomica.atom.iterator.AtomIteratorLeafFilteredType;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPump;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.meter.MeterEnergy;
import etomica.graphics.DisplayPlot;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMD;
import etomica.species.SpeciesSpheresMono;

public class ZeoliteSimStart extends IntegratorActionAdapter{

    public ZeoliteSimStart(ZeoliteSimulation sim,zeoliteSimGraphic graphic) {
        this(sim.integrator,false,sim.getInterval());
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
        	//Turn off BoxDisplay
        	//graphic.remove(graphic.display);
        	
        	
        	//Set Max Steps Based on SimulationParameters.xls
        	sim.getController2().setMaxSteps(3270000/2);
        	sim.getController2().setSleepPeriod(0);
        	((IntegratorMD)integrator).setThermostatInterval(327000);
        	//Keeping another graphic of the total energy drift
        	MeterEnergy eMeter = new MeterEnergy(((IntegratorBox)integrator).getPotentialMaster(), sim.box);
        	AccumulatorHistory energyHistory = new AccumulatorHistory(new HistoryCollapsingDiscard());
        	energyHistory.getHistory().setHistoryLength(sim.getInterval()*500);
        	//AccumulatorAverage enAcc = new AccumulatorAverage(sim);
        	//enAcc.setPushInterval(20);
        	//DataFork enFork = new DataFork(new DataSink[]{energyHistory, enAcc});
        	DataPump energyPump = new DataPump(eMeter, energyHistory);
        	IntegratorListenerAction pumpListener = new IntegratorListenerAction(energyPump);
        	pumpListener.setInterval(10);
            sim.integrator.getEventManager().addListener(pumpListener);
        	energyHistory.setPushInterval(50);
            //XXX BORK!  This class needs to be rewritten
        	//sim.register(eMeter,energyPump);
        	
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
            MSDCoordWriter coordWriter = new MSDCoordWriter(sim.integrator, sim.box, filename, interval);
            coordWriter.setIterator(new AtomIteratorLeafFilteredType(sim.box, sp.getLeafType()));
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
    private SpeciesSpheresMono sp;
    private zeoliteSimGraphic graphic;
}
