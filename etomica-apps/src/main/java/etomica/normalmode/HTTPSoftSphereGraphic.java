/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.IAction;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataPumpListener;
import etomica.data.DataSplitter;
import etomica.data.DataTag;
import etomica.data.histogram.HistogramCollapsing;
import etomica.graphics.DeviceBox;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.modifier.ModifierGeneral;
import etomica.space3d.Space3D;
import etomica.util.ParameterBase;

/**
 * Applet to illustrate HTTP method
 * 
 * @author Tai Boon Tan
 */
public class HTTPSoftSphereGraphic extends SimulationGraphic {

    public HTTPSoftSphereGraphic(final HTTPSoftSphereSimFasterer sim, double otherTemperature) {
        
    	super(sim, SimulationGraphic.TABBED_PANE);
        this.sim = sim;
        
    	meter = new MeterBoltzmannHTTPFasterer(sim.potentialMaster, sim.species, sim, sim.coordinateDefinition, sim.rc);
        meter.setLatticeEnergy(sim.latticeEnergy);
        meter.setTemperature(sim.integrator.getTemperature());
        meter.setOtherTemperature(otherTemperature);
        meter.setConstraint(sim.p1Constraint);
            
        DataSplitter  histogramSplitter = new DataSplitter();
        DataPumpListener boltzmannPump = new DataPumpListener(meter, histogramSplitter, 50);
        sim.integrator.getEventManager().addListener(boltzmannPump);
            
        DisplayPlot boltzmannPlot = new DisplayPlot();
        final AccumulatorHistogram[] accumulatorHistogram = new AccumulatorHistogram[5];
        
        for(int i=0; i<5; i++){
        	accumulatorHistogram[i] = new AccumulatorHistogram(new HistogramCollapsing());
        	histogramSplitter.setDataSink(i, accumulatorHistogram[i]);
        	accumulatorHistogram[i].addDataSink(boltzmannPlot.getDataSet().makeDataSink());
        	accumulatorHistogram[i].setPushInterval(100);
        }
        
        boltzmannPlot.setLabel("Reduced Energy Histogram");
        
        boltzmannPlot.setLegend(new DataTag[]{accumulatorHistogram[0].getTag()}, "beta1*u");
        boltzmannPlot.setLegend(new DataTag[]{accumulatorHistogram[1].getTag()}, "beta2*u");
        boltzmannPlot.setLegend(new DataTag[]{accumulatorHistogram[2].getTag()}, "beta2*u'");
      
        boltzmannPlot.setLegend(new DataTag[]{accumulatorHistogram[3].getTag()}, "(beta2 - beta1)*u");
        boltzmannPlot.setLegend(new DataTag[]{accumulatorHistogram[4].getTag()}, "(beta2*u' - beta1*u");

        add(boltzmannPlot);
        
        DeviceBox t1Box = new DeviceBox(sim.getController());
        t1Box.setModifier(new ModifierGeneral(this, "t1"));
        t1Box.setLabel("Temperature 1");
        
        DeviceBox t2Box = new DeviceBox(sim.getController());
        t2Box.setModifier(new ModifierGeneral(this, "t2"));
        t2Box.setLabel("Temperature 2");
        
        add(t1Box);
        add(t2Box);

        getController().getResetAveragesButton().setPostAction(new IAction(){
        	public void actionPerformed(){
        		for (int i=0; i<5; i++){
        			accumulatorHistogram[i].reset();
        		}
        	}
        });
        
    }
    
    public void setT1(double temperature){
    	sim.integrator.setTemperature(temperature);
    	meter.setTemperature(temperature);
    }
    
    public double getT1(){
    	return sim.integrator.getTemperature();
    }
    
    public void setT2(double temperature){
    	meter.setOtherTemperature(temperature);
    }
    
    public double getT2(){
    	return meter.otherTemperature;
    }
    
    public static void main(String[] args) {
    	SimOverlapTPSSParams params = new SimOverlapTPSSParams();
    	
        int numMolecules = params.numMolecules;
        double temperature = params.temperature;
        double otherTemperature = params.otherTemperature;
        
        HTTPSoftSphereSimFasterer sim = new HTTPSoftSphereSimFasterer(Space3D.getInstance(), numMolecules, temperature);
        HTTPSoftSphereGraphic simGraphic = new HTTPSoftSphereGraphic(sim, otherTemperature);
        simGraphic.makeAndDisplayFrame("HTTP Method - Soft Sphere FCC Crystal");
    }
    
    protected HTTPSoftSphereSimFasterer sim;
    protected MeterBoltzmannHTTPFasterer meter;
    
    public static class SimOverlapTPSSParams extends ParameterBase{
        int numMolecules = 32;
        double temperature = 0.2;
        double otherTemperature = 0.4;
   }
}
