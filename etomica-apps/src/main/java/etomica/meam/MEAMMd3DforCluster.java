/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.meam;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataLogger;
import etomica.data.DataPump;
import etomica.data.DataTableWriter;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorListenerAction;

/**
 * Molecular-Dynamics Simulation Using the Modified Embedded-Atom Method 
 * (MEAM) Potential.  
 * 
 * The MEAM potential is intended for use with metallic and covalently-bonded
 * solid systems.
 * 
 * The MEAM potential for an atom is built using terms describing parts of the
 * relationships between the atom and each of its neighbors, the number of which 
 * is determined by a cutoff and/or screening function.  Each type of pair-
 * wise term is summed over all the neighbors, and then used in expressions 
 * describing the embedding energy and the repulsive energy of the atom.   
 * Effectively, the MEAM potential is a many-body potential.  
 * 
 * This class was adapted from LjMd3D.java by K.R. Schadel and A. Schultz in July 
 * 2005.  Intitially, it employed a version of the embedded-atom method potential, 
 * and was later adapted in February 2006 to use the modified embedded-atom method
 * potential.
 */
 
public class MEAMMd3DforCluster {

    public static void main(String[] args) {
    	MEAMMd3D sim = new MEAMMd3D();
    	double temperature = sim.integrator.getTemperature();
    	if (args.length > 0) {
    		temperature = Double.parseDouble(args[0]);
    		sim.integrator.setTemperature(temperature);
    	}
    	
    	sim.activityIntegrate.setMaxSteps(3000);
    	sim.getController().actionPerformed();
    	
    	MeterPotentialEnergy energyMeter = new MeterPotentialEnergy(sim.potentialMaster);
    	MeterKineticEnergy kineticMeter = new MeterKineticEnergy();
    	
    	energyMeter.setBox(sim.box);
    	kineticMeter.setBox(sim.box);
        
        AccumulatorAverage accumulatorAveragePE = new AccumulatorAverageFixed(50);
    	AccumulatorAverage accumulatorAverageKE = new AccumulatorAverageFixed(50);
    	
    	DataPump energyManager = new DataPump(energyMeter,accumulatorAveragePE);   	
    	DataPump kineticManager = new DataPump(kineticMeter, accumulatorAverageKE);
    	
    	DataLogger dataPELogger = new DataLogger();
    	dataPELogger.setFileName("PEData_"+temperature+"K");
    	dataPELogger.setAppending(true);
    	dataPELogger.setDataSink(new DataTableWriter());
    	
    	DataLogger dataKELogger = new DataLogger();
    	dataKELogger.setFileName("KEData_"+temperature+"K");
    	dataKELogger.setAppending(true);
    	dataKELogger.setDataSink(new DataTableWriter());
    	
    	sim.getController().getEventManager().addListener(dataPELogger);
    	sim.getController().getEventManager().addListener(dataKELogger);
    	
    	accumulatorAveragePE.addDataSink(dataPELogger, new StatType[]{accumulatorAveragePE.MOST_RECENT});
    	accumulatorAverageKE.addDataSink(dataKELogger, new StatType[]{accumulatorAverageKE.MOST_RECENT});
        
    	accumulatorAveragePE.setPushInterval(1);
    	accumulatorAverageKE.setPushInterval(1);

        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(energyManager));
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(kineticManager));
        
        sim.getController().addAction(sim.activityIntegrate);
        sim.activityIntegrate.setMaxSteps(1000000);
        sim.getController().actionPerformed();
       
    	double cvPE = ((DataDouble)((DataGroup)accumulatorAveragePE.getData()).
    			getData(accumulatorAveragePE.STANDARD_DEVIATION.index)).x;
    	double systemTemp = sim.integrator.getTemperature();
    	cvPE /= systemTemp;
    	cvPE *= cvPE/sim.box.getMoleculeList().getMoleculeCount();
    	
    	double cvKE = ((DataDouble)((DataGroup)accumulatorAverageKE.getData()).
    			getData(accumulatorAverageKE.STANDARD_DEVIATION.index)).x;
      	cvKE /= systemTemp;
    	cvKE *= cvKE/sim.box.getMoleculeList().getMoleculeCount();
    
    	System.out.println("PE contribution to Cv/k is  "+ cvPE + " (simulation units)");
    	System.out.println("KE contribution to Cv/k is  "+ cvKE + " (simulation units)");
    	System.out.println("Cv/k is  "+ (cvKE+cvPE) + " (simulation units)");
    }
    
}
