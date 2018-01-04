/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.paracetamol;

import etomica.action.WriteConfiguration;
import etomica.lattice.crystal.PrimitiveMonoclinic;
import etomica.integrator.IntegratorListenerAction;
import etomica.normalmode.MeterNormalMode;
import etomica.normalmode.WaveVectorFactorySimple;
import etomica.normalmode.WriteS;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.units.Kelvin;

/**
 * 
 * Three-dimensional soft-sphere MC simulation for paracetamol molecule
 *  using the potential model from DL_MULTI package
 *  with inclusion of multipole ewald-summation
 * 
 * MC simulation of Paracetamol molecules in Form I (Monoclinic) with 
 *  tabulation of the collective-coordinate S-matrix.
 * No graphic display
 *  
 * Monoclinic Crystal
 * 
 * @author Tai Tan
 *
 */
public class MCParacetamolMonoclinicDLMULTISimCalcS extends Simulation{
	private static final long serialVersionUID = 1L;

	public MCParacetamolMonoclinicDLMULTISimCalcS(Space _space) {
		super(_space);
	}

    public static void main(String[] args) {
    	
    	int numMolecules = 96;
    	double temperature = Kelvin.UNIT.toSim(123);
    	long simSteps = 10000;
    	int simType = 2;
      
        String filename = "SimCalcS_Monoclinic_"+ Kelvin.UNIT.fromSim(temperature)+"K";
        if (args.length > 0) {
            filename = args[0];
        }
        if (args.length > 1) {
            simSteps = Long.parseLong(args[1]);
        }
        if (args.length > 2) {
            temperature = Kelvin.UNIT.toSim(Double.parseDouble(args[2]));
        }
        if (args.length > 3) {
            simType = Integer.parseInt(args[3]);
        }
        
        System.out.println("Running "+ "Sim Calculate-S Monoclinic Paracetamol simulation");
        System.out.println("Type of simulation: "+simType);
        System.out.println(numMolecules + " molecules " +" and temperature "+ Kelvin.UNIT.fromSim(temperature) +"K");
        System.out.println(simSteps+ " steps");
        System.out.println("output data to " + filename);
        
    	MCParacetamolMonoclinicDLMULTI sim = 
        	new MCParacetamolMonoclinicDLMULTI(Space.getInstance(3), numMolecules, temperature, simType, new int[] {2,3,4});

        sim.actionIntegrate.setMaxSteps(simSteps);
        PrimitiveMonoclinic primitive = sim.primitive;
        
        //Set up Normal-Mode Meter
        MeterNormalMode meterNormalMode = new MeterNormalMode();
        meterNormalMode.setCoordinateDefinition(sim.coordDef);
        WaveVectorFactorySimple waveVectorFactory = new WaveVectorFactorySimple(primitive, sim.getSpace());
       
        meterNormalMode.setWaveVectorFactory(waveVectorFactory);
        meterNormalMode.setBox(sim.box);
       
        IntegratorListenerAction meterNormalModeListener = new IntegratorListenerAction(meterNormalMode);
        meterNormalModeListener.setInterval(300);
        sim.integrator.getEventManager().addListener(meterNormalModeListener);
    
        //Write S-Vectors
        WriteS sWriter = new WriteS(sim.getSpace());
        sWriter.setFilename(filename);
        sWriter.setMeter(meterNormalMode);
        sWriter.setWaveVectorFactory(waveVectorFactory);
        sWriter.setTemperature(temperature);
        sWriter.setOverwrite(true);
       
        IntegratorListenerAction sWriterListener = new IntegratorListenerAction(sWriter);
        sWriterListener.setInterval(1000);
        sim.integrator.getEventManager().addListener(sWriterListener);
       
        sim.getController().actionPerformed();
        
        WriteConfiguration writeConfig = new WriteConfiguration(sim.getSpace());
        writeConfig.setConfName("FinalCoord_SimCalcS_Monoclinic_"+Kelvin.UNIT.fromSim(temperature)+"K");
        writeConfig.setBox(sim.box);
        writeConfig.setDoApplyPBC(false);
        writeConfig.actionPerformed();
        
    }//end of main
    
}//end of class
