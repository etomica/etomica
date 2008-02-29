package etomica.paracetamol;

import etomica.action.WriteConfiguration;
import etomica.lattice.crystal.PrimitiveOrthorhombic;
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
 * MC simulation of Paracetamol molecules in Form II (Orthorhombic) with 
 *  tabulation of the collective-coordinate S-matrix.
 * No graphic display
 * 
 * Orthorhombic Crystal
 * 
 * @author Tai Tan
 *
 */
public class MCParacetamolOrthorhombicDLMULTISimCalcS extends Simulation {

	private static final long serialVersionUID = 1L;
	
	MCParacetamolOrthorhombicDLMULTISimCalcS(Space _space) {
		super(_space);
	}

	public static void main(String[] args) {
    	
    	int numMolecules = 32;
    	double temperature = Kelvin.UNIT.toSim(10);
    	long simSteps = 1000;
    	int simType = 2;
    	
        String filename = "SimCalcS_Orthorhombic_"+ Kelvin.UNIT.fromSim(temperature)+"K";
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

        System.out.println("Running "+ " Sim Calculate-S Orthorhombic Paracetamol simulation");
        System.out.println("Type of simulation: "+simType);
        System.out.println(numMolecules + " molecules " +" and temperature "+ Kelvin.UNIT.fromSim(temperature) +"K");
        System.out.println(simSteps+ " steps");
        System.out.println("output data to " + filename);
        
    	MCParacetamolOrthorhombicDLMULTI sim = 
        	new MCParacetamolOrthorhombicDLMULTI(Space.getInstance(3), numMolecules, temperature, simType, new int[] {1,2,2});
        
        sim.actionIntegrate.setMaxSteps(simSteps);
        PrimitiveOrthorhombic primitive = sim.primitive;
       
        //Set up Normal-Mode Meter
        MeterNormalMode meterNormalMode = new MeterNormalMode();
        meterNormalMode.setCoordinateDefinition(sim.coordDef);
        WaveVectorFactorySimple waveVectorFactory = new WaveVectorFactorySimple(primitive, sim.getSpace());
       
        meterNormalMode.setWaveVectorFactory(waveVectorFactory);
        meterNormalMode.setBox(sim.box);
       
        sim.integrator.addIntervalAction(meterNormalMode);
        sim.integrator.setActionInterval(meterNormalMode, 300);
    
        //Write S-Vectors
        WriteS sWriter = new WriteS(sim.getSpace());
        sWriter.setFilename(filename);
        sWriter.setMeter(meterNormalMode);
        sWriter.setWaveVectorFactory(waveVectorFactory);
        sWriter.setTemperature(temperature);
        sWriter.setOverwrite(true);
       
        sim.integrator.addIntervalAction(sWriter);
        sim.integrator.setActionInterval(sWriter, 1000);
       
        sim.getController().actionPerformed();
        
        WriteConfiguration writeConfig = new WriteConfiguration(sim.getSpace());
        writeConfig.setConfName("FinalCoord_SimCalcS_Orthorhombic_"+Kelvin.UNIT.fromSim(temperature)+"K");
        writeConfig.setBox(sim.box);
        writeConfig.setDoApplyPBC(false);
        writeConfig.actionPerformed();
        
    }//end of main

}//end of class