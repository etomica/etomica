/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.paracetamol;

import etomica.action.PDBWriter;
import etomica.action.WriteConfiguration;
import etomica.data.DataLogger;
import etomica.data.DataPump;
import etomica.data.DataTableWriter;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorListenerAction;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.units.Kelvin;

/**
 * 
 * Three-dimensional soft-sphere MC simulation for paracetamol molecule
 *  using the potential model from DL_MULTI package
 *  with inclusion of multipole ewald-summation
 *  
 *  To get the lowest minimum cofiguration with MCMoveHarmonicStep
 * 
 * Orthorhombic Crystal
 * 
 * @author Tai Tan
 *
 */
public class MCParacetamolOrthorhombicDLMULTIEquilibration extends Simulation {

	private static final long serialVersionUID = 1L;

	public MCParacetamolOrthorhombicDLMULTIEquilibration(Space _space) {
		super(_space);
	}

    public static void main(String[] args) {
    	
    	int numMolecules = 192;
    	double temperature = Kelvin.UNIT.toSim(10);
    	long simSteps = 100;
    	int simType = 1;
    	
        String filename = "Paracetamol_Orthorhombic_"+ Kelvin.UNIT.fromSim(temperature)+"K";
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

        System.out.println("Running "+ "Orthorhombic Paracetamol Equilibration simulation");
        System.out.println("Type of simulation: "+simType);
        System.out.println(numMolecules + " molecules " +" and temperature "+ Kelvin.UNIT.fromSim(temperature) +"K");
        System.out.println(simSteps+ " steps");
        System.out.println("output data to " + filename);
       
    	Space sp = Space.getInstance(3);
    	
    	etomica.paracetamol.MCParacetamolOrthorhombicDLMULTI sim = 
        	new MCParacetamolOrthorhombicDLMULTI(sp, numMolecules, temperature, simType, new int[] {2,3,4});
       

       sim.actionIntegrate.setMaxSteps(simSteps);
       MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster);
       meterPE.setBox(sim.box);
       
       DataLogger dataLoggerPE = new DataLogger();
       dataLoggerPE.setWriteInterval(5);
       dataLoggerPE.setFileName("Paracetamol_Orthorhombic_"+Kelvin.UNIT.fromSim(temperature)+"K");
       dataLoggerPE.setDataSink(new DataTableWriter());
       dataLoggerPE.setAppending(true);
       dataLoggerPE.setCloseFileEachTime(true);
       
       DataPump PEpump = new DataPump(meterPE, dataLoggerPE);
       sim.getController().getEventManager().addListener(dataLoggerPE);
       
       IntegratorListenerAction PEpumpListener = new IntegratorListenerAction(PEpump);
       PEpumpListener.setInterval(1);
       sim.integrator.getEventManager().addListener(PEpumpListener);
       
       
       sim.getController().actionPerformed();
       
       WriteConfiguration writeConfig = new WriteConfiguration(sp);
       writeConfig.setConfName("FinalCoord_Paracetamol_Orthorhombic_"+Kelvin.UNIT.fromSim(temperature)+"K");
       writeConfig.setBox(sim.box);
       writeConfig.setDoApplyPBC(false);
       writeConfig.actionPerformed();
        
       PDBWriter pdbWriter = new PDBWriter(sim.box);
       pdbWriter.setFileName("GraphicView_Orthorhombic_"+ Kelvin.UNIT.fromSim(temperature)+"K.pdb");
       pdbWriter.actionPerformed();
        
    }//end of main

}//end of class
