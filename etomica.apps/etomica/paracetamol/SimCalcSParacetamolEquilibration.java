package etomica.paracetamol;

import etomica.action.WriteConfiguration;
import etomica.data.DataLogger;
import etomica.data.DataPump;
import etomica.data.DataTableWriter;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.units.Kelvin;

/**
 * MC simulation of 3D Paracetamol molecules in Form II (Orthorhombic) with tabulation of the
 * collective-coordinate S-matrix. No graphic display of simulation.
 */
public class SimCalcSParacetamolEquilibration extends Simulation {

    

    /**
     * @param args
     */
    public static void main(String[] args) {

        int nA = 192;
        double temperature = Kelvin.UNIT.toSim(10);
  
        long simSteps = 1000000;

        // parse arguments

        if (args.length > 1) {
            simSteps = Long.parseLong(args[1]);
        }
        if (args.length > 2) {
            temperature = Kelvin.UNIT.toSim(Double.parseDouble(args[2]));
        }
        String filename = "Paracetamol_FormII_"+ Kelvin.UNIT.fromSim(temperature)+"K";
        if (args.length > 0) {
            filename = args[0];
        }

        System.out.println("Running "+ " Orthorhombic Paracetamol simulation");
        System.out.println(nA + " atoms " +" and temperature "+ Kelvin.UNIT.fromSim(temperature));
        System.out.println(simSteps+ " steps");
        System.out.println("output data to " + filename);

        // construct simulation
        SimCalcSParacetamol sim = new SimCalcSParacetamol(Space.getInstance(3), nA, temperature);
        
        sim.activityIntegrate.setMaxSteps(simSteps);
        
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.integrator.getPotential());
        meterPE.setBox(sim.box);
        
        DataLogger dataLoggerPE = new DataLogger();
        dataLoggerPE.setWriteInterval(1);
        dataLoggerPE.setFileName("Paracetamol_Form_II_@_" + Kelvin.UNIT.fromSim(temperature));
        dataLoggerPE.setAppending(true);
        dataLoggerPE.setDataSink(new DataTableWriter());
        
        DataPump PEpump = new DataPump(meterPE, dataLoggerPE);
        sim.getController().getEventManager().addListener(dataLoggerPE);

        sim.integrator.addIntervalAction(PEpump);
        sim.integrator.setActionInterval(PEpump, 1000);
        
        
        sim.getController().actionPerformed();
        WriteConfiguration writeConfig = new WriteConfiguration();
        writeConfig.setConfName("Coord_Paracetamol_FormII_"+Kelvin.UNIT.fromSim(temperature)+"_K");
        writeConfig.setBox(sim.box);
        writeConfig.setDoApplyPBC(false);
        writeConfig.actionPerformed();
        
    }
    
    public PotentialMaster potentialMaster;
    private static final long serialVersionUID = 1L;
}