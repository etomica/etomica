package etomica.dimer;

import etomica.simulation.Simulation;
import etomica.space.Space;


/**
 * Simulation using Henkelman's Dimer method to find a saddle point for
 * an adatom of Sn on a surface, modeled with MEAM.
 * 
 * @author msellers
 *
 */

public class EnergyMapMEAMadatomCluster extends Simulation{

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "EnergyMapMEAMadatom";

    public EnergyMapMEAMadatomCluster(Space space) {
    	super(space);
    }

    public static void main(String[] args){
        double height;
        height = Double.parseDouble(args[0]);
        String fileTail = ""+height;

    	final EnergyMapMEAMadatom sim = new EnergyMapMEAMadatom(height, fileTail);
    	
    	
    	sim.activityIntegrateMAP.setMaxSteps(1);
        sim.getController().actionPerformed();

    	
    }
}