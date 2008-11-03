package etomica.kmc;

import etomica.api.IVector;
import etomica.simulation.Simulation;
import etomica.space.ISpace;

public class SimKMCmaster extends Simulation{

    /**
     * 
     */
    private static final long serialVersionUID = 1L;

    public SimKMCmaster(ISpace space) {
        super(space);
    }

    public static void main(String[] args){
        String fileName = args[0];
        final String APP_NAME = "SimKMCmaster";

        final SimKMCMEAMadatom sim = new SimKMCMEAMadatom();
        IVector vect = sim.getSpace().makeVector();
        vect.setX(0, 9.8);
        vect.setX(1, -0.2);
        vect.setX(2, -0.2);
        
        sim.setMovableAtoms(100.0, vect);
        
        sim.setPotentialListAtoms();
        
        sim.initializeConfiguration("initialStart");
        
        sim.integratorKMCCluster(295.0, 500);
        sim.integratorKMCCluster.setInitialStateConditions(1, 2);
        
        sim.getController().actionPerformed();
    }
}
