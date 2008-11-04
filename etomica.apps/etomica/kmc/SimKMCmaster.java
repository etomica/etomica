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
        double temp = Double.parseDouble(args[0]);
        int steps = Integer.parseInt(args[1]);
        int totalSearch = Integer.parseInt(args[2]);
        final String APP_NAME = "SimKMCmaster";

        final SimKMCLJadatom sim = new SimKMCLJadatom();
        IVector vect = sim.getSpace().makeVector();
        vect.setX(0, 3.5);
        vect.setX(1, 0.0);
        vect.setX(2, 0.0);
        
        
        sim.setMovableAtoms(2.0, vect);
        
        sim.initializeConfiguration("initialStart");
        
        sim.integratorKMCCluster(temp, steps, totalSearch);
        sim.integratorKMCCluster.setInitialStateConditions(-539.543484823175, 3.1145942027562522E72);
        sim.getController().actionPerformed();
    }
}
