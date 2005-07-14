package etomica.models.hexane;

import etomica.Default;
import etomica.Phase;
import etomica.Simulation;
import etomica.Space;
import etomica.SpeciesSpheresMono;
import etomica.data.meter.MeterPressureHard;
import etomica.integrator.IntegratorHard;
import etomica.nbr.list.PotentialMasterNbr;
import etomica.potential.Potential2;
import etomica.space3d.Space3D;

/**
 * Hard sphere simulation of hexane molecule behavior.
 * 
 * @author nancycribbin
 */

public class HSMD3DHexane extends Simulation {

    public IntegratorHard integrator;
    public SpeciesSpheresMono species, species2;
    public Phase phase;
    public Potential2 potential;
    
    public HSMD3DHexane(Space space, int numAtoms){
        // use custom bit lengths to allow for more "molecules"        
        super(space, true, new PotentialMasterNbr(space), new int[]{1,4,4,21,1,1});
        
        double neighborRangeFac = 1.6;
        Default.makeLJDefaults();
        
        Default.BOX_SIZE = 14.4573*Math.pow(numAtoms/2000.0,1.0/3.0);
        ((PotentialMasterNbr)potentialMaster).setNCells((int)(Default.BOX_SIZE/neighborRangeFac));
        
        
        
        
        
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        int numAtoms = 4000;
        if (args.length > 0) {
            numAtoms = Integer.valueOf(args[0]).intValue();
        }
        HSMD3DHexane sim = new HSMD3DHexane(new Space3D(), numAtoms);

        MeterPressureHard pMeter = new MeterPressureHard(sim.space,sim.integrator);
        pMeter.setPhase(sim.phase);
        
        sim.getController().actionPerformed();
        
        double Z = pMeter.getDataAsScalar()*sim.phase.volume()/(sim.phase.moleculeCount()*sim.integrator.getTemperature());
        System.out.println("Z="+Z);
        
        // compressibility factor for this system should be 5.22
        if (Math.abs(Z-5.22) > 0.03) {
            System.exit(1);
        }
    }
}
