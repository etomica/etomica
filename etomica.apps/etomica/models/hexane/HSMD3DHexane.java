package etomica.models.hexane;

import etomica.data.meter.MeterPressureHard;
import etomica.integrator.IntegratorHard;
import etomica.nbr.list.PotentialMasterList;
import etomica.phase.Phase;
import etomica.potential.Potential2;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.Default;

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
        super(space, true, new int[]{1,4,4,21,1,1}, new Default());
        
        PotentialMasterList potentialMaster = new PotentialMasterList(space);
        
        double neighborRangeFac = 1.6;
        defaults.makeLJDefaults();
        
        defaults.boxSize = 14.4573*Math.pow(numAtoms/2000.0,1.0/3.0);
        potentialMaster.setCellRange(1);
        potentialMaster.setRange(neighborRangeFac*defaults.atomSize);
        
        
        
        
        
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        int numAtoms = 4000;
        if (args.length > 0) {
            numAtoms = Integer.valueOf(args[0]).intValue();
        }
        HSMD3DHexane sim = new HSMD3DHexane(Space3D.getInstance(), numAtoms);

        MeterPressureHard pMeter = new MeterPressureHard(sim.getSpace());
        pMeter.setIntegrator(sim.integrator);
        
        sim.getController().actionPerformed();
        
        double Z = pMeter.getDataAsScalar()*sim.phase.volume()/(sim.phase.moleculeCount()*sim.integrator.getTemperature());
        System.out.println("Z="+Z);
        
        // compressibility factor for this system should be 5.22
        if (Math.abs(Z-5.22) > 0.03) {
            System.exit(1);
        }
    }
}
