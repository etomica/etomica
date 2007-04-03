package etomica.tests;

import etomica.action.activity.ActivityIntegrate;
import etomica.config.ConfigurationFile;
import etomica.data.meter.MeterPressureHard;
import etomica.integrator.IntegratorHard;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.phase.Phase;
import etomica.potential.P2HardSphere;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

/**
 * Simple hard-sphere molecular dynamics simulation in 3D.
 * Initial configurations at http://gordon.eng.buffalo.edu/etomica/tests/
 * @author David Kofke
 */
 
public class TestHSMD3D extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public IntegratorHard integrator;
    public SpeciesSpheresMono species, species2;
    public Phase phase;

    public TestHSMD3D(Space space, int numAtoms) {
        // use custom bit lengths to allow for more "molecules"
        super(space, true, new PotentialMasterList(space)); //, new int[] {1, 4, 4, 21, 1, 1}, new Default());
        
        double neighborRangeFac = 1.6;
        defaults.makeLJDefaults();
        // makes eta = 0.35
        defaults.boxSize = 14.4573*Math.pow((numAtoms/2000.0),1.0/3.0);
        ((PotentialMasterList)potentialMaster).setCellRange(1);
        ((PotentialMasterList)potentialMaster).setRange(neighborRangeFac*defaults.atomSize);
        integrator = new IntegratorHard(this);
        NeighborListManager nbrManager = ((PotentialMasterList)potentialMaster).getNeighborManager();
        nbrManager.setRange(defaults.atomSize*neighborRangeFac);
        nbrManager.getPbcEnforcer().setApplyToMolecules(false);
        integrator.addListener(nbrManager);
        integrator.setTimeStep(0.01);
        integrator.setIsothermal(true);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(this,integrator);
        getController().addAction(activityIntegrate);
        activityIntegrate.setMaxSteps(20000000/numAtoms);
        species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);
        getSpeciesManager().removeSpecies(species);
        species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);
        species2 = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species2);

        potentialMaster.addPotential(new P2HardSphere(this),new Species[]{species,species});

        potentialMaster.addPotential(new P2HardSphere(this),new Species[]{species,species2});

        potentialMaster.addPotential(new P2HardSphere(this),new Species[]{species2,species2});
        
        phase = new Phase(this);
        phase.getAgent(species).setNMolecules(numAtoms);
        phase.getAgent(species2).setNMolecules(numAtoms/100);
        integrator.setPhase(phase);
        ConfigurationFile config = new ConfigurationFile("HSMD3D"+Integer.toString(numAtoms));
        config.initializeCoordinates(phase);
        
//        WriteConfiguration writeConfig = new WriteConfiguration("foo",phase,1);
//        integrator.addIntervalListener(writeConfig);
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        int numAtoms = 500;
        if (args.length > 0) {
            numAtoms = Integer.valueOf(args[0]).intValue();
        }
        TestHSMD3D sim = new TestHSMD3D(Space3D.getInstance(), numAtoms);

        MeterPressureHard pMeter = new MeterPressureHard(sim.space);
        pMeter.setIntegrator(sim.integrator);
        
        sim.getController().actionPerformed();
        
        double Z = pMeter.getDataAsScalar()*sim.phase.volume()/(sim.phase.moleculeCount()*sim.integrator.getTemperature());
        System.out.println("Z="+Z);
        
        // compressibility factor for this system should be 5.22
        if (Double.isNaN(Z) || Math.abs(Z-5.22) > 0.03) {
            System.exit(1);
        }
    }
}
