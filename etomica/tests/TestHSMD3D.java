package etomica.tests;

import etomica.ConfigurationFile;
import etomica.Default;
import etomica.IntegratorHard;
import etomica.P2HardSphere;
import etomica.Phase;
import etomica.Potential2;
import etomica.Simulation;
import etomica.Space;
import etomica.Space3D;
import etomica.Species;
import etomica.SpeciesSpheresMono;
import etomica.action.activity.ActivityIntegrate;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.NeighborCriterionSimple;
import etomica.nbr.PotentialMasterNbr;
import etomica.simulations.HSMD2D;

/**
 * Simple hard-sphere molecular dynamics simulation in 3D.
 *
 * @author David Kofke
 */
 
public class TestHSMD3D extends Simulation {
    
    public IntegratorHard integrator;
    public SpeciesSpheresMono species, species2;
    public Phase phase;
    public Potential2 potential;

    public TestHSMD3D(Space space, int numAtoms) {
        super(space, new PotentialMasterNbr(space));
        
        double neighborRangeFac = 1.6;
        Default.makeLJDefaults();
        Default.ATOM_SIZE = 1.0;
        // makes eta = 0.35
        Default.BOX_SIZE = 14.4573*Math.pow((numAtoms/2000.0),1.0/3.0);
        ((PotentialMasterNbr)potentialMaster).setNCells((int)(Default.BOX_SIZE/neighborRangeFac));
        integrator = new IntegratorHard(potentialMaster);
        integrator.addIntervalListener(((PotentialMasterNbr)potentialMaster).getNeighborManager());
        integrator.setTimeStep(0.01);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        activityIntegrate.setMaxSteps(200000/numAtoms);
        species = new SpeciesSpheresMono(this);
        species2 = new SpeciesSpheresMono(this);
        species.setNMolecules(numAtoms);
        species2.setNMolecules(numAtoms/100);
        phase = new Phase(space);
        potential = new P2HardSphere(space);

        NeighborCriterion criterion = new NeighborCriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());
        ((PotentialMasterNbr)potentialMaster).setSpecies(potential,new Species[]{species,species},criterion);
        criterion = new NeighborCriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());
        ((PotentialMasterNbr)potentialMaster).setSpecies(potential,new Species[]{species2,species2},criterion);
        criterion = new NeighborCriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());
        ((PotentialMasterNbr)potentialMaster).setSpecies(potential,new Species[]{species,species2},criterion);
        
        phase.setConfiguration(new ConfigurationFile(space,Integer.toString(numAtoms)));
        phase.speciesMaster.addSpecies(species);
        phase.speciesMaster.addSpecies(species2);
        integrator.addPhase(phase);

//        WriteConfiguration writeConfig = new WriteConfiguration("foo",phase,1);
//        integrator.addIntervalListener(writeConfig);
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        int numAtoms = 2000;
        if (args.length > 0) {
            numAtoms = Integer.valueOf(args[0]).intValue();
        }
        TestHSMD3D sim = new TestHSMD3D(new Space3D(), numAtoms);
        sim.getController().actionPerformed();
    }
}
