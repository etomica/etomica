package etomica.tests;

import etomica.ConfigurationFile;
import etomica.Default;
import etomica.Phase;
import etomica.Simulation;
import etomica.Space;
import etomica.Species;
import etomica.SpeciesSpheresMono;
import etomica.action.activity.ActivityIntegrate;
import etomica.data.meter.MeterPressureHard;
import etomica.integrator.IntegratorHard;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.NeighborCriterionSimple;
import etomica.nbratom.CriterionSpecies;
import etomica.nbratom.NeighborManager;
import etomica.nbratom.PotentialMasterNbr;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
import etomica.space3d.Space3D;

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
        // use custom bit lengths to allow for more "molecules"
        super(space, new PotentialMasterNbr(space), new int[] {1, 4, 4, 21, 1, 1});
        
        double neighborRangeFac = 1.6;
        Default.makeLJDefaults();
        // makes eta = 0.35
        Default.BOX_SIZE = 14.4573*Math.pow((numAtoms/2000.0),1.0/3.0);
        ((PotentialMasterNbr)potentialMaster).setNCells((int)(Default.BOX_SIZE/neighborRangeFac));
        integrator = new IntegratorHard(potentialMaster);
        NeighborManager nbrManager = ((PotentialMasterNbr)potentialMaster).getNeighborManager();
        nbrManager.setRange(Default.ATOM_SIZE*1.6);
        integrator.addIntervalListener(nbrManager);
        integrator.setTimeStep(0.01);
        integrator.setIsothermal(true);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        activityIntegrate.setMaxSteps(20000000/numAtoms);
        species = new SpeciesSpheresMono(this);
        species2 = new SpeciesSpheresMono(this);
        species.setNMolecules(numAtoms);
        species2.setNMolecules(numAtoms/100);

        potential = new P2HardSphere(space);
        NeighborCriterion nbrCriterion = new NeighborCriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());
        CriterionSpecies criterion = new CriterionSpecies(nbrCriterion);
        criterion.setIntraSpecies(true);
        potential.setCriterion(criterion);
        potentialMaster.setSpecies(potential,new Species[]{species,species});
        nbrManager.addCriterion(nbrCriterion);
        species.getFactory().getType().getNbrManagerAgent().addCriterion(nbrCriterion);

        potential = new P2HardSphere(space);
        nbrCriterion = new NeighborCriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());
        criterion = new CriterionSpecies(nbrCriterion);
        criterion.setIntraSpecies(false);
        potential.setCriterion(criterion);
        potentialMaster.setSpecies(potential,new Species[]{species,species2});
        nbrManager.addCriterion(nbrCriterion);
        species.getFactory().getType().getNbrManagerAgent().addCriterion(nbrCriterion);
        species2.getFactory().getType().getNbrManagerAgent().addCriterion(nbrCriterion);

        potential = new P2HardSphere(space);
        nbrCriterion = new NeighborCriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());
        criterion = new CriterionSpecies(nbrCriterion);
        criterion.setIntraSpecies(true);
        potential.setCriterion(criterion);
        potentialMaster.setSpecies(potential,new Species[]{species2,species2});
        nbrManager.addCriterion(nbrCriterion);
        species2.getFactory().getType().getNbrManagerAgent().addCriterion(nbrCriterion);
        
        phase = new Phase(this);
        integrator.addPhase(phase);
        phase.setConfiguration(new ConfigurationFile(space,"HSMD3D"+Integer.toString(numAtoms)));
        
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

        MeterPressureHard pMeter = new MeterPressureHard(sim.integrator);
        
        sim.getController().actionPerformed();
        
        double Z = pMeter.getDataAsScalar(sim.phase)*sim.phase.volume()/(sim.phase.moleculeCount()*sim.integrator.temperature());
        System.out.println("Z="+Z);
        
        // compressibility factor for this system should be 5.22
        if (Math.abs(Z-5.22) > 0.03) {
            System.exit(1);
        }
    }
}
