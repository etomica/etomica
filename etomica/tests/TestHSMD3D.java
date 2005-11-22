package etomica.tests;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.config.ConfigurationFile;
import etomica.data.meter.MeterPressureHard;
import etomica.integrator.IntegratorHard;
import etomica.nbr.CriterionSimple;
import etomica.nbr.CriterionSpecies;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterNbr;
import etomica.phase.Phase;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.util.Default;

/**
 * Simple hard-sphere molecular dynamics simulation in 3D.
 * Initial configurations at http://gordon.eng.buffalo.edu/etomica/tests/
 * @author David Kofke
 */
 
public class TestHSMD3D extends Simulation {
    
    public IntegratorHard integrator;
    public SpeciesSpheresMono species, species2;
    public Phase phase;
    public Potential2 potential;

    public TestHSMD3D(Space space, int numAtoms) {
        // use custom bit lengths to allow for more "molecules"
        super(space, true, new PotentialMasterNbr(space, 1.6), new int[] {1, 4, 4, 21, 1, 1}, new Default());
        
        double neighborRangeFac = 1.6;
        defaults.makeLJDefaults();
        // makes eta = 0.35
        defaults.boxSize = 14.4573*Math.pow((numAtoms/2000.0),1.0/3.0);
        ((PotentialMasterNbr)potentialMaster).setNCells((int)(defaults.boxSize/neighborRangeFac));
        integrator = new IntegratorHard(this);
        NeighborListManager nbrManager = ((PotentialMasterNbr)potentialMaster).getNeighborManager();
        nbrManager.setRange(defaults.atomSize*1.6);
        nbrManager.getPbcEnforcer().setApplyToMolecules(false);
        integrator.addListener(nbrManager);
        integrator.setTimeStep(0.01);
        integrator.setIsothermal(true);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(this,integrator);
        getController().addAction(activityIntegrate);
        activityIntegrate.setMaxSteps(20000000/numAtoms);
        species = new SpeciesSpheresMono(this);
        species2 = new SpeciesSpheresMono(this);
        species.setNMolecules(numAtoms);
        species2.setNMolecules(numAtoms/100);

        potential = new P2HardSphere(this);
        NeighborCriterion nbrCriterion = new CriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());
        CriterionSpecies criterion = new CriterionSpecies(nbrCriterion, species, species);
        potential.setCriterion(criterion);
        potentialMaster.addPotential(potential,new Species[]{species,species});
        nbrManager.addCriterion(nbrCriterion,new AtomType[]{species.getFactory().getType(),species.getFactory().getType()});

        potential = new P2HardSphere(this);
        nbrCriterion = new CriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());
        criterion = new CriterionSpecies(nbrCriterion, species, species2);
        potential.setCriterion(criterion);
        potentialMaster.addPotential(potential,new Species[]{species,species2});
        //only need to call since all the criteria are equivalent
        nbrManager.addCriterion(nbrCriterion,new AtomType[]{species.getFactory().getType(),species2.getFactory().getType()});

        potential = new P2HardSphere(this);
        nbrCriterion = new CriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());
        criterion = new CriterionSpecies(nbrCriterion, species2, species2);
        potential.setCriterion(criterion);
        potentialMaster.addPotential(potential,new Species[]{species2,species2});
        nbrManager.addCriterion(nbrCriterion,new AtomType[]{species2.getFactory().getType(),species2.getFactory().getType()});
        
        phase = new Phase(this);
        integrator.setPhase(phase);
        new ConfigurationFile(space,"HSMD3D"+Integer.toString(numAtoms)).initializeCoordinates(phase);
        
//        WriteConfiguration writeConfig = new WriteConfiguration("foo",phase,1);
//        integrator.addIntervalListener(writeConfig);
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        int numAtoms = 4000;
        if (args.length > 0) {
            numAtoms = Integer.valueOf(args[0]).intValue();
        }
        TestHSMD3D sim = new TestHSMD3D(Space3D.getInstance(), numAtoms);

        MeterPressureHard pMeter = new MeterPressureHard(sim.space,sim.integrator);
        pMeter.setPhase(sim.phase);
        
        sim.getController().actionPerformed();
        
        double Z = pMeter.getDataAsScalar()*sim.phase.volume()/(sim.phase.moleculeCount()*sim.integrator.getTemperature());
        System.out.println("Z="+Z);
        
        // compressibility factor for this system should be 5.22
        if (Double.isNaN(Z) || Math.abs(Z-5.22) > 0.03) {
            System.exit(1);
        }
    }
}
