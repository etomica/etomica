package etomica.simulation.prototypes;
import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.atom.AtomTypeSphere;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorGEMC;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.lattice.SpaceLattice;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.species.SpeciesSpheresRotating;

/**
 * Simple Gibbs-ensemble Monte Carlo simulation of rotating molecules.
 */
//in present form uses just a LJ potential, so orientation is irrelevant
public class GEMCWithRotation extends Simulation {
    
    private static final long serialVersionUID = 1L;

    
    public GEMCWithRotation(Space space) {
        super(space, false);
        double sigma = 1.2;
        IPotentialMaster potentialMaster = new PotentialMaster(space);
        integrator = new IntegratorGEMC(getRandom());
        integrator.setEventInterval(400);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        activityIntegrate.setSleepPeriod(1);
	    
	    species = new SpeciesSpheresRotating(this);
        getSpeciesManager().addSpecies(species);
        ((AtomTypeSphere)species.getMoleculeType().getChildTypes()[0]).setDiameter(sigma);

	    box1 = new Box(this, space);
        addBox(box1);
        box1.setNMolecules(species, 200);
        
	    IntegratorMC integratorMC1 = new IntegratorMC(this, potentialMaster);
        integratorMC1.setBox(box1);
        integratorMC1.setTemperature(0.420);
        MCMoveManager moveManager = integratorMC1.getMoveManager();
        moveManager.addMCMove(new MCMoveRotate(potentialMaster, getRandom()));
        moveManager.addMCMove(new MCMoveAtom(this, potentialMaster));
        integrator.addIntegrator(integratorMC1);
        

	    box2 = new Box(this, space);
        addBox(box2);
        box2.setNMolecules(species, 200);
        IntegratorMC integratorMC2 = new IntegratorMC(this, potentialMaster);
        integratorMC2.setBox(box2);
        integratorMC2.setTemperature(0.420);
        moveManager = integratorMC2.getMoveManager();
        moveManager.addMCMove(new MCMoveRotate(potentialMaster, getRandom()));
        moveManager.addMCMove(new MCMoveAtom(this, potentialMaster));
        // GEMC integrator adds volume and molecule exchange moves once
        // it has 2 integrators
        integrator.addIntegrator(integratorMC2);
        
        SpaceLattice lattice;
        if (space.D() == 2) {
            lattice = new LatticeOrthorhombicHexagonal();
        }
        else {
            lattice = new LatticeCubicFcc();
        }
        ConfigurationLattice config = new ConfigurationLattice(lattice, space);
        config.initializeCoordinates(box1);
        config.initializeCoordinates(box2);
            
	    potential = new P2LennardJones(space);
	    potential.setSigma(sigma);

        potentialMaster.addPotential(potential,new IAtomType[] {species.getLeafType(), species.getLeafType()});

        integratorMC1.addIntervalAction(new BoxImposePbc(box1, space));
        integratorMC2.addIntervalAction(new BoxImposePbc(box2, space));

	    box2.setDensity(0.1);
    }
    
    public IBox box1, box2;
    public IntegratorGEMC integrator;
    public SpeciesSpheresRotating species;
    public P2LennardJones potential;
}