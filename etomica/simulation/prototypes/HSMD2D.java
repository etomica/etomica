package etomica.simulation.prototypes;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorHard;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
import etomica.simulation.Simulation;
import etomica.space2d.Space2D;
import etomica.species.SpeciesSpheresMono;

/**
 * Simple hard-sphere molecular dynamics simulation in 2D.
 *
 * @author David Kofke
 */
 
public class HSMD2D extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public IntegratorHard integrator;
    public SpeciesSpheresMono species1, species2;
    public Box box;
    public Potential2 potential11;
    public Potential2 potential12;
    public Potential2 potential22;

    public HSMD2D() {
        super(Space2D.getInstance(), true);
        PotentialMasterList potentialMaster = new PotentialMasterList(this);
//        super(space, new PotentialMaster(space));//,IteratorFactoryCell.instance));
        double sigma = 0.38;

        double neighborRangeFac = 1.6;
        potentialMaster.setRange(neighborRangeFac*sigma);

        integrator = new IntegratorHard(this, potentialMaster);
        integrator.setIsothermal(false);
        integrator.setTimeStep(0.01);

        potentialMaster.setRange(sigma*1.6);

        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);
        species1 = new SpeciesSpheresMono(this);
	    species2 = new SpeciesSpheresMono(this);
	    AtomTypeSphere leafType1 = (AtomTypeSphere)species1.getLeafType();
        AtomTypeSphere leafType2 = (AtomTypeSphere)species2.getLeafType();
        leafType1.setDiameter(sigma);
        leafType2.setDiameter(sigma);
        getSpeciesManager().addSpecies(species1);
        getSpeciesManager().addSpecies(species2);
        potential11 = new P2HardSphere(space, sigma, false);
        potential12 = new P2HardSphere(space, sigma, false);
        potential22 = new P2HardSphere(space, sigma, false);
        
        potentialMaster.addPotential(potential11,new AtomType[]{leafType1, leafType1});

        potentialMaster.addPotential(potential12,new AtomType[]{leafType2, leafType2});

        potentialMaster.addPotential(potential22,new AtomType[]{leafType1, leafType2});

        box = new Box(this);
        addBox(box);
        box.setNMolecules(species1, 512);
        box.setNMolecules(species2, 5);
        NeighborListManager nbrManager = potentialMaster.getNeighborManager(box);
        integrator.addIntervalAction(nbrManager);
        integrator.addNonintervalListener(nbrManager);
        new ConfigurationLattice(new LatticeOrthorhombicHexagonal()).initializeCoordinates(box);
        integrator.setBox(box);
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        HSMD2D sim = new HSMD2D();
		sim.getController().actionPerformed();
    }
    
}