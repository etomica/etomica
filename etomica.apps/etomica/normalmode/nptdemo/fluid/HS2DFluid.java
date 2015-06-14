package etomica.normalmode.nptdemo.fluid;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorHard;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.CoordinateDefinition;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2Spherical;
import etomica.simulation.Simulation;
import etomica.space2d.Space2D;
import etomica.species.SpeciesSpheresMono;

/**
 * Simple hard-sphere molecular dynamics simulation in 2D.
 *
 * @author David Kofke
 */
 
public class HS2DFluid extends Simulation {
    
    public ActivityIntegrate ai;
    public IntegratorHard integrator;
    public SpeciesSpheresMono species1;
    public IBox box;
    public Potential2Spherical potential;
    public PotentialMasterList potentialMaster;
    public CoordinateDefinition coordinateDefinition;
    public double pressure;

    public HS2DFluid() {
        super(Space2D.getInstance());
        potentialMaster = new PotentialMasterList(this, space);
//        super(space, new PotentialMaster(space));//,IteratorFactoryCell.instance));
        double sigma = 1;

        double neighborRangeFac = 1.6;
        potentialMaster.setRange(neighborRangeFac*sigma);

        integrator = new IntegratorHard(this, potentialMaster, space);
        integrator.setIsothermal(false);
        integrator.setTimeStep(0.01);

        ai = new ActivityIntegrate(integrator);
        getController().addAction(ai);
        species1 = new SpeciesSpheresMono(this, space);
        species1.setIsDynamic(true);
	    IAtomType leafType1 = species1.getLeafType();
        addSpecies(species1);
        potential = new P2HardSphere(space, sigma, false);
        
        potentialMaster.addPotential(potential,new IAtomType[]{leafType1, leafType1});

        box = new Box(space);
        addBox(box);
        double rho = 0.2;
        box.setNMolecules(species1, 50);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(rho);
        inflater.actionPerformed();
        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));
        
        new ConfigurationLattice(new LatticeOrthorhombicHexagonal(space), space).initializeCoordinates(box);
        integrator.setBox(box);
    }

}
