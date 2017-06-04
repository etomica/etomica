package etomica.normalmode.nptdemo;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.integrator.IntegratorHard;
import etomica.lattice.crystal.BasisOrthorhombicHexagonal;
import etomica.lattice.crystal.PrimitiveOrthorhombicHexagonal;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinitionLeaf;
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
 
public class HSNPT2DSim extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public ActivityIntegrate ai;
    public IntegratorHard integrator;
    public SpeciesSpheresMono species1;
    public Box box;
    public Potential2 potential;
    public PotentialMasterList potentialMaster;
    public CoordinateDefinition coordinateDefinition;
    public double pressure;

    public HSNPT2DSim() {
        super(Space2D.getInstance());
        potentialMaster = new PotentialMasterList(this, space);
//        super(space, new PotentialMaster(space));//,IteratorFactoryCell.instance));
        double sigma = 1;

        double neighborRangeFac = 1.4;
        potentialMaster.setRange(neighborRangeFac*sigma);

        integrator = new IntegratorHard(this, potentialMaster, space);
        integrator.setIsothermal(false);
        integrator.setTimeStep(0.05);

        potentialMaster.setRange(sigma*1.6);

        ai = new ActivityIntegrate(integrator);
        getController().addAction(ai);
        species1 = new SpeciesSpheresMono(this, space);
        species1.setIsDynamic(true);
        AtomType leafType1 = species1.getLeafType();
        addSpecies(species1);
        potential = new P2HardSphere(space, sigma, false);

        potentialMaster.addPotential(potential, new AtomType[]{leafType1, leafType1});

        box = new Box(space);
        addBox(box);
        int nx = 10;
        int ny = 6;
        double rho = 1.0;
        box.setNMolecules(species1, nx*ny*2);
        double bx = 1;
        double by = Math.sqrt(3);
        double v1 = Math.sqrt(3);
        double v2 = 2/rho;
        bx *= nx*Math.sqrt(v2/v1);
        by *= ny*Math.sqrt(v2/v1);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{bx,by}));
        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));
        
        coordinateDefinition = new CoordinateDefinitionLeaf(box, new PrimitiveOrthorhombicHexagonal(space, bx/nx), new BasisOrthorhombicHexagonal(), space);
        coordinateDefinition.initializeCoordinates(new int[]{nx, ny});
        integrator.setBox(box);
        
        potentialMaster.getNeighborManager(box).reset();
        integrator.getEventManager().removeListener(potentialMaster.getNeighborManager(box));
    }

}
