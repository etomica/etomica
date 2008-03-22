package etomica.modules.sam;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomSet;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.atom.AtomTypeSphere;
import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConformationChain3D;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayCanvas;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.potential.P2HardMoleculeMonatomic;
import etomica.potential.P2Harmonic;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialGroup;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheres;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;
import etomica.units.Pixel;

/**
 * Self-assembled monolayer module.
 * @author Andrew Schultz
 *
 */
public class Sam extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public SpeciesSpheres species;
    public IBox box;
    public IntegratorVelocityVerlet integrator;
    public P2HardMoleculeMonatomic potentialWrapper;
    public ActivityIntegrate activityIntegrate;
    
    public Sam() {
        super(Space.getInstance(3));
        PotentialMaster potentialMaster = new PotentialMaster(space); //List(this, 2.0);

        int nCellX = 10;
        int nCellZ = 10;
        double sizeCellX = 4;
        double sizeCellZ = 4;
        int chainLength = 5;

        double surfaceSigma = 4.0;
        double sigma = 4.0;
        double bondL = 3.0;
        double epsilon = Kelvin.UNIT.toSim(100);

        //controller and integrator

	    //species and potentials
	    species = new SpeciesSpheres(this, chainLength, new ElementSimple(this), new ConformationChain3D(space, new IVector[]{space.makeVector(new double[]{0, bondL, 0})}), space);
//	    ((ElementSimple)species.getLeafType().getElement()).setMass(Dalton.UNIT.toSim(space.D() == 3 ? 131 : 40));
	    ((AtomTypeSphere)species.getLeafType()).setDiameter(sigma);
        getSpeciesManager().addSpecies(species);

        //construct box
	    box = new Box(new BoundaryRectangularSlit(this, 1, space), space);
        addBox(box);
        IVector dim = space.makeVector();
        dim.E(new double[]{sizeCellX*nCellX, chainLength*bondL, sizeCellZ*nCellZ});
        box.setDimensions(dim);
        box.setNMolecules(species, nCellX*nCellZ);

        SpeciesSpheresMono speciesSurface = new SpeciesSpheresMono(this, space);
        ((AtomTypeSphere)speciesSurface.getLeafType()).setDiameter(surfaceSigma);
        ((ElementSimple)speciesSurface.getLeafType().getElement()).setMass(Double.POSITIVE_INFINITY);
        getSpeciesManager().addSpecies(speciesSurface);

        ConfigurationSAM config = new ConfigurationSAM(this, space, species, speciesSurface);
        config.setBasisMolecules(new BasisMonatomic(space));
        config.setBasisSurface(new BasisMonatomic(space));
        config.setCellSizeX(4);
        config.setCellSizeZ(4);
        config.setNCellsX(8);
        config.setNCellsZ(8);
        config.setYOffset(4);

        config.initializeCoordinates(box);
        
        P2LennardJones p2 = new P2LennardJones(space, sigma, epsilon);
        potentialMaster.addPotential(p2, new IAtomTypeLeaf[]{species.getLeafType(), species.getLeafType()});
        potentialMaster.addPotential(p2, new IAtomTypeLeaf[]{speciesSurface.getLeafType(), species.getLeafType()});
        P2Harmonic p2Bond = new P2Harmonic(space, 10000, bondL);
        AtomsetIteratorBasisDependent bondIterator = ApiBuilder.makeAdjacentPairIterator();
        AtomsetIteratorBasisDependent nonbondedIterator = ApiBuilder.makeNonAdjacentPairIterator();
        PotentialGroup pIntra = potentialMaster.makePotentialGroup(1);
        pIntra.addPotential(p2Bond, bondIterator);
        pIntra.addPotential(p2, nonbondedIterator);
        potentialMaster.addPotential(pIntra, new ISpecies[]{species});
        ApiTether apiTether = new ApiTether(species);
        apiTether.setBox(box);
        IAtomSet polymerMolecules = box.getMoleculeList(species);
        IAtomSet surfaceMolecules = box.getMoleculeList(speciesSurface);
        int nMolecules = box.getNMolecules(species);
        for (int i=0; i<nMolecules; i++) {
            apiTether.setBondedSurfaceAtom((IMolecule)polymerMolecules.getAtom(i), (IAtomLeaf)((IMolecule)surfaceMolecules.getAtom(i)).getChildList().getAtom(0));
        }
        P2Harmonic p2SurfaceBond = new P2Harmonic(space, 10000, bondL+1);
        potentialMaster.addPotential(p2SurfaceBond, apiTether, null);

        integrator = new IntegratorVelocityVerlet(potentialMaster, random, 0.005, 100, space);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        integrator.setBox(box);
    }
    
    public static void main(String[] args) {
        Sam sim = new Sam();
        SimulationGraphic simGraphic = new SimulationGraphic(sim, Space3D.getInstance());
        simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(15));
        sim.integrator.setActionInterval(simGraphic.getPaintAction(sim.box), 1);
        ((DisplayBoxCanvasG3DSys)simGraphic.getDisplayBox(sim.box).canvas).setDrawBoundary(DisplayCanvas.DRAW_BOUNDARY_NONE);
        simGraphic.makeAndDisplayFrame();
    }
}
