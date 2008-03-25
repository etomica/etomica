package etomica.modules.sam;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomSet;
import etomica.api.IAtomType;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IPotentialMaster;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.atom.AtomTypeSphere;
import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.atom.iterator.Atomset4IteratorIndexList;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.Conformation;
import etomica.config.ConformationChainZigZag;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayCanvas;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.paracetamol.ApiIndexList;
import etomica.potential.P2Harmonic;
import etomica.potential.P2LennardJones;
import etomica.potential.P3BondAngle;
import etomica.potential.P4BondTorsion;
import etomica.potential.PotentialGroup;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.virial.SpeciesAlkane;

/**
 * Self-assembled monolayer module.
 * @author Andrew Schultz
 *
 */
public class Sam extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public SpeciesAlkane species;
    public IBox box;
    public IntegratorVelocityVerlet integrator;
    public ActivityIntegrate activityIntegrate;
    public P1WCAWall wallPotential;
    public ConfigurationSAM config;
    public P2LennardJones p2CH2, p2CH3, p2CH2CH3;
    public P2Harmonic p2Bond;
    public P3BondAngle p3Bond;
    public P4BondTorsion p4Bond;
    public PotentialGroup p1Intra;
    public int chainLength;
    
    public Sam() {
        super(Space.getInstance(3));
        PotentialMaster potentialMaster = new PotentialMaster(space); //List(this, 2.0);

        int nCellX = 10;
        int nCellZ = 10;
        double sizeCellX = 4;
        double sizeCellZ = 4;
        chainLength = 7;

        double surfaceSigma = 4.0;
        double sigma = 4.0;
        double bondL = 3.0;
        double epsilon = Kelvin.UNIT.toSim(100);

        //controller and integrator

	    //species and potentials
//	    species = new SpeciesSpheres(this, chainLength, new ElementSimple(this), new ConformationChain3D(space, new IVector[]{space.makeVector(new double[]{0, bondL, 0})}), space);
//	    ((AtomTypeSphere)species.getLeafType()).setDiameter(sigma);
        species = new SpeciesAlkane(this, space, chainLength);
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

        bondL = 1.54;
        double bondTheta = Math.PI*114/180;
        IVector vector1 = space.makeVector();
        vector1.setX(0, Math.cos(bondTheta/2)*bondL);
        vector1.setX(1, Math.sin(bondTheta/2)*bondL);
        IVector vector2 = space.makeVector();
        vector2.setX(0, -Math.cos(bondTheta/2)*bondL);
        vector2.setX(1, Math.sin(bondTheta/2)*bondL);
        Conformation conformation = new ConformationChainZigZag(space, vector1, vector2);
        species.setConformation(conformation);

        config = new ConfigurationSAM(this, space, species, speciesSurface);
        config.setBasisMolecules(new BasisMonatomic(space));
        config.setBasisSurface(new BasisMonatomic(space));
        config.setCellSizeX(4);
        config.setCellSizeZ(4);
        config.setNCellsX(4);
        config.setNCellsZ(4);
        config.setYOffset(3.3);

        config.initializeCoordinates(box);

        integrator = new IntegratorVelocityVerlet(potentialMaster, random, 0.001, 0, space);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        integrator.setBox(box);

//        P2LennardJones p2 = new P2LennardJones(space, sigma, epsilon);
//        potentialMaster.addPotential(p2, new IAtomTypeLeaf[]{species.getLeafType(), secies.getLeafType()});
//        potentialMaster.addPotential(p2, new IAtomTypeLeaf[]{speciesSurface.getLeafType(), species.getLeafType()});
//        P2Harmonic p2Bond = new P2Harmonic(space, 10000, bondL);
//        AtomsetIteratorBasisDependent bondIterator = ApiBuilder.makeAdjacentPairIterator();
//        AtomsetIteratorBasisDependent nonbondedIterator = ApiBuilder.makeNonAdjacentPairIterator();
//        PotentialGroup pIntra = potentialMaster.makePotentialGroup(1);
//        pIntra.addPotential(p2Bond, bondIterator);
//        pIntra.addPotential(p2, nonbondedIterator);
//        potentialMaster.addPotential(pIntra, new ISpecies[]{species});

        IAtomTypeLeaf typeCH2 = species.getCH2Type();
        IAtomTypeLeaf typeCH3 = species.getCH3Type();
        
        double sigmaCH2 = 3.95;
        double sigmaCH3 = 3.75;

        ((AtomTypeSphere)species.getCH2Type()).setDiameter(sigmaCH2);
        ((AtomTypeSphere)species.getCH3Type()).setDiameter(sigmaCH3);

        double epsilonCH2 = Kelvin.UNIT.toSim(47);
        double epsilonCH3 = Kelvin.UNIT.toSim(98.0);
        double epsilonCH2CH3 = Math.sqrt(epsilonCH2*epsilonCH3);
        p2CH2 = new P2LennardJones(space, sigmaCH2, epsilonCH2);
        p2CH3 = new P2LennardJones(space, sigmaCH3, epsilonCH3);
        p2CH2CH3 = new P2LennardJones(space, 0.5*(sigmaCH2+sigmaCH3), epsilonCH2CH3);
        
        PotentialGroup pGroup = potentialMaster.makePotentialGroup(2);
        pGroup.addPotential(p2CH2, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeCH2, typeCH2}));
        pGroup.addPotential(p2CH2CH3, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeCH2, typeCH3}));
        pGroup.addPotential(p2CH2CH3, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeCH3, typeCH2}));
        pGroup.addPotential(p2CH3, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeCH3, typeCH3}));
        potentialMaster.addPotential(pGroup, new ISpecies[]{species,species});
        p1Intra = potentialMaster.makePotentialGroup(1);
        potentialMaster.addPotential(p1Intra, new ISpecies[]{species});

        p2Bond = new P2Harmonic(space, 10000, bondL);
        p3Bond = new P3BondAngle(space);
        p3Bond.setAngle(Math.PI*114.0/180.0);
        p3Bond.setEpsilon(Kelvin.UNIT.toSim(62500));
        p4Bond = new P4BondTorsion(space, Kelvin.UNIT.toSim(355.03), Kelvin.UNIT.toSim(-68.19), Kelvin.UNIT.toSim(791.32));
        setChainLength(chainLength);

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

        wallPotential = new P1WCAWall(space, 1, 10, 1000);
        wallPotential.setWallPosition(20);
        potentialMaster.addPotential(wallPotential, new IAtomTypeLeaf[]{species.getCH2Type()});
        potentialMaster.addPotential(wallPotential, new IAtomTypeLeaf[]{species.getCH3Type()});
        
        P2LennardJones p2Surface = new P2LennardJones(space, 4.0, Kelvin.UNIT.toSim(50));
        potentialMaster.addPotential(p2Surface, new IAtomTypeLeaf[]{speciesSurface.getLeafType(), species.getCH2Type()});

    }
    
    public void setChainLength(int newChainLength) {
        if (newChainLength < 7) {
            throw new RuntimeException("too short!");
        }
        chainLength = newChainLength;
        IPotentialMaster potentialMaster = integrator.getPotential();
        potentialMaster.removePotential(p1Intra);
        p1Intra = integrator.getPotential().makePotentialGroup(1);
        potentialMaster.addPotential(p1Intra, new ISpecies[]{species});

        AtomsetIteratorBasisDependent bondIterator = ApiBuilder.makeAdjacentPairIterator();
        p1Intra.addPotential(p2Bond, bondIterator);

        int[][] triplets = new int[chainLength-2][3];
        for (int i=0; i<chainLength-2; i++) {
            triplets[i][0] = i;
            triplets[i][1] = i+1;
            triplets[i][2] = i+2;
        }
        p1Intra.addPotential(p3Bond, new Atomset3IteratorIndexList(triplets));

        int[][] quads = new int[chainLength-3][4];
        for (int i=0; i<chainLength-3; i++) {
            quads[i][0] = i;
            quads[i][1] = i+1;
            quads[i][2] = i+2;
            quads[i][3] = i+3;
        }
        p1Intra.addPotential(p4Bond, new Atomset4IteratorIndexList(quads));

        
        p1Intra.addPotential(p2CH3,new ApiIndexList(new int[][]{{0,chainLength-1}}));

        int[][] pairs = new int[2*(chainLength-5)][2];
        for (int i=0; i<chainLength-5; i++) {
            pairs[2*i][0] = 0;
            pairs[2*i][1] = chainLength-2-i;
            pairs[2*i+1][0] = chainLength-1;
            pairs[2*i+1][1] = i+1;
        }
        p1Intra.addPotential(p2CH3,new ApiIndexList(pairs));

        pairs = new int[(chainLength-6)*(chainLength-5)/2][2];
        int k = 0;
        for (int i=1; i<chainLength-5; i++) {
            for (int j=i+4; j<chainLength-1; j++) {
                pairs[k][0] = i;
                pairs[k][1] = j;
                k++;
            }
        }
        p1Intra.addPotential(p2CH2,new ApiIndexList(pairs));
        
    }

    public static void main(String[] args) {
        Sam sim = new Sam();
        SimulationGraphic simGraphic = new SimulationGraphic(sim, Space3D.getInstance());
        simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(15));
        sim.integrator.setActionInterval(simGraphic.getPaintAction(sim.box), 10);
        ((DisplayBoxCanvasG3DSys)simGraphic.getDisplayBox(sim.box).canvas).setDrawBoundary(DisplayCanvas.DRAW_BOUNDARY_NONE);
        simGraphic.makeAndDisplayFrame();
    }
}
