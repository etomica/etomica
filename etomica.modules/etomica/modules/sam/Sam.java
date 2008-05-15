package etomica.modules.sam;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IAtomType;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IAtomTypeSphere;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IPotentialMaster;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPositionFirstAtom;
import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.atom.iterator.Atomset4IteratorIndexList;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConformationChainZigZag;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayCanvas;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.crystal.Basis;
import etomica.potential.P2Harmonic;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedSwitched;
import etomica.potential.P3BondAngle;
import etomica.potential.P4BondTorsion;
import etomica.potential.PotentialGroup;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Space;
import etomica.space3d.IOrientation3D;
import etomica.space3d.Space3D;
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
    public SpeciesAlkaneThiol species;
    public SpeciesSpheresMono speciesSurface;
    public IBox box;
    public IntegratorVelocityVerlet integrator;
    public ActivityIntegrate activityIntegrate;
    public P1WCAWall wallPotential;
    public ConfigurationSAM config;
    public P2LennardJones p2CH2, p2CH3, p2CH2CH3, p2S, p2SCH2;
    public P2SoftSphericalTruncatedSwitched p2CH2t, p2CH3t, p2CH2CH3t, p2St, p2SCH2t;
    public P2Harmonic p2BondCC, p2BondCS;
    public P3BondAngle p3Bond;
    public P4BondTorsion p4BondCCCC, p4BondCCCS;
    public PotentialGroup p1Intra;
    public ApiTether3 apiTether;
    public int chainLength;
    public double chainTheta, chainPsi, chainPhi;
    public double secondaryChainTheta, secondaryChainPsi, secondaryChainPhi;
    public double bondL_CC;
    public double bondTheta;

    
    public Sam() {
        super(Space.getInstance(3));
        PotentialMaster potentialMaster = new PotentialMaster(space); //List(this, 2.0);

        int nCellX = 4;
        int nCellZ = 2;
        // gold has FCC unit cell, a=4.0782A
        double sizeCellZ = 4.0782/Math.sqrt(2)*3; //Math.sqrt(3)*sizeCellX;
        double sizeCellX = sizeCellZ/Math.sqrt(3);
        chainLength = 14;

        double surfaceSigma = 3.0;

        //controller and integrator

	    //species and potentials
//	    species = new SpeciesSpheres(this, chainLength, new ElementSimple(this), new ConformationChain3D(space, new IVector[]{space.makeVector(new double[]{0, bondL, 0})}), space);
//	    ((AtomTypeSphere)species.getLeafType()).setDiameter(sigma);
        species = new SpeciesAlkaneThiol(this, space, chainLength-1);
        getSpeciesManager().addSpecies(species);

        //construct box
	    box = new Box(new BoundaryRectangularSlit(this, 1, space), space);
        addBox(box);
        IVector dim = space.makeVector();
        dim.E(new double[]{sizeCellX*nCellX, chainLength*2, sizeCellZ*nCellZ});
        box.setDimensions(dim);
        box.setNMolecules(species, nCellX*nCellZ);

        speciesSurface = new SpeciesSpheresMono(this, space);
        ((IAtomTypeSphere)speciesSurface.getLeafType()).setDiameter(surfaceSigma);
        ((ElementSimple)speciesSurface.getLeafType().getElement()).setMass(Double.POSITIVE_INFINITY);
        speciesSurface.setPositionDefinition(new AtomPositionFirstAtom());
        getSpeciesManager().addSpecies(speciesSurface);

        bondL_CC = 1.54;
        double bondL_CS = 1.82;
        bondTheta = Math.PI*114/180;
        ConformationChainZigZag conformation = new ConformationChainZigZag(space);
        species.setConformation(conformation);
        chainTheta = 0;
        chainPsi = 0;
        chainPhi = 0;

        config = new ConfigurationSAM(this, space, species, speciesSurface);
        Basis alkaneBasis = new Basis(new IVector[]{space.makeVector(new double[]{1.0/6.0,0,1.0/6.0}), ((Space)space).makeVector(new double[]{2.0/3.0, 0, 2.0/3.0})});
        Basis surfaceBasis = new Basis(new IVector[]{
                space.makeVector(new double[]{2.0/6.0, 0, 0}),
                space.makeVector(new double[]{5.0/6.0, 0, 1.0/6.0}),
                space.makeVector(new double[]{2.0/6.0, 0, 2.0/6.0}),
                space.makeVector(new double[]{5.0/6.0, 0, 3.0/6.0}),
                space.makeVector(new double[]{2.0/6.0, 0, 4.0/6.0}),
                space.makeVector(new double[]{5.0/6.0, 0, 5.0/6.0})});
        config.setBasisMolecules(alkaneBasis);
        config.setBasisSurface(surfaceBasis);
        config.setCellSizeX(sizeCellX);
        config.setCellSizeZ(sizeCellZ);
        config.setNCellsX(nCellX);
        config.setNCellsZ(nCellZ);
        config.setSurfaceYOffset(2);

        updateConformation();
        ConformationChainZigZag secondaryConformation = new ConformationChainZigZag(space);
        species.setConformation(secondaryConformation);
        updateConformation();
        species.setConformation(conformation);
        config.setSecondaryConformation(secondaryConformation);

        config.initializeCoordinates(box);

        integrator = new IntegratorVelocityVerlet(potentialMaster, random, 0.005, Kelvin.UNIT.toSim(300), space);
        integrator.setIsothermal(true);
        integrator.setThermostatInterval(500);
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
        IAtomTypeLeaf typeS = species.getSulfurType();
        
        double sigmaCH2 = 3.95;
        double sigmaCH3 = 3.75;
        double sigmaSulfur = 3.62;

        ((IAtomTypeSphere)typeCH2).setDiameter(sigmaCH2);
        ((IAtomTypeSphere)typeCH3).setDiameter(sigmaCH3);
        ((IAtomTypeSphere)typeS).setDiameter(sigmaSulfur);

        double epsilonCH2 = Kelvin.UNIT.toSim(46);
        double epsilonCH3 = Kelvin.UNIT.toSim(98);
        double epsilonSulfur = Kelvin.UNIT.toSim(232);
        double epsilonCH2CH3 = Math.sqrt(epsilonCH2*epsilonCH3);
        double epsilonCH2Sulfur = Math.sqrt(epsilonCH2*epsilonSulfur);
        // sulfur and CH3 will never be close
        //double epsilonCH3Sulfur = Math.sqrt(epsilonCH3*epsilonSulfur);
        double rCut = box.getBoundary().getDimensions().x(0)*0.5;
        if (0.5*box.getBoundary().getDimensions().x(2) < rCut) {
            rCut = 0.5*box.getBoundary().getDimensions().x(2);
        }
        p2CH2 = new P2LennardJones(space, sigmaCH2, epsilonCH2);
        p2CH3 = new P2LennardJones(space, sigmaCH3, epsilonCH3);
        p2S = new P2LennardJones(space, sigmaSulfur, epsilonSulfur);
        p2CH2CH3 = new P2LennardJones(space, 0.5*(sigmaCH2+sigmaCH3), epsilonCH2CH3);
        p2SCH2 = new P2LennardJones(space, 0.5*(sigmaSulfur+sigmaCH2), epsilonCH2Sulfur);
        p2CH2t = new P2SoftSphericalTruncatedSwitched(space, p2CH2, rCut);
        p2CH3t = new P2SoftSphericalTruncatedSwitched(space, p2CH2, rCut);
        p2CH2CH3t = new P2SoftSphericalTruncatedSwitched(space, p2CH2, rCut);
        p2St = new P2SoftSphericalTruncatedSwitched(space, p2S, rCut);
        p2SCH2t = new P2SoftSphericalTruncatedSwitched(space, p2SCH2, rCut);
        
        PotentialGroup pGroup = potentialMaster.makePotentialGroup(2);
        pGroup.addPotential(p2CH2t, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeCH2, typeCH2}));
        pGroup.addPotential(p2CH3t, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeCH3, typeCH3}));
        pGroup.addPotential(p2St, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeS, typeS}));
        pGroup.addPotential(p2SCH2t, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeS, typeCH2}));
        pGroup.addPotential(p2CH2CH3t, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeCH2, typeCH3}));
        pGroup.addPotential(p2CH2CH3t, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeCH3, typeCH2}));
        potentialMaster.addPotential(pGroup, new ISpecies[]{species,species});
        p1Intra = potentialMaster.makePotentialGroup(1);
        potentialMaster.addPotential(p1Intra, new ISpecies[]{species});

        p2BondCC = new P2Harmonic(space, 10000, bondL_CC);
        p2BondCS = new P2Harmonic(space, 10000, bondL_CS);
        // bond angle potential is the same for CCC and CCS
        p3Bond = new P3BondAngle(space);
        p3Bond.setAngle(Math.PI*114.0/180.0);
        p3Bond.setEpsilon(Kelvin.UNIT.toSim(62500));
        p4BondCCCC = new P4BondTorsion(space, 0, Kelvin.UNIT.toSim(355.03), Kelvin.UNIT.toSim(-68.19), Kelvin.UNIT.toSim(791.32));
        p4BondCCCS = new P4BondTorsion(space, Kelvin.UNIT.toSim(-251.06), Kelvin.UNIT.toSim(428.73), Kelvin.UNIT.toSim(-111.85), Kelvin.UNIT.toSim(441.27));
        setChainLength(chainLength);

        apiTether = new ApiTether3(this, species);
        apiTether.setBox(box);
        findTetherBonds();
        P2Harmonic p2SurfaceBond = new P2Harmonic(space, 10000, 2);
        potentialMaster.addPotential(p2SurfaceBond, apiTether, null);

        wallPotential = new P1WCAWall(space, 1, 4, 1000);
        wallPotential.setWallPosition(box.getBoundary().getDimensions().x(1));
        potentialMaster.addPotential(wallPotential, new IAtomTypeLeaf[]{species.getCH2Type()});
        potentialMaster.addPotential(wallPotential, new IAtomTypeLeaf[]{species.getCH3Type()});
        
        P2LennardJones p2Surface = new P2LennardJones(space, 3.0, Kelvin.UNIT.toSim(50));
        potentialMaster.addPotential(p2Surface, new IAtomTypeLeaf[]{speciesSurface.getLeafType(), species.getCH2Type()});
        potentialMaster.addPotential(p2Surface, new IAtomTypeLeaf[]{speciesSurface.getLeafType(), species.getSulfurType()});
        potentialMaster.addPotential(p2Surface, new IAtomTypeLeaf[]{speciesSurface.getLeafType(), species.getCH3Type()});
    }

    protected void updateConformation() {
        double bondTheta0 = chainTheta + .5*(Math.PI - bondTheta);
        IVector vector1 = ((ConformationChainZigZag)species.getConformation()).getFirstVector();
        vector1.setX(0, Math.cos(chainPsi)*Math.sin(bondTheta0)*bondL_CC);
        vector1.setX(1, Math.cos(bondTheta0)*bondL_CC);
        vector1.setX(2, Math.sin(chainPsi)*Math.sin(bondTheta0)*bondL_CC);
        double bondTheta2 = bondTheta0 - (Math.PI - bondTheta);
        IVector vector2 = ((ConformationChainZigZag)species.getConformation()).getSecondVector();
        vector2.setX(0, Math.cos(chainPsi)*(Math.sin(bondTheta2))*bondL_CC);
        vector2.setX(1, Math.cos(bondTheta2)*bondL_CC);
        vector2.setX(2, Math.sin(chainPsi)*(Math.sin(bondTheta2))*bondL_CC);
        
        IVector vector0 = space.makeVector();
        vector0.Ev1Pv2(vector1, vector2);
        IOrientation3D orientation = (IOrientation3D)space.makeOrientation();
        orientation.setDirection(vector1);
        IVector vector0Axis = space.makeVector();
        vector0Axis.Ea1Tv1(1.0/Math.sqrt(vector0.squared()), vector0);
        orientation.rotateBy(chainPhi, vector0Axis);
        vector1.Ea1Tv1(Math.sqrt(vector1.squared()), orientation.getDirection());
        vector2.Ev1Mv2(vector0, vector1);
        
        IMolecule molecule = species.makeMolecule();
        IVector moleculePos = space.makeVector();
        moleculePos.E(molecule.getType().getPositionDefinition().position(molecule));
        IVector sulfurPosition = ((IAtomPositioned)molecule.getChildList().getAtom(0)).getPosition();
        sulfurPosition.ME(moleculePos);
        molecule = null;
//        sulfurPosition.setX(1, 0);
        sulfurPosition.TE(-1);
        sulfurPosition.setX(1, sulfurPosition.x(1)+2.5);
        
        config.setMoleculeOffset(sulfurPosition);
    }
    
    public void setChainTheta(double newTheta) {
        chainTheta = newTheta;
        updateConformation();
    }
    
    public double getChainTheta() {
        return chainTheta;
    }
    
    public void setChainPsi(double newPsi) {
        chainPsi = newPsi;
        updateConformation();
    }
    
    public double getChainPsi() {
        return chainPsi;
    }
    
    public void setChainPhi(double newPhi) {
        chainPhi = newPhi;
        updateConformation();
    }
    
    public double getChainPhi() {
        return chainPhi;
    }
    
    protected void updateSecondaryConformation() {
        double bondTheta0 = secondaryChainTheta + .5*(Math.PI - bondTheta);
        IVector vector1 = ((ConformationChainZigZag)config.getSecondaryConformation()).getFirstVector();
        vector1.setX(0, Math.cos(secondaryChainPsi)*Math.sin(bondTheta0)*bondL_CC);
        vector1.setX(1, Math.cos(bondTheta0)*bondL_CC);
        vector1.setX(2, Math.sin(secondaryChainPsi)*Math.sin(bondTheta0)*bondL_CC);
        double bondTheta2 = bondTheta0 - (Math.PI - bondTheta);
        IVector vector2 = ((ConformationChainZigZag)config.getSecondaryConformation()).getSecondVector();
        vector2.setX(0, Math.cos(secondaryChainPsi)*(Math.sin(bondTheta2))*bondL_CC);
        vector2.setX(1, Math.cos(bondTheta2)*bondL_CC);
        vector2.setX(2, Math.sin(secondaryChainPsi)*(Math.sin(bondTheta2))*bondL_CC);

        IVector vector0 = space.makeVector();
        vector0.Ev1Pv2(vector1, vector2);
        IOrientation3D orientation = (IOrientation3D)space.makeOrientation();
        orientation.setDirection(vector1);
        IVector vector0Axis = space.makeVector();
        vector0Axis.Ea1Tv1(1.0/Math.sqrt(vector0.squared()), vector0);
        orientation.rotateBy(secondaryChainPhi, vector0Axis);
        vector1.Ea1Tv1(Math.sqrt(vector1.squared()), orientation.getDirection());
        vector2.Ev1Mv2(vector0, vector1);
    }
        
    public void setSecondaryChainTheta(double newSecondaryTheta) {
        secondaryChainTheta = newSecondaryTheta;
        updateSecondaryConformation();
    }
    
    public double getSecondaryChainTheta() {
        return secondaryChainTheta;
    }
    
    public void setSecondaryChainPsi(double newSecondaryPsi) {
        secondaryChainPsi = newSecondaryPsi;
        updateSecondaryConformation();
    }
    
    public double getSecondaryChainPsi() {
        return secondaryChainPsi;
    }
    
    public void setSecondaryChainPhi(double newSecondaryPhi) {
        secondaryChainPhi = newSecondaryPhi;
        updateSecondaryConformation();
    }
    
    public double getSecondaryChainPhi() {
        return secondaryChainPhi;
    }
    
    public void findTetherBonds() {
        IAtomSet polymerMolecules = box.getMoleculeList(species);
        IAtomSet surfaceMolecules = box.getMoleculeList(speciesSurface);
        int nMolecules = box.getNMolecules(species);
        double maxDistance = 3.5*3.5;
        IVector dr = space.makeVector();
        IBoundary boundary = box.getBoundary();
        for (int i=0; i<nMolecules; i++) {
            AtomArrayList bondedSurfaceAtoms = new AtomArrayList(3);
            IAtomPositioned sulfur = (IAtomPositioned)((IMolecule)polymerMolecules.getAtom(i)).getChildList().getAtom(0);
            for (int j=0; j<surfaceMolecules.getAtomCount(); j++) {
                IAtomPositioned gold = (IAtomPositioned)((IMolecule)surfaceMolecules.getAtom(j)).getChildList().getAtom(0);
                dr.Ev1Mv2(sulfur.getPosition(), gold.getPosition());
                boundary.nearestImage(dr);
                if (dr.squared() < maxDistance) {
                    bondedSurfaceAtoms.add(gold);
                }
            }
            if (bondedSurfaceAtoms.getAtomCount() != 3) {
                throw new RuntimeException("only found "+bondedSurfaceAtoms.getAtomCount()+" bonded atoms");
            }
            apiTether.setBondedSurfaceAtoms((IMolecule)polymerMolecules.getAtom(i), bondedSurfaceAtoms);
        }
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

        int[][] pairs = new int[1][2];
        pairs[0][0] = 0;
        pairs[0][1] = 1;
        p1Intra.addPotential(p2BondCS, new ApiIndexList(pairs));
        pairs = new int[chainLength-2][2];
        for (int i=1; i<chainLength-1; i++) {
            pairs[i-1][0] = i;
            pairs[i-1][1] = i+1;
        }
        p1Intra.addPotential(p2BondCC, new ApiIndexList(pairs));

        // CCC and CCS bond is the same
        int[][] triplets = new int[chainLength-2][3];
        for (int i=0; i<chainLength-2; i++) {
            triplets[i][0] = i;
            triplets[i][1] = i+1;
            triplets[i][2] = i+2;
        }
        p1Intra.addPotential(p3Bond, new Atomset3IteratorIndexList(triplets));

        int[][] quads = new int[][]{{0,1,2,3}};
        p1Intra.addPotential(p4BondCCCS, new Atomset4IteratorIndexList(quads));
        quads = new int[chainLength-3][4];
        for (int i=0; i<chainLength-3; i++) {
            quads[i][0] = i;
            quads[i][1] = i+1;
            quads[i][2] = i+2;
            quads[i][3] = i+3;
        }
        p1Intra.addPotential(p4BondCCCC, new Atomset4IteratorIndexList(quads));

        pairs = new int[chainLength-5][2];
        for (int i=1; i<chainLength-4; i++) {
            pairs[i-1][0] = chainLength-1;
            pairs[i-1][1] = i;
        }
        p1Intra.addPotential(p2CH2CH3t,new ApiIndexList(pairs));

        pairs = new int[chainLength-5][2];
        for (int i=4; i<chainLength-1; i++) {
            pairs[i-4][0] = 0;
            pairs[i-4][1] = i;
        }
        p1Intra.addPotential(p2SCH2t,new ApiIndexList(pairs));

        pairs = new int[(chainLength-6)*(chainLength-5)/2][2];
        int k = 0;
        for (int i=1; i<chainLength-5; i++) {
            for (int j=i+4; j<chainLength-1; j++) {
                pairs[k][0] = i;
                pairs[k][1] = j;
                k++;
            }
        }
        p1Intra.addPotential(p2CH2t,new ApiIndexList(pairs));
    }

    public static void main(String[] args) {
        Sam sim = new Sam();
        SimulationGraphic simGraphic = new SimulationGraphic(sim, Space3D.getInstance());
        simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(15));
//        sim.integrator.setActionInterval(simGraphic.getPaintAction(sim.box), 10);
        ((DisplayBoxCanvasG3DSys)simGraphic.getDisplayBox(sim.box).canvas).setDrawBoundary(DisplayCanvas.DRAW_BOUNDARY_NONE);
        simGraphic.makeAndDisplayFrame();
    }
}
