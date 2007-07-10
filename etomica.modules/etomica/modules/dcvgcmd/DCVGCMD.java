package etomica.modules.dcvgcmd;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeGroup;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.chem.elements.ElementSimple;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.DataSourceGroup;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterProfile;
import etomica.data.types.DataDoubleArray;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.CriterionPositionWall;
import etomica.nbr.CriterionType;
import etomica.nbr.PotentialMasterHybrid;
import etomica.nbr.list.NeighborListManager;
import etomica.box.Box;
import etomica.potential.P2WCA;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;

/**
 * 
 * Dual-control-volume grand-canonical molecular dynamics simulation.
 *
 */
public class DCVGCMD extends Simulation {

    private static final long serialVersionUID = 1L;
    public IntegratorDCVGCMD integratorDCV;
    public P2WCA potential;
    public P2WCA potential1;
    public P1WCAWall potentialwall;
    public P1WCAWall potentialwall1;
    public P1WCAPorousWall potentialwallPorousA, potentialwallPorousA1;
    public P1WCAPorousWall potentialwallPorousB, potentialwallPorousB1;
    public SpeciesSpheresMono species;
    public SpeciesSpheresMono species1;
    public SpeciesTube speciesTube;
    public Box box;
    public DataSourceGroup fluxMeters;
    public MeterFlux meterFlux0, meterFlux1, meterFlux2, meterFlux3;
    public MeterTemperature thermometer;
    public MeterNMolecules density1;
    public MeterNMolecules density2;
    public MeterProfile profile1;
    public MeterProfile profile2;
    public AccumulatorAverage accumulator1;
    public AccumulatorAverage accumulator2;
    public DataPump profile1pump, profile2pump;
    public AccumulatorAverage fluxAccumulator;
    public IVector poreCenter;
    public ActivityIntegrate activityIntegrate;
    public ConfigurationLatticeTube config;
    
    //Constructor
    public DCVGCMD() {
        this(Space3D.getInstance());
    }

    private DCVGCMD(Space space) {
        //Instantiate classes
        super(space, true, new int[] {9,12,11});
        PotentialMasterHybrid potentialMaster = new PotentialMasterHybrid(this, 5.2);
        double mass = 40.;
        double sigma = 3.0;
        double epsilon = 119.8;
        //Default.makeLJDefaults();
        //Default.BOX_SIZE = 14.0;

        species = new SpeciesSpheresMono(this);
        species1 = new SpeciesSpheresMono(this);
        speciesTube = new SpeciesTube(this, 20, 40);
        getSpeciesManager().addSpecies(species);
        getSpeciesManager().addSpecies(species1);
        getSpeciesManager().addSpecies(speciesTube);
        AtomType tubetype = ((AtomTypeGroup) speciesTube.getMoleculeType()).getChildTypes()[0];
        AtomType speciestype = species.getMoleculeType();
        AtomType speciestype1 = species1.getMoleculeType();
        ((ElementSimple)((AtomTypeLeaf)speciestype).getElement()).setMass(mass);
        ((ElementSimple)((AtomTypeLeaf)speciestype1).getElement()).setMass(mass);
        ((AtomTypeSphere)speciestype).setDiameter(sigma);
        ((AtomTypeSphere)speciestype1).setDiameter(sigma);
        ((AtomTypeSphere)tubetype).setDiameter(sigma);
        
        double neighborRangeFac = 1.4;
        potentialMaster.setCellRange(1);
        potentialMaster.setRange(neighborRangeFac * sigma);

        //0-0 intraspecies interaction
        potential = new P2WCA(space, sigma, epsilon);
        potentialMaster.addPotential(potential, new Species[] {species, species});
        
        //1-1 intraspecies interaction
        P2WCA potential11 = new P2WCA(space, sigma, epsilon);
        potentialMaster.addPotential(potential11, new Species[] {species1, species1});

        //0-1 interspecies interaction
        potential1 = new P2WCA(space, sigma, epsilon);
        potentialMaster.addPotential(potential1, new Species[] { species1, species });

        P2WCA potentialTubeAtom = new P2WCA(space, sigma, epsilon);
        potentialMaster.addPotential(potentialTubeAtom,new AtomType[] { tubetype, speciestype});
        
        P2WCA potentialTubeAtom1 = new P2WCA(space, sigma, epsilon);
        potentialMaster.addPotential(potentialTubeAtom1,new AtomType[] { tubetype, speciestype1});

        double neighborRangeFacHalf = (1.0+neighborRangeFac)*0.5;
        
        potentialwall = new P1WCAWall(space, sigma, epsilon);
        CriterionPositionWall criterionWall = new CriterionPositionWall(this);
        criterionWall.setInteractionRange(potentialwall.getRange());
        criterionWall.setNeighborRange(neighborRangeFacHalf*potentialwall.getRange());
        criterionWall.setWallDim(2);
        potentialMaster.addPotential(potentialwall, new Species[] { species });
        potentialMaster.getPotentialMasterList().setCriterion(potentialwall, new CriterionType(criterionWall, speciestype));
        
        potentialwall1 = new P1WCAWall(space, sigma, epsilon);
        CriterionPositionWall criterionWall1 = new CriterionPositionWall(this);
        criterionWall1.setInteractionRange(potentialwall1.getRange());
        criterionWall1.setNeighborRange(neighborRangeFacHalf*potentialwall1.getRange());
        criterionWall1.setWallDim(2);
        potentialMaster.addPotential(potentialwall1, new Species[] { species1 });
        potentialMaster.getPotentialMasterList().setCriterion(potentialwall1, new CriterionType(criterionWall1, speciestype1));

        potentialwallPorousA = new P1WCAPorousWall(space, sigma, epsilon);
        CriterionPositionWall criterionWallA = new CriterionPositionWall(this);
        criterionWallA.setInteractionRange(potentialwallPorousA.getRange());
        criterionWallA.setNeighborRange(neighborRangeFacHalf*potentialwallPorousA.getRange());
        criterionWallA.setWallDim(2);
        potentialMaster.addPotential(potentialwallPorousA, new Species[] { species });
        potentialMaster.getPotentialMasterList().setCriterion(potentialwallPorousA, new CriterionType(criterionWallA, speciestype));
        
        potentialwallPorousA1 = new P1WCAPorousWall(space, sigma, epsilon);
        CriterionPositionWall criterionWallA1 = new CriterionPositionWall(this);
        criterionWallA1.setInteractionRange(potentialwallPorousA1.getRange());
        criterionWallA1.setNeighborRange(neighborRangeFacHalf*potentialwallPorousA1.getRange());
        criterionWallA1.setWallDim(2);
        potentialMaster.addPotential(potentialwallPorousA1, new Species[] { species1 });
        potentialMaster.getPotentialMasterList().setCriterion(potentialwallPorousA1, new CriterionType(criterionWallA1, speciestype1));
        
        potentialwallPorousB = new P1WCAPorousWall(space, sigma, epsilon);
        CriterionPositionWall criterionWallB = new CriterionPositionWall(this);
        criterionWallB.setInteractionRange(potentialwallPorousB.getRange());
        criterionWallB.setNeighborRange(neighborRangeFacHalf*potentialwallPorousB.getRange());
        criterionWallB.setWallDim(2);
        potentialMaster.addPotential(potentialwallPorousB, new Species[] { species });
        potentialMaster.getPotentialMasterList().setCriterion(potentialwallPorousB, new CriterionType(criterionWallB, speciestype));
        
        potentialwallPorousB1 = new P1WCAPorousWall(space, sigma, epsilon);
        CriterionPositionWall criterionWallB1 = new CriterionPositionWall(this);
        criterionWallB1.setInteractionRange(potentialwallPorousB.getRange());
        criterionWallB1.setNeighborRange(neighborRangeFacHalf*potentialwallPorousB.getRange());
        criterionWallB1.setWallDim(2);
        potentialMaster.addPotential(potentialwallPorousB1, new Species[] { species1 });
        potentialMaster.getPotentialMasterList().setCriterion(potentialwallPorousB1, new CriterionType(criterionWallB1, speciestype1));


        box = new Box(this);
        addBox(box);
        box.setNMolecules(species, 20);
        box.setNMolecules(species1, 20);
        box.setNMolecules(speciesTube, 1);
        
        double temperature = Kelvin.UNIT.toSim(500.);
        integratorDCV = new IntegratorDCVGCMD(potentialMaster, temperature, species,
                species1);
        final IntegratorVelocityVerlet integratorMD = new IntegratorVelocityVerlet(this, potentialMaster);
        final IntegratorMC integratorMC = new IntegratorMC(this, potentialMaster);
        integratorDCV.setBox(box);

        /***/
        potentialMaster.setRange(potential.getRange() * neighborRangeFac);
        integratorMC.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());


        activityIntegrate = new ActivityIntegrate(integratorDCV);
        getController().addAction(activityIntegrate);

        integratorDCV.setIntegrators(integratorMC, integratorMD, getRandom());
        integratorMD.setIsothermal(false);
        integratorMD.setMeterTemperature(new MeterTemperature(speciesTube));
        //integrator.setSleepPeriod(1);
        integratorMD.setTimeStep(0.007);
        //integrator.setInterval(10);
        final NeighborListManager nbrManager = potentialMaster.getNeighborManager(box);
        integratorMD.addIntervalAction(nbrManager);
        integratorMD.addNonintervalListener(nbrManager);
        activityIntegrate.setDoSleep(true);
        box.setBoundary(new BoundaryRectangularSlit(this, 2));
//        box.setBoundary(new BoundaryRectangularPeriodic(space));
        box.setDimensions(new Vector3D(40, 40, 80));
        // Crystal crystal = new Crystal(new PrimitiveTetragonal(space, 20,
        // 40),new BasisMonatomic(3));
        double length = 0.25;
        config = new ConfigurationLatticeTube(new LatticeCubicFcc(), length);
        config.initializeCoordinates(box);

        //position of hole in porous-wall potential
        poreCenter = space.makeVector();
//        poreCenter.Ea1Tv1(0.5, box.getBoundary().getDimensions());
        IVector[] poreCentersVector = new IVector[] { poreCenter };
        potentialwallPorousA.setPoreCenters(poreCentersVector);
        potentialwallPorousA1.setPoreCenters(poreCentersVector);
        potentialwallPorousB.setPoreCenters(poreCentersVector);
        potentialwallPorousB1.setPoreCenters(poreCentersVector);

        //radius of hole in porous-wall potential
        double poreRadius = 1.05 * ((ConformationTube)((AtomTypeGroup)speciesTube.getMoleculeType())
                .getConformation()).tubeRadius;
        potentialwallPorousA.setPoreRadius(poreRadius);
        potentialwallPorousA1.setPoreRadius(poreRadius);
        potentialwallPorousB.setPoreRadius(poreRadius);
        potentialwallPorousB1.setPoreRadius(poreRadius);

        //place porous-wall potentials; put just past the edges of the tube
        double zA = (-0.5 + length + 0.05) * box.getBoundary().getDimensions().x(2);
        double zB = ( 0.5 - length - 0.05) * box.getBoundary().getDimensions().x(2);
        potentialwallPorousA.setZ(zA);
        potentialwallPorousA1.setZ(zA);
        potentialwallPorousB.setZ(zB);
        potentialwallPorousB1.setZ(zB);
        criterionWallA.setWallPosition(zA);
        criterionWallA1.setWallPosition(zA);
        criterionWallB.setWallPosition(zB);
        criterionWallB1.setWallPosition(zB);
        criterionWallA.setBoundaryWall(false);
        criterionWallA1.setBoundaryWall(false);
        criterionWallB.setBoundaryWall(false);
        criterionWallB1.setBoundaryWall(false);

        MyMCMove[] moves = integratorDCV.mcMoves();
        meterFlux0 = new MeterFlux(moves[0], integratorDCV);
        meterFlux1 = new MeterFlux(moves[1], integratorDCV);
        meterFlux2 = new MeterFlux(moves[2], integratorDCV);
        meterFlux3 = new MeterFlux(moves[3], integratorDCV);
        meterFlux0.setBox(box);
        meterFlux1.setBox(box);
        meterFlux2.setBox(box);
        meterFlux3.setBox(box);
        fluxMeters = new DataSourceGroup(new MeterFlux[] { meterFlux0, meterFlux1, meterFlux2, meterFlux3 });

        thermometer = new MeterTemperature(speciesTube);
        thermometer.setBox(box);

        density1 = new MeterNMolecules();
        density2 = new MeterNMolecules();
        density1.setBox(box);
        density2.setBox(box);
        density1.setSpecies(species);
        density2.setSpecies(species1);
        profile1 = new MeterProfile(space);
        profile2 = new MeterProfile(space);
        profile1.setDataSource(density1);
        profile1.setBox(box);
        profile1.setProfileVector(new Vector3D(0.0, 0.0, 1.0));
        profile2.setDataSource(density2);
        profile2.setBox(box);
        profile2.setProfileVector(new Vector3D(0.0, 0.0, 1.0));

        accumulator1 = new AccumulatorAverage();
        profile1pump = new DataPump(profile1, accumulator1);
        integratorDCV.addIntervalAction(profile1pump);

        accumulator2 = new AccumulatorAverage();
        profile2pump = new DataPump(profile2, accumulator2);
        integratorDCV.addIntervalAction(profile2pump);

        potentialMaster.getNbrCellManager(box).assignCellAll();
    }
}
