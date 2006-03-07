package etomica.modules.dcvgcmd;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomType;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.DataSourceGroup;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterProfile;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntervalActionAdapter;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.CriterionPositionWall;
import etomica.nbr.CriterionSimple;
import etomica.nbr.CriterionSpecies;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.PotentialMasterHybrid;
import etomica.nbr.list.NeighborListManager;
import etomica.phase.Phase;
import etomica.potential.P2WCA;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;
import etomica.util.Default;

/**
 * 
 * Dual-control-volume grand-canonical molecular dynamics simulation.
 *
 */
public class DCVGCMD extends Simulation {

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
    public Phase phase;
    public DataSourceGroup fluxMeters;
    public MeterFlux meterFlux0, meterFlux1, meterFlux2, meterFlux3;
    public MeterTemperature thermometer;
    public MeterNMolecules density1;
    public MeterNMolecules density2;
    public MeterProfile profile1;
    public MeterProfile profile2;
    public AccumulatorAverage accumulator1;
    public AccumulatorAverage accumulator2;
    public AccumulatorAverage fluxAccumulator;
    public Vector poreCenter;
    public ActivityIntegrate activityIntegrate;
    
    //Constructor
    public DCVGCMD() {
        this(Space3D.getInstance());
    }

    private DCVGCMD(Space space) {
        //Instantiate classes
        super(space, true, new PotentialMasterHybrid(space,5.2),new int[] {1,4,4,12,11,0},new Default());
        defaults.atomMass = 40.;
        defaults.atomSize = 3.0;
        defaults.potentialWell = 119.8;
        //Default.makeLJDefaults();
        //Default.BOX_SIZE = 14.0;

        species = new SpeciesSpheresMono(this);
        species1 = new SpeciesSpheresMono(this);
        speciesTube = new SpeciesTube(this, 20, 40);
        AtomType tubetype = ((AtomFactoryHomo) speciesTube.moleculeFactory())
                .getChildFactory().getType();
        AtomType speciestype = species.moleculeFactory().getType();
        AtomType speciestype1 = species1.moleculeFactory().getType();
        
        PotentialMasterHybrid potentialMasterHybrid = (PotentialMasterHybrid) potentialMaster;
        double neighborRangeFac = 1.4;
        final NeighborListManager nbrManager = potentialMasterHybrid.getNeighborManager();
        nbrManager.getPbcEnforcer().setApplyToMolecules(false);
        potentialMasterHybrid.setCellRange(1);
        potentialMasterHybrid.setRange(neighborRangeFac * defaults.atomSize);

        //0-0 intraspecies interaction
        potential = new P2WCA(this);
        NeighborCriterion nbrCriterion = new CriterionSimple(this,potential.getRange(),neighborRangeFac*potential.getRange());
        CriterionSpecies criterion = new CriterionSpecies(nbrCriterion, species, species);
        potential.setCriterion(criterion);
//        nbrManager.addCriterion(nbrCriterion,new AtomType[]{species.getFactory().getType()});
        potentialMaster.addPotential(potential, new Species[] {species, species});
        
        //1-1 intraspecies interaction
        P2WCA potential11 = new P2WCA(this);
        nbrCriterion = new CriterionSimple(this,potential.getRange(),neighborRangeFac*potential11.getRange());
        criterion = new CriterionSpecies(nbrCriterion, species1, species1);
        potential11.setCriterion(criterion);
//        nbrManager.addCriterion(nbrCriterion,new AtomType[]{species1.getFactory().getType()});
        potentialMaster.addPotential(potential11, new Species[] {species1, species1});

        //0-1 interspecies interaction
        potential1 = new P2WCA(this);
        nbrCriterion = new CriterionSimple(this,potential.getRange(),neighborRangeFac*potential1.getRange());
        criterion = new CriterionSpecies(nbrCriterion, species1, species);
        potential1.setCriterion(criterion);
//        nbrManager.addCriterion(nbrCriterion, new AtomType[]{species.getFactory().getType(),species1.getFactory().getType()});
        potentialMaster.addPotential(potential1, new Species[] { species1, species });

        P2WCA potentialTubeAtom = new P2WCA(this);
        potentialMaster.addPotential(potentialTubeAtom,new AtomType[] { tubetype, speciestype});
        nbrCriterion = new CriterionSimple(this,potentialTubeAtom.getRange(),neighborRangeFac*potentialTubeAtom.getRange());
        criterion = new CriterionSpecies(nbrCriterion, speciesTube, species);
        potentialTubeAtom.setCriterion(criterion);
//        nbrManager.addCriterion(nbrCriterion,new AtomType[]{species.getFactory().getType(),((AtomFactoryHomo)speciesTube.getFactory()).getChildFactory().getType()});
        
        P2WCA potentialTubeAtom1 = new P2WCA(this);
        potentialMaster.addPotential(potentialTubeAtom1,new AtomType[] { tubetype, speciestype1});
        nbrCriterion = new CriterionSimple(this,potentialTubeAtom1.getRange(),neighborRangeFac*potentialTubeAtom.getRange());
        criterion = new CriterionSpecies(nbrCriterion, speciesTube, species1);
        potentialTubeAtom1.setCriterion(criterion);
//        nbrManager.addCriterion(nbrCriterion,new AtomType[]{species1.getFactory().getType()});

        double neighborRangeFacHalf = (1.0+neighborRangeFac)*0.5;
        
        potentialwall = new P1WCAWall(this);
        CriterionPositionWall criterionWall = new CriterionPositionWall(this);
        criterionWall.setInteractionRange(potentialwall.getRange());
        criterionWall.setNeighborRange(neighborRangeFacHalf*potentialwall.getRange());
        criterionWall.setWallDim(2);
        potentialwall.setCriterion(criterionWall);
        potentialMaster.addPotential(potentialwall, new Species[] { species });
//        nbrManager.addCriterion(criterionWall,new AtomType[]{species.getFactory().getType()});
        
        potentialwall1 = new P1WCAWall(this);
        CriterionPositionWall criterionWall1 = new CriterionPositionWall(this);
        criterionWall1.setInteractionRange(potentialwall1.getRange());
        criterionWall1.setNeighborRange(neighborRangeFacHalf*potentialwall1.getRange());
        criterionWall1.setWallDim(2);
        potentialwall1.setCriterion(criterionWall1);
        potentialMaster.addPotential(potentialwall1, new Species[] { species1 });
//        nbrManager.addCriterion(criterionWall1,new AtomType[]{species1.getFactory().getType()});

        potentialwallPorousA = new P1WCAPorousWall(this);
        CriterionPositionWall criterionWallA = new CriterionPositionWall(this);
        criterionWallA.setInteractionRange(potentialwallPorousA.getRange());
        criterionWallA.setNeighborRange(neighborRangeFacHalf*potentialwallPorousA.getRange());
        criterionWallA.setWallDim(2);
        potentialwallPorousA.setCriterion(criterionWallA);
        potentialMaster.addPotential(potentialwallPorousA, new Species[] { species });
//        nbrManager.addCriterion(criterionWallA,new AtomType[]{species.getFactory().getType()});
        
        potentialwallPorousA1 = new P1WCAPorousWall(this);
        CriterionPositionWall criterionWallA1 = new CriterionPositionWall(this);
        criterionWallA1.setInteractionRange(potentialwallPorousA1.getRange());
        criterionWallA1.setNeighborRange(neighborRangeFacHalf*potentialwallPorousA1.getRange());
        criterionWallA1.setWallDim(2);
        potentialwallPorousA1.setCriterion(criterionWallA1);
        potentialMaster.addPotential(potentialwallPorousA1, new Species[] { species1 });
//        nbrManager.addCriterion(criterionWallA1,new AtomType[]{species1.getFactory().getType()});
        
        potentialwallPorousB = new P1WCAPorousWall(this);
        CriterionPositionWall criterionWallB = new CriterionPositionWall(this);
        criterionWallB.setInteractionRange(potentialwallPorousB.getRange());
        criterionWallB.setNeighborRange(neighborRangeFacHalf*potentialwallPorousB.getRange());
        criterionWallB.setWallDim(2);
        potentialwallPorousB.setCriterion(criterionWallB);
        potentialMaster.addPotential(potentialwallPorousB, new Species[] { species });
//        nbrManager.addCriterion(criterionWallB,new AtomType[]{species.getFactory().getType()});
        
        potentialwallPorousB1 = new P1WCAPorousWall(this);
        CriterionPositionWall criterionWallB1 = new CriterionPositionWall(this);
        criterionWallB1.setInteractionRange(potentialwallPorousB.getRange());
        criterionWallB1.setNeighborRange(neighborRangeFacHalf*potentialwallPorousB.getRange());
        criterionWallB1.setWallDim(2);
        potentialwallPorousB1.setCriterion(criterionWallB1);
        potentialMaster.addPotential(potentialwallPorousB1, new Species[] { species1 });
//        nbrManager.addCriterion(criterionWallB1,new AtomType[]{species1.getFactory().getType()});


        species.setNMolecules(20);
        species1.setNMolecules(20);
        phase = new Phase(this);

        double temperature = Kelvin.UNIT.toSim(500.);
        integratorDCV = new IntegratorDCVGCMD(potentialMaster, temperature, species,
                species1);
        final IntegratorVelocityVerlet integratorMD = new IntegratorVelocityVerlet(this);
        final IntegratorMC integratorMC = new IntegratorMC(this);
        integratorDCV.setPhase(phase);

        /***/
        nbrManager.setRange(potential.getRange() * neighborRangeFac);
        integratorMC.getMoveEventManager().addListener(potentialMasterHybrid.getNbrCellManager(phase).makeMCMoveListener());
        potentialMasterHybrid.updateTypeList(phase);


        activityIntegrate = new ActivityIntegrate(this,integratorDCV);
        getController().addAction(activityIntegrate);

        integratorDCV.setIntegrators(integratorMC, integratorMD);
        integratorMD.setIsothermal(false);
        integratorMD.setMeterTemperature(new MeterTemperature(speciesTube));
        //integrator.setSleepPeriod(1);
        integratorMD.setTimeStep(0.007);
        //integrator.setInterval(10);
        integratorMD.addListener(nbrManager);
        activityIntegrate.setDoSleep(true);
        phase.setBoundary(new BoundaryRectangularSlit(this, 2));
//        phase.setBoundary(new BoundaryRectangularPeriodic(space));
        phase.setDimensions(new Vector3D(40, 40, 80));
        // Crystal crystal = new Crystal(new PrimitiveTetragonal(space, 20,
        // 40),new BasisMonatomic(3));
        double length = 0.25;
        ConfigurationLatticeTube config = new ConfigurationLatticeTube(
                new LatticeCubicFcc(), length, speciesTube);
        config.initializeCoordinates(phase);

        //position of hole in porous-wall potential
        poreCenter = space.makeVector();
        poreCenter.Ea1Tv1(0.5, phase.getBoundary().getDimensions());
        Vector[] poreCentersVector = new Vector[] { poreCenter };
        potentialwallPorousA.setPoreCenters(poreCentersVector);
        potentialwallPorousA1.setPoreCenters(poreCentersVector);
        potentialwallPorousB.setPoreCenters(poreCentersVector);
        potentialwallPorousB1.setPoreCenters(poreCentersVector);

        //radius of hole in porous-wall potential
        double poreRadius = 1.05 * ((ConformationTube) speciesTube.getFactory()
                .getConformation()).tubeRadius;
        potentialwallPorousA.setPoreRadius(poreRadius);
        potentialwallPorousA1.setPoreRadius(poreRadius);
        potentialwallPorousB.setPoreRadius(poreRadius);
        potentialwallPorousB1.setPoreRadius(poreRadius);

        //place porous-wall potentials; put just past the edges of the tube
        double zA = (length + 0.05) * phase.getBoundary().getDimensions().x(2);
        double zB = (1.0 - length - 0.05) * phase.getBoundary().getDimensions().x(2);
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
        meterFlux0.setPhase(phase);
        meterFlux1.setPhase(phase);
        meterFlux2.setPhase(phase);
        meterFlux3.setPhase(phase);
        fluxMeters = new DataSourceGroup(new MeterFlux[] { meterFlux0,
                meterFlux1, meterFlux2, meterFlux3 });
//        fluxAccumulator = new AccumulatorAverage();
//        DataPump fluxPump = new DataPump(fluxMeters, fluxAccumulator);
//        IntervalActionAdapter fluxInterval = new IntervalActionAdapter(
//                fluxPump, integratorDCV);

        thermometer = new MeterTemperature(speciesTube);
        thermometer.setPhase(phase);

        density1 = new MeterNMolecules();
        density2 = new MeterNMolecules();
        density1.setPhase(phase);
        density2.setPhase(phase);
        density1.setSpecies(species);
        density2.setSpecies(species1);
        profile1 = new MeterProfile(space);
        profile2 = new MeterProfile(space);
        profile1.setDataSource(density1);
        profile1.setPhase(phase);
        profile1.setProfileVector(new Vector3D(0.0, 0.0, 1.0));
        profile2.setDataSource(density2);
        profile2.setPhase(phase);
        profile2.setProfileVector(new Vector3D(0.0, 0.0, 1.0));

        accumulator1 = new AccumulatorAverage(this);
        DataPump profile1pump = new DataPump(profile1, accumulator1);
        new IntervalActionAdapter(profile1pump, integratorDCV);

        accumulator2 = new AccumulatorAverage(this);
        DataPump profile2pump = new DataPump(profile2, accumulator2);
        new IntervalActionAdapter(profile2pump, integratorDCV);

        ((PotentialMasterHybrid)potentialMaster).getNbrCellManager(phase).assignCellAll();
//remove for nbrlist        integrator.addIntervalListener(new PhaseImposePbc(phase));
    } //End of constructor

} //End of DCVGCMD class
