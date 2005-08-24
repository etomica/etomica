package etomica.modules.dcvgcmd;

import etomica.Default;
import etomica.Phase;
import etomica.Simulation;
import etomica.Species;
import etomica.SpeciesSpheresMono;
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
import etomica.nbr.CriterionSimple;
import etomica.nbr.CriterionSpecies;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.PotentialCalculationAgents;
import etomica.nbr.PotentialMasterHybrid;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.nbr.list.NeighborListManager;
import etomica.potential.P2WCA;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.units.Kelvin;

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
        super(space, true, new PotentialMasterHybrid(space));
        Default.ATOM_MASS = 40.;
        Default.ATOM_SIZE = 3.0;
        Default.POTENTIAL_WELL = 119.8;
        //Default.makeLJDefaults();
        //Default.BOX_SIZE = 14.0;

        species = new SpeciesSpheresMono(this);
        species1 = new SpeciesSpheresMono(this);
        speciesTube = new SpeciesTube(this, 20, 40);
        AtomType tubetype = ((AtomFactoryHomo) speciesTube.moleculeFactory())
                .childFactory().getType();
        AtomType speciestype = species.moleculeFactory().getType();
        AtomType speciestype1 = species1.moleculeFactory().getType();
        
        PotentialMasterHybrid potentialMasterHybrid = (PotentialMasterHybrid) potentialMaster;
        double neighborRangeFac = 1.4;
        final NeighborListManager nbrManager = potentialMasterHybrid.getNeighborManager();
        nbrManager.getPbcEnforcer().setApplyToMolecules(false);
        int nCells = (int) (40 / (neighborRangeFac * Default.ATOM_SIZE));
        potentialMasterHybrid.setNCells(nCells);

        //0-0 intraspecies interaction
        potential = new P2WCA(space);
        NeighborCriterion nbrCriterion = new CriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());
        CriterionSpecies criterion = new CriterionSpecies(nbrCriterion, species, species);
        potential.setCriterion(criterion);
        nbrManager.addCriterion(nbrCriterion,new AtomType[]{species.getFactory().getType()});
        potentialMaster.setSpecies(potential, new Species[] {species, species});
        
        //1-1 intraspecies interaction
        P2WCA potential11 = new P2WCA(space);
        nbrCriterion = new CriterionSimple(space,potential.getRange(),neighborRangeFac*potential11.getRange());
        criterion = new CriterionSpecies(nbrCriterion, species1, species1);
        potential11.setCriterion(criterion);
        nbrManager.addCriterion(nbrCriterion,new AtomType[]{species1.getFactory().getType()});
        potentialMaster.setSpecies(potential11, new Species[] {species1, species1});

        //0-1 interspecies interaction
        potential1 = new P2WCA(space);
        nbrCriterion = new CriterionSimple(space,potential.getRange(),neighborRangeFac*potential1.getRange());
        criterion = new CriterionSpecies(nbrCriterion, species1, species);
        potential1.setCriterion(criterion);
        nbrManager.addCriterion(nbrCriterion, new AtomType[]{species.getFactory().getType(),species1.getFactory().getType()});
        potentialMaster.setSpecies(potential1, new Species[] { species1, species });

        P2WCA potentialTubeAtom = new P2WCA(space);
        potentialMaster.addPotential(potentialTubeAtom,new AtomType[] { tubetype, speciestype});
        nbrCriterion = new CriterionSimple(space,potentialTubeAtom.getRange(),neighborRangeFac*potentialTubeAtom.getRange());
        criterion = new CriterionSpecies(nbrCriterion, speciesTube, species);
        potentialTubeAtom.setCriterion(criterion);
        nbrManager.addCriterion(nbrCriterion,new AtomType[]{species.getFactory().getType(),((AtomFactoryHomo)speciesTube.getFactory()).childFactory().getType()});
        
        P2WCA potentialTubeAtom1 = new P2WCA(space);
        potentialMaster.addPotential(potentialTubeAtom1,new AtomType[] { tubetype, speciestype1});
        nbrCriterion = new CriterionSimple(space,potentialTubeAtom1.getRange(),neighborRangeFac*potentialTubeAtom.getRange());
        criterion = new CriterionSpecies(nbrCriterion, speciesTube, species1);
        potentialTubeAtom1.setCriterion(criterion);
        nbrManager.addCriterion(nbrCriterion,new AtomType[]{species1.getFactory().getType()});

        potentialwall = new P1WCAWall(space);
        potentialMaster.setSpecies(potentialwall, new Species[] { species });

        potentialwall1 = new P1WCAWall(space);
        potentialMaster.setSpecies(potentialwall1, new Species[] { species1 });

        potentialwallPorousA = new P1WCAPorousWall(space);
        potentialMaster.setSpecies(potentialwallPorousA, new Species[] { species });
        
        potentialwallPorousA1 = new P1WCAPorousWall(space);
        potentialMaster.setSpecies(potentialwallPorousA1, new Species[] { species1 });
        
        potentialwallPorousB = new P1WCAPorousWall(space);
        potentialMaster.setSpecies(potentialwallPorousB, new Species[] { species });
        
        potentialwallPorousB1 = new P1WCAPorousWall(space);
        potentialMaster.setSpecies(potentialwallPorousB1, new Species[] { species1 });


        species.setNMolecules(20);
        species1.setNMolecules(20);
        phase = new Phase(this);
        
        integratorDCV = new IntegratorDCVGCMD(potentialMaster, species,
                species1);
        Default.AUTO_REGISTER = false;
        final IntegratorVelocityVerlet integrator = new IntegratorVelocityVerlet(
                potentialMaster, space);
        final IntegratorMC integratorMC = new IntegratorMC(potentialMaster);
        Default.AUTO_REGISTER = true;
        integratorDCV.addPhase(phase);

        /***/
        integratorDCV.addListener(nbrManager);
        nbrManager.setRange(potential.getRange() * neighborRangeFac);
        integratorMC.addMCMoveListener(potentialMasterHybrid.getNbrCellManager(phase).makeMCMoveListener());
        potentialMasterHybrid.calculate(phase, new PotentialCalculationAgents(potentialMasterHybrid));


        activityIntegrate = new ActivityIntegrate(
                integratorDCV);
        getController().addAction(activityIntegrate);

        //make MC integrator next
        integratorDCV.setIntegrators(integratorMC, integrator);
        integratorDCV.setTemperature(Kelvin.UNIT.toSim(500.));
        integrator.setIsothermal(false);
        integrator.setMeterTemperature(new MeterTemperature(speciesTube));
        //integrator.setSleepPeriod(1);
        integrator.setTimeStep(0.01);
        //integrator.setInterval(10);
        activityIntegrate.setDoSleep(true);
        phase.setBoundary(new BoundaryRectangularSlit(space, 2));
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
        poreCenter.Ea1Tv1(0.5, phase.boundary().dimensions());
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
        double zA = (length + 0.05) * phase.boundary().dimensions().x(2);
        double zB = (1.0 - length - 0.05) * phase.boundary().dimensions().x(2);
        potentialwallPorousA.setZ(zA);
        potentialwallPorousA1.setZ(zA);
        potentialwallPorousB.setZ(zB);
        potentialwallPorousB1.setZ(zB);

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

        accumulator1 = new AccumulatorAverage();
        DataPump profile1pump = new DataPump(profile1, accumulator1);
        new IntervalActionAdapter(profile1pump, integratorDCV);

        accumulator2 = new AccumulatorAverage();
        DataPump profile2pump = new DataPump(profile2, accumulator2);
        new IntervalActionAdapter(profile2pump, integratorDCV);

//remove for nbrlist        integrator.addIntervalListener(new PhaseImposePbc(phase));
    } //End of constructor

} //End of DCVGCMD class
