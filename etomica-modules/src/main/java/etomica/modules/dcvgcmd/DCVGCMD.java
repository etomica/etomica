/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.dcvgcmd;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPump;
import etomica.data.DataSourceGroup;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterProfileByVolume;
import etomica.data.meter.MeterTemperature;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.integrator.IntegratorListenerAction;
import etomica.nbr.CriterionPositionWall;
import etomica.nbr.CriterionType;
import etomica.nbr.PotentialMasterHybrid;
import etomica.potential.P2WCA;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
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
    public SpeciesSpheresMono species1;
    public SpeciesSpheresMono species2;
    public SpeciesTube speciesTube;
    public Box box;
    public DataSourceGroup fluxMeters;
    public MeterFlux meterFlux0, meterFlux1, meterFlux2, meterFlux3;
    public MeterTemperature thermometer;
    public MeterNMolecules density1;
    public MeterNMolecules density2;
    public MeterProfileByVolume profile1;
    public MeterProfileByVolume profile2;
    public AccumulatorAverage accumulator1;
    public AccumulatorAverage accumulator2;
    public DataPump profile1pump, profile2pump;
    public AccumulatorAverage fluxAccumulator;
    public Vector poreCenter;
    public ActivityIntegrate activityIntegrate;
    public ConfigurationLatticeTube config;
    
    //Constructor
    public DCVGCMD() {
        this(Space3D.getInstance());
    }

    private DCVGCMD(Space _space) {
        //Instantiate classes
        super(_space);
        PotentialMasterHybrid potentialMaster = new PotentialMasterHybrid(this, 5.2, space);
        double mass = 40.;
        double sigma = 3.0;
        double epsilon = 119.8;
        //Default.makeLJDefaults();
        //Default.BOX_SIZE = 14.0;

        species1 = new SpeciesSpheresMono(this, space);
        species1.setIsDynamic(true);
        species2 = new SpeciesSpheresMono(this, space);
        species2.setIsDynamic(true);
        speciesTube = new SpeciesTube(20, 40, space);
        addSpecies(species1);
        addSpecies(species2);
        addSpecies(speciesTube);
        AtomType tubetype = speciesTube.getLeafType();
        AtomType speciestype = species1.getLeafType();
        AtomType speciestype1 = species2.getLeafType();
        ((ElementSimple)speciestype.getElement()).setMass(mass);
        ((ElementSimple)speciestype1.getElement()).setMass(mass);
        
        double neighborRangeFac = 1.4;
        potentialMaster.setCellRange(1);
        potentialMaster.setRange(neighborRangeFac * sigma);

        //0-0 intraspecies interaction
        potential = new P2WCA(space, sigma, epsilon);
        potentialMaster.addPotential(potential, new AtomType[]{species1.getLeafType(), species1.getLeafType()});
        
        //1-1 intraspecies interaction
        P2WCA potential11 = new P2WCA(space, sigma, epsilon);
        potentialMaster.addPotential(potential11, new AtomType[]{species2.getLeafType(), species2.getLeafType()});

        //0-1 interspecies interaction
        potential1 = new P2WCA(space, sigma, epsilon);
        potentialMaster.addPotential(potential1, new AtomType[]{species2.getLeafType(), species1.getLeafType()});

        P2WCA potentialTubeAtom = new P2WCA(space, sigma, epsilon);
        potentialMaster.addPotential(potentialTubeAtom, new AtomType[]{tubetype, speciestype});
        
        P2WCA potentialTubeAtom1 = new P2WCA(space, sigma, epsilon);
        potentialMaster.addPotential(potentialTubeAtom1, new AtomType[]{tubetype, speciestype1});

        double neighborRangeFacHalf = (1.0+neighborRangeFac)*0.5;
        
        potentialwall = new P1WCAWall(space, 2, sigma, epsilon);
        CriterionPositionWall criterionWall = new CriterionPositionWall(this);
        criterionWall.setInteractionRange(potentialwall.getRange());
        criterionWall.setNeighborRange(neighborRangeFacHalf*potentialwall.getRange());
        criterionWall.setWallDim(2);
        potentialMaster.addPotential(potentialwall, new AtomType[]{species1.getLeafType()});
        potentialMaster.getPotentialMasterList().setCriterion(potentialwall, new CriterionType(criterionWall, speciestype));
        
        potentialwall1 = new P1WCAWall(space, 2, sigma, epsilon);
        CriterionPositionWall criterionWall1 = new CriterionPositionWall(this);
        criterionWall1.setInteractionRange(potentialwall1.getRange());
        criterionWall1.setNeighborRange(neighborRangeFacHalf*potentialwall1.getRange());
        criterionWall1.setWallDim(2);
        potentialMaster.addPotential(potentialwall1, new AtomType[]{species2.getLeafType()});
        potentialMaster.getPotentialMasterList().setCriterion(potentialwall1, new CriterionType(criterionWall1, speciestype1));

        potentialwallPorousA = new P1WCAPorousWall(space, sigma, epsilon);
        CriterionPositionWall criterionWallA = new CriterionPositionWall(this);
        criterionWallA.setInteractionRange(potentialwallPorousA.getRange());
        criterionWallA.setNeighborRange(neighborRangeFacHalf*potentialwallPorousA.getRange());
        criterionWallA.setWallDim(2);
        potentialMaster.addPotential(potentialwallPorousA, new AtomType[]{species1.getLeafType()});
        potentialMaster.getPotentialMasterList().setCriterion(potentialwallPorousA, new CriterionType(criterionWallA, speciestype));
        
        potentialwallPorousA1 = new P1WCAPorousWall(space, sigma, epsilon);
        CriterionPositionWall criterionWallA1 = new CriterionPositionWall(this);
        criterionWallA1.setInteractionRange(potentialwallPorousA1.getRange());
        criterionWallA1.setNeighborRange(neighborRangeFacHalf*potentialwallPorousA1.getRange());
        criterionWallA1.setWallDim(2);
        potentialMaster.addPotential(potentialwallPorousA1, new AtomType[]{species2.getLeafType()});
        potentialMaster.getPotentialMasterList().setCriterion(potentialwallPorousA1, new CriterionType(criterionWallA1, speciestype1));
        
        potentialwallPorousB = new P1WCAPorousWall(space, sigma, epsilon);
        CriterionPositionWall criterionWallB = new CriterionPositionWall(this);
        criterionWallB.setInteractionRange(potentialwallPorousB.getRange());
        criterionWallB.setNeighborRange(neighborRangeFacHalf*potentialwallPorousB.getRange());
        criterionWallB.setWallDim(2);
        potentialMaster.addPotential(potentialwallPorousB, new AtomType[]{species1.getLeafType()});
        potentialMaster.getPotentialMasterList().setCriterion(potentialwallPorousB, new CriterionType(criterionWallB, speciestype));
        
        potentialwallPorousB1 = new P1WCAPorousWall(space, sigma, epsilon);
        CriterionPositionWall criterionWallB1 = new CriterionPositionWall(this);
        criterionWallB1.setInteractionRange(potentialwallPorousB.getRange());
        criterionWallB1.setNeighborRange(neighborRangeFacHalf*potentialwallPorousB.getRange());
        criterionWallB1.setWallDim(2);
        potentialMaster.addPotential(potentialwallPorousB1, new AtomType[]{species2.getLeafType()});
        potentialMaster.getPotentialMasterList().setCriterion(potentialwallPorousB1, new CriterionType(criterionWallB1, speciestype1));


        box = new Box(space);
        box.setBoundary(new BoundaryRectangularSlit(2, space));
        box.getBoundary().setBoxSize(new Vector3D(40, 40, 80));
        addBox(box);
        box.setNMolecules(species1, 20);
        box.setNMolecules(species2, 20);
        box.setNMolecules(speciesTube, 1);
        
        double temperature = Kelvin.UNIT.toSim(500.);
        integratorDCV = new IntegratorDCVGCMD(potentialMaster, temperature, space,
        		                              species1, species2);
        final IntegratorVelocityVerlet integratorMD = new IntegratorVelocityVerlet(this, potentialMaster, space);
        final IntegratorMC integratorMC = new IntegratorMC(this, potentialMaster);
        integratorDCV.setBox(box);

        /***/
        potentialMaster.setRange(potential.getRange() * neighborRangeFac);
        integratorMC.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());


        activityIntegrate = new ActivityIntegrate(integratorDCV);
        getController().addAction(activityIntegrate);

        integratorDCV.setIntegrators(integratorMC, integratorMD, getRandom());
        integratorMD.setIsothermal(false);

        thermometer = new MeterTemperature(this, box, space.D());

        integratorMD.setMeterTemperature(thermometer);
        //integrator.setSleepPeriod(1);
        integratorMD.setTimeStep(0.005);
        //integrator.setInterval(10);
        integratorMD.getEventManager().addListener(potentialMaster.getNeighborManager(box));
        // Crystal crystal = new Crystal(new PrimitiveTetragonal(space, 20,
        // 40),new BasisMonatomic(3));
        double length = 0.25;
        config = new ConfigurationLatticeTube(new LatticeCubicFcc(space), length, space);
        config.setSpeciesSpheres(new SpeciesSpheresMono[]{species1, species2});
        config.setSpeciesTube(speciesTube);
        config.initializeCoordinates(box);

        //position of hole in porous-wall potential
        poreCenter = space.makeVector();
//        poreCenter.Ea1Tv1(0.5, box.getBoundary().getBoxSize());
        Vector[] poreCentersVector = new Vector[] { poreCenter };
        potentialwallPorousA.setPoreCenters(poreCentersVector);
        potentialwallPorousA1.setPoreCenters(poreCentersVector);
        potentialwallPorousB.setPoreCenters(poreCentersVector);
        potentialwallPorousB1.setPoreCenters(poreCentersVector);

        //radius of hole in porous-wall potential
        double poreRadius = 1.05 * ((ConformationTube)speciesTube.getConformation()).tubeRadius;
        potentialwallPorousA.setPoreRadius(poreRadius);
        potentialwallPorousA1.setPoreRadius(poreRadius);
        potentialwallPorousB.setPoreRadius(poreRadius);
        potentialwallPorousB1.setPoreRadius(poreRadius);

        //place porous-wall potentials; put just past the edges of the tube
        double zA = (-0.5 + length + 0.05) * box.getBoundary().getBoxSize().getX(2);
        double zB = ( 0.5 - length - 0.05) * box.getBoundary().getBoxSize().getX(2);
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

        density1 = new MeterNMolecules();
        density2 = new MeterNMolecules();
        density1.setBox(box);
        density2.setBox(box);
        density1.setSpecies(species1);
        density2.setSpecies(species2);
        profile1 = new MeterProfileByVolume(space);
        profile2 = new MeterProfileByVolume(space);
        profile1.setDataSource(density1);
        profile1.setBox(box);
        profile1.setProfileDim(2);
        profile2.setDataSource(density2);
        profile2.setBox(box);
        profile2.setProfileDim(2);

        accumulator1 = new AccumulatorAverageFixed();
        profile1pump = new DataPump(profile1, accumulator1);
        integratorDCV.getEventManager().addListener(new IntegratorListenerAction(profile1pump));

        accumulator2 = new AccumulatorAverageFixed();
        profile2pump = new DataPump(profile2, accumulator2);
        integratorDCV.getEventManager().addListener(new IntegratorListenerAction(profile2pump));

        potentialMaster.getNbrCellManager(box).assignCellAll();
    }
}
