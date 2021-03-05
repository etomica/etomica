/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.dcvgcmd;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomType;
import etomica.atom.IAtomList;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.box.Box;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceGroup;
import etomica.data.meter.*;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.CriterionPositionWall;
import etomica.nbr.PotentialMasterHybrid;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.units.Kelvin;
import etomica.util.random.RandomMersenneTwister;

/**
 * Dual-control-volume grand-canonical molecular dynamics simulation.
 * TraPPE propane and propene
 */
public class DCVGCMD extends Simulation {

    public IntegratorDCVGCMD integratorDCV;
    public P1WCAWall potentialwall;
    public P1WCAWall potentialwall1;
    public SpeciesGeneral propane;
    public SpeciesGeneral propene;
    public SpeciesGeneral membrane;
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
    public DataPumpListener profile1pump, profile2pump;
    public MeterMoleculeTemperature temperature1;
    public MeterMoleculeTemperature temperature2;
    public MeterProfileByAtoms profileTemperature1;
    public MeterProfileByAtoms profileTemperature2;
    public DataPumpListener profile1TemperaturePump, profile2TemperaturePump;
    public Vector poreCenter;
    public P2LennardJones p2MM;
    public double sigma;

    //Constructor
    public DCVGCMD() {
        this(Space3D.getInstance());
    }

    private DCVGCMD(Space _space) {
        //Instantiate classes
        super(_space);

        setRandom(new RandomMersenneTwister(1));
        double massC = Carbon.INSTANCE.getMass();
        double massH = Hydrogen.INSTANCE.getMass();
        ElementSimple CH3 = new ElementSimple("CH3", massC + 3 * massH);
        ElementSimple CH2 = new ElementSimple("CH2", massC + 2 * massH);
        ElementSimple CH = new ElementSimple("CH", massC + massH);
        AtomType propaneCH3 = new AtomType(CH3, "propaneCH3");
        AtomType propaneCH2 = new AtomType(CH2, "propaneCH2");
        double thetaPropane = 114 * Math.PI / 180;
        propane = new SpeciesBuilder(space)
                .setDynamic(true)
                .addAtom(propaneCH3, Vector.of(-1.54, 0, 0))
                .addAtom(propaneCH2, Vector.of(0, 0, 0))
                .addAtom(propaneCH3, Vector.of(1.54 * Math.cos(Math.PI - thetaPropane), 1.54 * Math.sin(thetaPropane), 0))
                .build();

        AtomType propeneCH3 = new AtomType(CH3, "propeneCH3");
        AtomType propeneCH = new AtomType(CH, "propeneCH");
        AtomType propeneCH2 = new AtomType(CH2, "propeneCH2");
        double thetaPropene = 119.7 * Math.PI / 180;
        propene = new SpeciesBuilder(space)
                .setDynamic(true)
                .addAtom(propeneCH3, Vector.of(-1.54, 0, 0))
                .addAtom(propeneCH, Vector.of(0, 0, 0))
                .addAtom(propeneCH2, Vector.of(1.33 * Math.cos(Math.PI - thetaPropene), 1.33 * Math.sin(thetaPropene), 0))
                .build();

        boolean membraneFixed = false;
        AtomType atomTypeMembrane = AtomType.simple("M", membraneFixed ? Double.POSITIVE_INFINITY : 12);
        membrane = SpeciesGeneral.monatomic(space, atomTypeMembrane, true);

        addSpecies(propane);
        addSpecies(propene);
        addSpecies(membrane);

        box = this.makeBox(new BoundaryRectangularSlit(2, space));
        double Lxy = 40;
        double Lz = 100;
        double zFraction = 0.1;
        box.getBoundary().setBoxSize(new Vector3D(Lxy, Lxy, Lz));

        PotentialMasterHybrid potentialMaster = new PotentialMasterHybrid(this, 5.2, space);
        sigma = 2.0;
        double epsilon = 119.8;
        //Default.makeLJDefaults();
        //Default.BOX_SIZE = 14.0;

        double neighborRangeFac = 1.2;
        potentialMaster.setCellRange(1);
        double rc = 11;
        potentialMaster.setRange(neighborRangeFac * rc);

        // propane
        P2LennardJones p2CH3CH3 = new P2LennardJones(space, 3.75, 98);
        P2LennardJones p2CH3CH2 = new P2LennardJones(space, (3.75 + 3.95) / 2, Math.sqrt(98 * 46));
        P2LennardJones p2CH2CH2 = new P2LennardJones(space, 3.95, 46);

        // propene
        P2LennardJones p2CH2eCH2e = new P2LennardJones(space, 3.675, 85);
        P2LennardJones p2CHeCHe = new P2LennardJones(space, 3.73, 47);
        P2LennardJones p2CH3CH2e = new P2LennardJones(space, (3.75 + 3.675) / 2, Math.sqrt(85 * 98));
        P2LennardJones p2CH3CHe = new P2LennardJones(space, (3.75 + 3.73) / 2, Math.sqrt(47 * 98));
        P2LennardJones p2CH2eCHe = new P2LennardJones(space, (3.675 + 3.73) / 2, Math.sqrt(47 * 85));

        // propane - propene
        P2LennardJones p2CH2CHe = new P2LennardJones(space, (3.675 + 3.73) / 2, Math.sqrt(85 * 47));
        P2LennardJones p2CH2CH2e = new P2LennardJones(space, (3.95 + 3.675) / 2, Math.sqrt(85 * 46));

        p2MM = new P2LennardJones(space, sigma, epsilon);
        P2LennardJones p2MCH3 = new P2LennardJones(space, (sigma + 3.75) / 2, Math.sqrt(epsilon * 98));
        P2LennardJones p2MCH2 = new P2LennardJones(space, (sigma + 3.95) / 2, Math.sqrt(epsilon * 46));
        P2LennardJones p2MCH2e = new P2LennardJones(space, (sigma + 3.675) / 2, Math.sqrt(epsilon * 85));
        P2LennardJones p2MCHe = new P2LennardJones(space, (sigma + 3.73) / 2, Math.sqrt(epsilon * 47));

        P2SoftSphericalTruncatedForceShifted p2CH3CH3sf = new P2SoftSphericalTruncatedForceShifted(space, p2CH3CH3, rc);
        P2SoftSphericalTruncatedForceShifted p2CH3CH2sf = new P2SoftSphericalTruncatedForceShifted(space, p2CH3CH2, rc);
        P2SoftSphericalTruncatedForceShifted p2CH2CH2sf = new P2SoftSphericalTruncatedForceShifted(space, p2CH2CH2, rc);
        P2SoftSphericalTruncatedForceShifted p2CH2eCH2esf = new P2SoftSphericalTruncatedForceShifted(space, p2CH2eCH2e, rc);
        P2SoftSphericalTruncatedForceShifted p2CHeCHesf = new P2SoftSphericalTruncatedForceShifted(space, p2CHeCHe, rc);
        P2SoftSphericalTruncatedForceShifted p2CH3CH2esf = new P2SoftSphericalTruncatedForceShifted(space, p2CH3CH2e, rc);
        P2SoftSphericalTruncatedForceShifted p2CH3CHesf = new P2SoftSphericalTruncatedForceShifted(space, p2CH3CHe, rc);
        P2SoftSphericalTruncatedForceShifted p2CH2eCHesf = new P2SoftSphericalTruncatedForceShifted(space, p2CH2eCHe, rc);
        P2SoftSphericalTruncatedForceShifted p2CH2CHesf = new P2SoftSphericalTruncatedForceShifted(space, p2CH2CHe, rc);
        P2SoftSphericalTruncatedForceShifted p2CH2CH2esf = new P2SoftSphericalTruncatedForceShifted(space, p2CH2CH2e, rc);

        P2SoftSphericalTruncatedForceShifted p2MMsf = new P2SoftSphericalTruncatedForceShifted(space, p2MM, rc);
        P2SoftSphericalTruncatedForceShifted p2MCH3sf = new P2SoftSphericalTruncatedForceShifted(space, p2MCH3, rc);
        P2SoftSphericalTruncatedForceShifted p2MCH2sf = new P2SoftSphericalTruncatedForceShifted(space, p2MCH2, rc);
        P2SoftSphericalTruncatedForceShifted p2MCH2esf = new P2SoftSphericalTruncatedForceShifted(space, p2MCH2e, rc);
        P2SoftSphericalTruncatedForceShifted p2MCHesf = new P2SoftSphericalTruncatedForceShifted(space, p2MCHe, rc);

        potentialMaster.addPotential(p2CH3CH3sf, new AtomType[]{propaneCH3, propaneCH3});
        potentialMaster.addPotential(p2CH3CH2sf, new AtomType[]{propaneCH3, propaneCH2});
        potentialMaster.addPotential(p2CH2CH2sf, new AtomType[]{propaneCH2, propaneCH2});

        potentialMaster.addPotential(p2CH3CH3sf, new AtomType[]{propeneCH3, propeneCH3});
        potentialMaster.addPotential(p2CH2eCH2esf, new AtomType[]{propeneCH2, propeneCH2});
        potentialMaster.addPotential(p2CHeCHesf, new AtomType[]{propeneCH, propeneCH});
        potentialMaster.addPotential(p2CH3CH2esf, new AtomType[]{propeneCH3, propeneCH2});
        potentialMaster.addPotential(p2CH3CHesf, new AtomType[]{propeneCH3, propeneCH});
        potentialMaster.addPotential(p2CH2eCHesf, new AtomType[]{propeneCH2, propeneCH});

        potentialMaster.addPotential(p2CH2CHesf, new AtomType[]{propaneCH2, propeneCH});
        potentialMaster.addPotential(p2CH2CH2esf, new AtomType[]{propaneCH2, propeneCH2});
        potentialMaster.addPotential(p2CH3CH2esf, new AtomType[]{propaneCH2, propeneCH3});
        potentialMaster.addPotential(p2CH3CHesf, new AtomType[]{propaneCH3, propeneCH});
        potentialMaster.addPotential(p2CH3CH2esf, new AtomType[]{propaneCH3, propeneCH2});
        potentialMaster.addPotential(p2CH3CH3sf, new AtomType[]{propaneCH3, propeneCH3});

        if (!membraneFixed) potentialMaster.addPotential(p2MMsf, new AtomType[]{atomTypeMembrane, atomTypeMembrane});

        potentialMaster.addPotential(p2MCH3sf, new AtomType[]{atomTypeMembrane, propaneCH3});
        potentialMaster.addPotential(p2MCH2sf, new AtomType[]{atomTypeMembrane, propaneCH2});
        potentialMaster.addPotential(p2MCH3sf, new AtomType[]{atomTypeMembrane, propeneCH3});
        potentialMaster.addPotential(p2MCH2esf, new AtomType[]{atomTypeMembrane, propeneCH2});
        potentialMaster.addPotential(p2MCHesf, new AtomType[]{atomTypeMembrane, propeneCH});

        P1HarmonicSite p1Membrane = new P1HarmonicSite(space);
        p1Membrane.setSpringConstant(10000);

        P2Harmonic p2BondCH3 = new P2Harmonic(space, 1000000, 1.54);
        P2Harmonic p2BondCHCH2 = new P2Harmonic(space, 1000000, 1.33);
        PotentialGroup p1Propane = potentialMaster.makePotentialGroup(1);
        p1Propane.addPotential(p2BondCH3, new ApiIndexList(new int[][]{{0, 1}, {1, 2}}));
        PotentialGroup p1Propene = potentialMaster.makePotentialGroup(1);
        p1Propene.addPotential(p2BondCH3, new ApiIndexList(new int[][]{{0, 1}}));
        p1Propene.addPotential(p2BondCHCH2, new ApiIndexList(new int[][]{{1, 2}}));

        P3BondAngle p3Propane = new P3BondAngle(space);
        p3Propane.setAngle(thetaPropane);
        p3Propane.setEpsilon(Kelvin.UNIT.toSim(62500));
        P3BondAngle p3Propene = new P3BondAngle(space);
        p3Propene.setAngle(thetaPropene);
        p3Propene.setEpsilon(Kelvin.UNIT.toSim(70420));
        p1Propane.addPotential(p3Propane, new Atomset3IteratorIndexList(new int[][]{{0, 1, 2}}));
        p1Propene.addPotential(p3Propene, new Atomset3IteratorIndexList(new int[][]{{0, 1, 2}}));

        potentialMaster.addPotential(p1Propane, new ISpecies[]{propane});
        potentialMaster.addPotential(p1Propene, new ISpecies[]{propene});

        double neighborRangeFacHalf = (1.0 + neighborRangeFac) * 0.5;

        potentialwall = new P1WCAWall(space, 2, sigma, epsilon);
        CriterionPositionWall criterionWall = new CriterionPositionWall(this);
        double wallCut = potentialwall.getRange() + 0.01 * Lz;
        criterionWall.setInteractionRange(wallCut);
        criterionWall.setNeighborRange(neighborRangeFacHalf * wallCut);
        criterionWall.setWallDim(2);
        potentialMaster.addPotential(potentialwall, new AtomType[]{propane.getTypeByName("propaneCH2")});
        potentialMaster.getPotentialMasterList().setCriterion1Body(potentialwall, propaneCH2, criterionWall);

        potentialwall1 = new P1WCAWall(space, 2, sigma, epsilon);
        CriterionPositionWall criterionWall1 = new CriterionPositionWall(this);
        criterionWall1.setInteractionRange(wallCut);
        criterionWall1.setNeighborRange(neighborRangeFacHalf * wallCut);
        criterionWall1.setWallDim(2);
        potentialMaster.addPotential(potentialwall1, new AtomType[]{propene.getTypeByName("propeneCH")});
        potentialMaster.getPotentialMasterList().setCriterion1Body(potentialwall1, propeneCH, criterionWall1);

        double temperature = Kelvin.UNIT.toSim(500.);
        integratorDCV = new IntegratorDCVGCMD(potentialMaster, temperature,
                propane, propene, box);
        integratorDCV.zFraction = zFraction;
        final IntegratorVelocityVerlet integratorMD = new IntegratorVelocityVerlet(potentialMaster, this.getRandom(), 0.05, 1.0, box);
        final IntegratorMC integratorMC = new IntegratorMC(this.getRandom(), potentialMaster, box);

        potentialMaster.setRange(3.95 * neighborRangeFac);
        integratorMC.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());

        getController().addActivity(new ActivityIntegrate(integratorDCV));

        integratorDCV.setIntegrators(integratorMC, integratorMD, getRandom());
        integratorMD.setIsothermal(false);

        thermometer = new MeterTemperature(getSpeciesManager(), box, space.D());

        integratorMD.setMeterTemperature(thermometer);
        //integrator.setSleepPeriod(1);
        integratorMD.setTimeStep(0.001);
        //integrator.setInterval(10);
        integratorMD.getEventManager().addListener(potentialMaster.getNeighborManager(box));
        // Crystal crystal = new Crystal(new PrimitiveTetragonal(space, 20,
        // 40),new BasisMonatomic(3));

        int nMembraneCell = 4;
        int nMembrane = nMembraneCell * nMembraneCell * nMembraneCell * 4;
        Box membraneBox = makeBox(new BoundaryRectangularPeriodic(space, Lxy));
        membraneBox.setNMolecules(membrane, nMembrane);
        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(membraneBox);
        box.setNMolecules(membrane, nMembrane);
        IAtomList pretendAtoms = membraneBox.getLeafList();
        IAtomList realAtoms = box.getLeafList();
        for (int i = 0; i < pretendAtoms.size(); i++) {
            realAtoms.get(i).getPosition().E(pretendAtoms.get(i).getPosition());
        }

        if (!membraneFixed) {
            AtomLeafAgentManager<Vector> siteManager = new AtomLeafAgentManager<>(iAtom -> {
                Vector v = space.makeVector();
                v.E(iAtom.getPosition());
                return v;
            }, box);
            p1Membrane.setAtomAgentManager(box, siteManager);
            potentialMaster.addPotential(p1Membrane, new AtomType[]{membrane.getLeafType()});
            integratorDCV.setDoThermalizeMembrane(true, membrane, random);
        }

        MyMCMove[] moves = integratorDCV.mcMoves();
        moves[0].setSpecies(propane, temperature, p1Propane);
        moves[1].setSpecies(propane, temperature, p1Propane);
        moves[2].setSpecies(propene, temperature, p1Propene);
        moves[3].setSpecies(propene, temperature, p1Propene);
        meterFlux0 = new MeterFlux(moves[0], integratorDCV);
        meterFlux1 = new MeterFlux(moves[1], integratorDCV);
        meterFlux2 = new MeterFlux(moves[2], integratorDCV);
        meterFlux3 = new MeterFlux(moves[3], integratorDCV);
        meterFlux0.setBox(box);
        meterFlux1.setBox(box);
        meterFlux2.setBox(box);
        meterFlux3.setBox(box);
        fluxMeters = new DataSourceGroup(new MeterFlux[]{meterFlux0, meterFlux1, meterFlux2, meterFlux3});

        density1 = new MeterNMolecules();
        density2 = new MeterNMolecules();
        density1.setBox(box);
        density2.setBox(box);
        density1.setSpecies(propane);
        density2.setSpecies(propene);
        profile1 = new MeterProfileByVolume(space);
        profile2 = new MeterProfileByVolume(space);
        profile1.setDataSource(density1);
        profile1.setBox(box);
        profile1.setProfileDim(2);
        profile2.setDataSource(density2);
        profile2.setBox(box);
        profile2.setProfileDim(2);

        accumulator1 = new AccumulatorAverageFixed();
        profile1pump = new DataPumpListener(profile1, accumulator1);
        integratorDCV.getEventManager().addListener(profile1pump);

        accumulator2 = new AccumulatorAverageFixed();
        profile2pump = new DataPumpListener(profile2, accumulator2);
        integratorDCV.getEventManager().addListener(profile2pump);

        temperature1 = new MeterMoleculeTemperature();
        temperature2 = new MeterMoleculeTemperature();
        temperature1.setBox(box);
        temperature2.setBox(box);
        profileTemperature1 = new MeterProfileByAtoms(space);
        profileTemperature2 = new MeterProfileByAtoms(space);
        profileTemperature1.setSpecies(propane);
        profileTemperature2.setSpecies(propene);
        profileTemperature1.setDataSource(temperature1);
        profileTemperature1.setBox(box);
        profileTemperature1.setProfileDim(2);
        profileTemperature2.setDataSource(temperature2);
        profileTemperature2.setBox(box);
        profileTemperature2.setProfileDim(2);

        profile1TemperaturePump = new DataPumpListener(profileTemperature1, null);
        integratorDCV.getEventManager().addListener(profile1TemperaturePump);

        profile2TemperaturePump = new DataPumpListener(profileTemperature2, null);
        integratorDCV.getEventManager().addListener(profile2TemperaturePump);

        potentialMaster.getNbrCellManager(box).assignCellAll();
    }

    public static void main(String[] args) {
        DCVGCMD sim = new DCVGCMD();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorDCV, 5000000));
    }
}
