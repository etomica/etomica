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
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.BondingInfo;
import etomica.potential.P2WCA;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeField;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.units.Kelvin;

/**
 * Dual-control-volume grand-canonical molecular dynamics simulation.
 */
public class DCVGCMD extends Simulation {

    public IntegratorDCVGCMD integratorDCV;
    public P2WCA potential;
    public P2WCA potential1;
    public P1WCAPorousWall potentialwallPorous;
    public SpeciesGeneral species1;
    public SpeciesGeneral species2;
    public SpeciesGeneral speciesTube;
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

    public ConfigurationLatticeTube config;

    //Constructor
    public DCVGCMD() {
        this(Space3D.getInstance());
    }

    private DCVGCMD(Space _space) {
        //Instantiate classes
        super(_space);

        species1 = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        species2 = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        speciesTube = new SpeciesBuilder(space)
                .addCount(AtomType.simple("T", Double.POSITIVE_INFINITY), 20 * 40)
                .withConformation(new ConformationTube(space, 20))
                .setDynamic(true)
                .build();
        addSpecies(species1);
        addSpecies(species2);
        addSpecies(speciesTube);

        box = this.makeBox(new BoundaryRectangularSlit(2, space));
        box.getBoundary().setBoxSize(new Vector3D(40, 40, 80));

        double neighborRangeFac = 1.4;
        double mass = 40.;
        double sigma = 3.0;
        double epsilon = 119.8;
        potential = new P2WCA(sigma, epsilon);
        PotentialMasterCell potentialMasterCell = new PotentialMasterCell(getSpeciesManager(), box, 1, BondingInfo.noBonding());
        PotentialMasterList potentialMasterList = new PotentialMasterList(getSpeciesManager(), box, 1, potential.getRange() * neighborRangeFac, BondingInfo.noBonding());
        PotentialComputeField pcField = new PotentialComputeField(getSpeciesManager(), box);

        AtomType tubetype = speciesTube.getAtomType(0);
        AtomType speciestype = species1.getLeafType();
        AtomType speciestype1 = species2.getLeafType();
        ((ElementSimple) speciestype.getElement()).setMass(mass);
        ((ElementSimple) speciestype1.getElement()).setMass(mass);

        //0-0 intraspecies interaction
        potentialMasterCell.setPairPotential(species1.getLeafType(), species1.getLeafType(), potential);
        potentialMasterList.setPairPotential(species1.getLeafType(), species1.getLeafType(), potential);

        //1-1 intraspecies interaction
        P2WCA potential11 = new P2WCA(sigma, epsilon);
        potentialMasterCell.setPairPotential(species2.getLeafType(), species2.getLeafType(), potential11);
        potentialMasterList.setPairPotential(species2.getLeafType(), species2.getLeafType(), potential11);

        //0-1 interspecies interaction
        potential1 = new P2WCA(sigma, epsilon);
        potentialMasterCell.setPairPotential(species2.getLeafType(), species1.getLeafType(), potential1);
        potentialMasterList.setPairPotential(species2.getLeafType(), species1.getLeafType(), potential1);

        P2WCA potentialTubeAtom = new P2WCA(sigma, epsilon);
        potentialMasterCell.setPairPotential(tubetype, speciestype, potentialTubeAtom);
        potentialMasterList.setPairPotential(tubetype, speciestype, potentialTubeAtom);

        P2WCA potentialTubeAtom1 = new P2WCA(sigma, epsilon);
        potentialMasterCell.setPairPotential(tubetype, speciestype1, potentialTubeAtom1);
        potentialMasterList.setPairPotential(tubetype, speciestype1, potentialTubeAtom1);

        double poreRadius = 1.05 * ((ConformationTube) speciesTube.getConformation()).tubeRadius;
        double length = 0.25;
        double z = 0.2 * box.getBoundary().getBoxSize().getX(2);
        potentialwallPorous = new P1WCAPorousWall(box, sigma, epsilon, poreRadius, z);
        pcField.setFieldPotential(species1.getLeafType(), potentialwallPorous);
        pcField.setFieldPotential(species2.getLeafType(), potentialwallPorous);

        PotentialComputeAggregate pcAggMC = new PotentialComputeAggregate(potentialMasterCell, pcField);
        PotentialComputeAggregate pcAggMD = new PotentialComputeAggregate(potentialMasterList, pcField);

        box.setNMolecules(species1, 20);
        box.setNMolecules(species2, 20);
        box.setNMolecules(speciesTube, 1);

        double temperature = Kelvin.UNIT.toSim(500.);
        integratorDCV = new IntegratorDCVGCMD(pcAggMC, temperature,
                species1, species2, box);
        final IntegratorVelocityVerlet integratorMD = new IntegratorVelocityVerlet(pcAggMD, this.getRandom(), 0.05, 1.0, box);
        final IntegratorMC integratorMC = new IntegratorMC(pcAggMC, this.getRandom(), 1.0, box);

        getController().addActivity(new ActivityIntegrate(integratorDCV));

        integratorDCV.setIntegrators(integratorMC, integratorMD, getRandom());
        integratorMD.setIsothermal(false);

        thermometer = new MeterTemperature(getSpeciesManager(), box, space.D());

        integratorMD.setMeterTemperature(thermometer);
        integratorMD.setTimeStep(0.005);
        //integrator.setInterval(10);
        config = new ConfigurationLatticeTube(new LatticeCubicFcc(space), length, space);
        config.setSpeciesSpheres(new SpeciesGeneral[]{species1, species2});
        config.setSpeciesTube(speciesTube);
        config.initializeCoordinates(box);

        MyMCMove[] moves = integratorDCV.mcMoves();
        meterFlux0 = new MeterFlux(moves[0], integratorMD);
        meterFlux1 = new MeterFlux(moves[1], integratorMD);
        meterFlux2 = new MeterFlux(moves[2], integratorMD);
        meterFlux3 = new MeterFlux(moves[3], integratorMD);
        meterFlux0.setBox(box);
        meterFlux1.setBox(box);
        meterFlux2.setBox(box);
        meterFlux3.setBox(box);
        fluxMeters = new DataSourceGroup(new MeterFlux[]{meterFlux0, meterFlux1, meterFlux2, meterFlux3});

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
    }

    public static void main(String[] args) {
        DCVGCMD sim = new DCVGCMD();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorDCV, 5000));
    }
}
