/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.tests;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.config.ConformationChainLinear;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.graphics.ColorSchemeRandomByMolecule;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.BondingInfo;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedForceShifted;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.units.Degree;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
public class TestLJMCDimerRigid extends Simulation {

    public IntegratorMC integrator;
    public SpeciesGeneral species;
    public Box box;
    public P2LennardJones potential;
    public MeterPotentialEnergyFromIntegrator energy;
    public AccumulatorAverageCollapsing avgEnergy;
    public DataPump pump;


    public TestLJMCDimerRigid(int moleculeSize, int totalAtoms, boolean cellListing) {
        super(Space3D.getInstance());

        species = new SpeciesBuilder(space)
                .addCount(AtomType.simpleFromSim(this), moleculeSize)
                .setDynamic(true)
                .withConformation(new ConformationChainLinear(space, 1, new double[]{Degree.UNIT.toSim(45), Degree.UNIT.toSim(45), 0}))
                .build();
        addSpecies(species);

        double sigma = 1.0;
        box = this.makeBox();
        PotentialMaster potentialMaster = cellListing
                ? new PotentialMasterCell(this.getSpeciesManager(), box, 2, BondingInfo.noBonding())
                : new PotentialMaster(this.getSpeciesManager(), box, BondingInfo.noBonding());

        integrator = new IntegratorMC(potentialMaster, this.getRandom(), 1.0, box);
        integrator.setTemperature(moleculeSize);
        integrator.setIsothermal(true);

        integrator.getMoveManager().addMCMove(new MCMoveMolecule(this.getRandom(), potentialMaster, box));

        box.setNMolecules(species, totalAtoms / moleculeSize);
        new BoxInflate(box, space, 0.9 / moleculeSize).actionPerformed();
        System.out.println("box size: "+box.getBoundary().getBoxSize());

        potential = new P2LennardJones(space, sigma, 1.0);
        AtomType leafType = species.getLeafType();
        P2SoftSphericalTruncatedForceShifted p2 = new P2SoftSphericalTruncatedForceShifted(space, potential, 3.0);
        potentialMaster.setPairPotential(leafType, leafType, p2);

        if (!cellListing) {
            BoxImposePbc imposepbc = new BoxImposePbc(space);
            imposepbc.setBox(box);
            integrator.getEventManager().addListener(new IntegratorListenerAction(imposepbc));
        }

        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box);
        integrator.reset();
        energy = new MeterPotentialEnergyFromIntegrator(integrator);
        System.out.println("u0: "+energy.getDataAsScalar());
        avgEnergy = new AccumulatorAverageCollapsing();
        avgEnergy.setPushInterval(10);
        pump = new DataPump(energy, avgEnergy);
        IntegratorListenerAction pumpListener = new IntegratorListenerAction(pump);
        pumpListener.setInterval(totalAtoms);
        integrator.getEventManager().addListener(pumpListener);
    }

    public static void main(String[] args) {
        final String APP_NAME = "LJMDDimer";
        final TestLJMCDimerRigid sim = new TestLJMCDimerRigid(2, 512, true);
//        long t0 = System.nanoTime();
//        sim.getController().actionPerformed();
//        long t1 = System.nanoTime();
//        System.out.println((t1 - t0) / 1e6);

        sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
        final SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME, 3);


        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
        simGraphic.getController().getDataStreamPumps().add(sim.pump);
        simGraphic.getDisplayBox(sim.box).setColorScheme(new ColorSchemeRandomByMolecule(sim.getSpeciesManager(), sim.box, sim.getRandom()));

        simGraphic.makeAndDisplayFrame(APP_NAME);

        DisplayTextBoxesCAE display = new DisplayTextBoxesCAE();
        display.setAccumulator(sim.avgEnergy);
        simGraphic.add(display);
    }
}
