/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedDensity.positionOrientation;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.config.ConformationChainLinear;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.dcvgcmd.P1WCAWall;
import etomica.graphics.ColorSchemeRandomByMolecule;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveMoleculeRotate;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.BondingInfo;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedForceShifted;
import etomica.potential.PotentialMaster;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeField;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.units.Degree;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.util.Arrays;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
public class LJMCDimer extends Simulation {

    public IntegratorMC integrator;
    public SpeciesGeneral species;
    public Box box;
    public P2LennardJones potential;

    public LJMCDimer(int D, int numMolecules, double density, double temperature) {
        super(Space.getInstance(D));

        double a = Degree.UNIT.toSim(45);
        double[] a2 = new double[D-1];
        Arrays.fill(a2, a);
        species = new SpeciesBuilder(space)
                .addCount(AtomType.simpleFromSim(this), 2)
                .setDynamic(true)
                .withConformation(new ConformationChainLinear(space, 1, a2))
                .build();
        addSpecies(species);

        double sigma = 1.0;
        box = this.makeBox(new BoundaryRectangularSlit(2, getSpace()));
        PotentialMaster potentialMaster = new PotentialMasterCell(this.getSpeciesManager(), box, 2, BondingInfo.noBonding());
        P1WCAWall wall = new P1WCAWall(box, 1, 1);
        PotentialComputeField pcField = new PotentialComputeField(getSpeciesManager(), box);
        AtomType leafType = species.getLeafType();
        pcField.setFieldPotential(leafType, wall);
        PotentialComputeAggregate pcAgg = new PotentialComputeAggregate(potentialMaster, pcField);

        integrator = new IntegratorMC(pcAgg, this.getRandom(), 1.0, box);
        integrator.setTemperature(temperature);

        integrator.getMoveManager().addMCMove(new MCMoveMolecule(this.getRandom(), pcAgg, box));
        integrator.getMoveManager().addMCMove(new MCMoveMoleculeRotate(this.getRandom(), pcAgg, box));

        box.setNMolecules(species, numMolecules);
        new BoxInflate(box, space, 1.2*density).actionPerformed();

        potential = new P2LennardJones(sigma, 1.0);
        P2SoftSphericalTruncatedForceShifted p2 = new P2SoftSphericalTruncatedForceShifted(potential, 3.0);
        potentialMaster.setPairPotential(leafType, leafType, p2);

        ConfigurationLattice configuration = new ConfigurationLattice(D==3 ? new LatticeCubicFcc(space) : new LatticeOrthorhombicHexagonal(space), space);
        configuration.initializeCoordinates(box);
        Vector l = box.getBoundary().getBoxSize();
        l.setX(2, l.getX(2)*1.2);
        box.getBoundary().setBoxSize(l);
    }

    public static void main(String[] args) {
        DimerParams params = new DimerParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            // override defaults here
            params.D = 2;
        }

        int numMolecules = params.numMolecules;
        double temperature = params.temperature;
        double density = params.density;
        int D = params.D;

        final LJMCDimer sim = new LJMCDimer(D, numMolecules, density, temperature);

        boolean graphics = true;

        MeterPotentialEnergyFromIntegrator energy = new MeterPotentialEnergyFromIntegrator(sim.integrator);

        if (graphics) {
            final String APP_NAME = "LJMDDimer";
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME, 3);

            AccumulatorAverageCollapsing avgEnergy = new AccumulatorAverageCollapsing();
            avgEnergy.setPushInterval(10);
            DataPumpListener pump = new DataPumpListener(energy, avgEnergy, numMolecules);
            sim.integrator.getEventManager().addListener(pump);

            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
            simGraphic.getController().getDataStreamPumps().add(pump);
            simGraphic.getDisplayBox(sim.box).setColorScheme(new ColorSchemeRandomByMolecule(sim.getSpeciesManager(), sim.box, sim.getRandom()));

            simGraphic.makeAndDisplayFrame(APP_NAME);

            DisplayTextBoxesCAE display = new DisplayTextBoxesCAE();
            display.setAccumulator(avgEnergy);
            simGraphic.add(display);
            return;
        }
    }

    public static class DimerParams extends ParameterBase {
        public int numMolecules = 256;
        public double temperature = 2;
        public double density = 0.4;
        public int D = 3;
    }
}
