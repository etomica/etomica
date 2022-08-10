/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.data.types.DataDouble;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.potential.*;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeField;
import etomica.simulation.Simulation;
import etomica.space.*;
import etomica.space1d.Space1D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.species.SpeciesManager;
import etomica.units.dimensions.Null;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import org.apache.commons.math3.analysis.function.Gaussian;

import java.util.ArrayList;
import java.util.List;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
public class SimQuantumAO extends Simulation {

    public SpeciesGeneral species;
    public Box box;

    public SimQuantumAO(Space space, int nBeads, double temperature, double omega, double k4) {
        super(space);
        SpeciesGeneral species = new SpeciesBuilder(space)
                .addCount(AtomType.simpleFromSim(this), nBeads)
                .build();
        addSpecies(species);
        box = this.makeBox(new BoundaryRectangularNonperiodic(space));
        box.setNMolecules(species, 1);
        box.getBoundary().setBoxSize(Vector.of(new double[]{1}));
        double mass = nBeads*species.getMass();
        //pm1 for sampling the x^4 quartic contribution
        PotentialComputeField pm1 = new PotentialComputeField(getSpeciesManager(), box);
        P1Anharmonic pQuartic = new P1Anharmonic(space, 0, k4);
        pm1.setFieldPotential(species.getLeafType(), pQuartic);

        //pm2 that uses the full PI potential, for data collection
        //spring P2 part (x_i-x_{i+1})^2
        PotentialMasterBonding pmBonding = new PotentialMasterBonding(getSpeciesManager(), box);
        double hbar = 1.0;
        double k2_kin = mass*temperature*temperature/(hbar*hbar);
        P2Harmonic p2Bond = new P2Harmonic(k2_kin, 0);
        List<int[]> pairs = new ArrayList<>();
        for (int i=0; i<nBeads; i++) {
            int[] p = new int[]{i,i+1};
            for (int j=0; j<p.length; j++) {
                if (p[j] >= nBeads) p[j] -= nBeads;
            }
            pairs.add(p);
        }
        pmBonding.setBondingPotentialPair(species, p2Bond, pairs);

        //potential p1 part
        double k2 = mass*omega*omega;
        P1Anharmonic p1ah = new P1Anharmonic(space, k2, k4);
        PotentialComputeField pcP1 = new PotentialComputeField(getSpeciesManager(), box);
        pcP1.setFieldPotential(species.getLeafType(), p1ah);

        //TOTAL
        PotentialCompute pm2 = new PotentialComputeAggregate(pmBonding, pcP1);


        IntegratorMC integratorMC = new IntegratorMC(pm1, random, temperature, box);
        MCMoveHO  atomMove = new MCMoveHO(space, random, temperature, omega);

//        atomMove.setLambdas(k_gauss);
//        atomMove.setTemperature(temperature);
        integratorMC.getMoveManager().addMCMove(atomMove);

        getController().addActivity(new ActivityIntegrate(integratorMC));



        // make an IntegratorMC and write a move with a doTrial that samples the harmonic part of your system


    }

    public static void main(String[] args) {

        OctaneParams params = new OctaneParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.numSteps = 35000000;
            params.graphics = false;
        }

        double temperature = params.temperature;
        double omega = params.omega;;
        double k4 = params.k4;;
        int nBeads = params.nBeads;
        boolean graphics = params.graphics;
        long numSteps = params.numSteps;
        System.out.println(numSteps+" steps");

        final SimQuantumAO sim = new SimQuantumAO(Space1D.getInstance(), nBeads, temperature, omega, k4);


    }

    public static class OctaneParams extends ParameterBase {
        public double temperature = 1;
        public int nBeads = 10;
        public boolean graphics = false;
        public double omega = 1.0; // m*w^2
        public double k4 = 1.0;
        public long numSteps = 200000;
    }
}
