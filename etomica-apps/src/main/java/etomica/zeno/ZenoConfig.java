/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.zeno;

import etomica.atom.AtomType;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.config.IConformation;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Compute quantities that ZENO computes for a single config read from a file.
 */
public class ZenoConfig {

    public static void main(String[] args) {
        long t1 = System.nanoTime();
        ParametersZENOConfig params = new ParametersZENOConfig();

        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.numAtoms = 641;
            params.confFile = "star";
        }

        Simulation sim = new Simulation(Space3D.getInstance());
        AtomType type = AtomType.simpleFromSim(sim);
        ISpecies species = new SpeciesBuilder(sim.getSpace())
                .addCount(type, params.numAtoms)
                .withConformation(new IConformation() {
                    @Override
                    public void initializePositions(IAtomList atomList) {}
                })
                .build();
        sim.addSpecies(species);
        sim.addBox(new Box(sim.getSpace()));
        sim.box().setNMolecules(species, 1);
        new ConfigurationFile(params.confFile).initializeCoordinates(sim.box());

//        boundingSphereRadius = 36.660339;
        double[] sigma = new double[]{1};
        MeterZENO meterIntrinsicViscosity = new MeterZENO(sim.box(), sim.getRandom(), sigma);
        meterIntrinsicViscosity.setNumWalks(1000000L);
        meterIntrinsicViscosity.actionPerformed();
        double viscosity = meterIntrinsicViscosity.getIntrinsicViscosity();
        double hydrodynamicRadius = meterIntrinsicViscosity.getHydrodynamicRadius();
        System.out.println("Intrinsic viscosity: " + viscosity);
        System.out.println("Hydrodynamic radius: " + hydrodynamicRadius);
        long t2 = System.nanoTime();
        System.out.println("time: " + (double)(t2 - t1) / (double)1.0E9F);
    }

    public static class ParametersZENOConfig extends ParameterBase {
        public int numAtoms;
        public String confFile;
    }
}

