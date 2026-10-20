/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.zeno;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.config.IConformation;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMoleculeRotate;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
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
            params.numAtoms = 1;
            params.confFile = "sphere";
            params.steps = 1;
            params.walks = 1000000;
            params.rotation = false;
            params.eigen = true;
            params.sigma = 1;
        }

        System.out.println("config: "+params.confFile);
        System.out.println("steps: "+params.steps);
        System.out.println("walks: "+params.walks);
        System.out.println("rotation: "+(params.rotation ? "on": "off"));
        System.out.println("eigenvalues for q: "+(params.eigen ? "on" : "off"));
        Simulation sim = new Simulation(Space3D.getInstance());
        AtomType type1 = AtomType.simpleFromSim(sim);
        AtomType type2 = AtomType.simpleFromSim(sim);
        ISpecies species = new SpeciesBuilder(sim.getSpace())
                .addCount(type1, 1)
//                .addCount(type2, 1)
                .withConformation(new IConformation() {
                    @Override
                    public void initializePositions(IAtomList atomList) {}
                })
                .build();
        sim.addSpecies(species);
        sim.addBox(new Box(sim.getSpace()));
        sim.box().setNMolecules(species, 1);
        new ConfigurationFile(params.confFile).initializeCoordinates(sim.box());
        MCMoveMoleculeRotate rotate = new MCMoveMoleculeRotate(sim.getRandom(), new PotentialComputeAggregate(),sim.box());
        IntegratorMC integrator = new IntegratorMC(new PotentialComputeAggregate(),sim.getRandom(),1,sim.box());
        if (params.rotation) integrator.getMoveManager().addMCMove(rotate);
        rotate.setStepSize(Math.PI);
//        boundingSphereRadius = 36.660339;
        double[] sigma = new double[]{1};
        MeterZENO meterIntrinsicViscosity = new MeterZENO(sim.box(), sim.getRandom(), sigma, false, WalkerExterior.Shape.CUBE);
        WalkerExterior w = meterIntrinsicViscosity.walker;
        System.err.println("Need to force rotation about boundingSphereCenter");
//        rotate.r0.E(w.boundingSphereCenter);
        System.out.println("bounding sphere at "+w.boundingSphereCenter+" radius "+w.boundingSphereRadius);
//        meterIntrinsicViscosity.setBoundingSphere(CenterOfMass.position(sim.box(), sim.box().getMoleculeList().get(0)), 4);
        meterIntrinsicViscosity.setBoundingSphere(new Vector3D(0,0,0), 2);
        meterIntrinsicViscosity.setEigenvaluesForPade(params.eigen);
        meterIntrinsicViscosity.setNumWalks(params.walks);
        //System.out.println("boundingSphereRadius: " + meterIntrinsicViscosity.boundingSphereRadius);
        integrator.getEventManager().addListener(new IntegratorListenerAction(meterIntrinsicViscosity));
       // meterIntrinsicViscosity.actionPerformed();
        ActivityIntegrate ai =new ActivityIntegrate(integrator,params.steps);
        if (false) {
            sim.getController().addActivity(ai);
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "ZENO");
            simGraphic.makeAndDisplayFrame();
            return;
        }
        sim.getController().runActivityBlocking(ai);
        double viscosity = meterIntrinsicViscosity.getIntrinsicViscosity();
        double hydrodynamicRadius = meterIntrinsicViscosity.getHydrodynamicRadius();
        double vol = meterIntrinsicViscosity.computeVolume();
        System.out.println("volume: "+ vol);
        System.out.println("Intrinsic viscosity: " + viscosity);
        System.out.println("Hydrodynamic radius: " + hydrodynamicRadius);
        long t2 = System.nanoTime();
        System.out.println("time: " + (double)(t2 - t1) / (double)1.0E9F);
    }

    public static class ParametersZENOConfig extends ParameterBase {
        public int numAtoms;
        public String confFile;
        public long steps = 1;
        public long walks = 1000000;
        public boolean rotation;
        public boolean eigen;
        public double sigma = 1;
    }
}

