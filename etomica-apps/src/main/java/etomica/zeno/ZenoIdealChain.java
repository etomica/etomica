/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.zeno;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.config.IConformation;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Compute quantities that ZENO computes for an ideal chain.
 * See <a href="https://dx.doi.org/10.1103/PhysRevE.81.021803">Influence
 * of variable hydrodynamic interaction strength on the transport properties
 * of coiled polymers</a>
 */
public class ZenoIdealChain {
    // for N=10, fig says Rh=3.98(4), eta=2.05(10)

    public static void main(String[] args) {
        long t1 = System.nanoTime();
        ParametersZENOConfig params = new ParametersZENOConfig();

        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.numAtoms = 100;
            params.steps = 100;
            params.walks = 100000;
            params.eigen = true;
            params.pForward = 1.0/6.0;
            params.pBack = 1.0/6.0;
        }

        System.out.println("number of atoms: "+params.numAtoms);
        System.out.println("steps: "+params.steps);
        System.out.println("walks: "+params.walks);
        System.out.println("eigenvalues for q: "+(params.eigen ? "on" : "off"));
        Simulation sim = new Simulation(Space3D.getInstance());
//        sim.setRandom(new RandomMersenneTwister(5));
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
        IntegratorMC integrator = new IntegratorMC(new PotentialComputeAggregate(),sim.getRandom(),1,sim.box());
        MCMoveIdealChain moveIdealChain = new MCMoveIdealChain(sim.box(), sim.getRandom(), params.pBack, params.pForward);
        integrator.getMoveManager().addMCMove(moveIdealChain);

        MeterZENO meterZENO = new MeterZENO(sim.box(), sim.getRandom(), new double[]{1}, false, WalkerExterior.Shape.CUBE);
        meterZENO.setEigenvaluesForPade(params.eigen);
        meterZENO.setNumWalks(params.walks);
        AccumulatorAverageFixed accZENO = new AccumulatorAverageFixed(1);
        DataPumpListener pumpZENO = new DataPumpListener(meterZENO, accZENO);
        integrator.getEventManager().addListener(pumpZENO);

        MeterZENO meterZENOeq21 = new MeterZENO(sim.box(), sim.getRandom(), new double[]{1}, true, WalkerExterior.Shape.CUBE);
        meterZENOeq21.setEigenvaluesForPade(true);
        meterZENOeq21.setNumWalks(params.walks);
        meterZENOeq21.setBoundingSphere(new Vector3D(), 25);
        integrator.getEventManager().addListener(new IntegratorListenerAction(meterZENOeq21));

        if (false) {
            ActivityIntegrate ai =new ActivityIntegrate(integrator);
            sim.getController().addActivity(ai);
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "ZENO");
            simGraphic.makeAndDisplayFrame();
            return;
        }
        ActivityIntegrate ai =new ActivityIntegrate(integrator,params.steps);
        sim.getController().runActivityBlocking(ai);

        IData dataAvg = accZENO.getData(accZENO.AVERAGE);
        IData dataErr = accZENO.getData(accZENO.ERROR);
        System.out.println("Hydrodynamic radius: " + dataAvg.getValue(0)+" "+dataErr.getValue(0));
        System.out.println("Practical conductivity: " + dataAvg.getValue(2)+" "+dataErr.getValue(2));
        System.out.println("Shape factor (q): " + dataAvg.getValue(3)+" "+dataErr.getValue(3));
        System.out.println("Practical viscosity: " + dataAvg.getValue(5)+" "+dataErr.getValue(5));
        System.out.println("Practical viscosity just W: " + dataAvg.getValue(6)+" "+dataErr.getValue(6));

        System.out.println();
        System.out.println("Eq 21 Intrinsic Viscosity: "+meterZENOeq21.getIntrinsicViscosityEq21());
        long t2 = System.nanoTime();
        System.out.println("time: " + (double)(t2 - t1) / (double)1.0E9F);
    }

    public static class ParametersZENOConfig extends ParameterBase {
        public int numAtoms;
        public long steps = 1000;
        public long walks = 1000;
        public boolean eigen;
        public double pBack, pForward;
        public boolean doEq21;
    }
}

