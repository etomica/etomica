/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.freeenergy.npath;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.*;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.chem.elements.Iron;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.graphics.ColorScheme;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.meam.P2EAM;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.PotentialCalculationEnergySum;
import etomica.potential.PotentialCalculationForceSum;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.ElectronVolt;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
 
public class SimFeFCC extends Simulation {
    
    public final PotentialMasterList potentialMaster;
    public final ActivityIntegrate ai;
    public IntegratorVelocityVerlet integrator;
    public SpeciesSpheresMono species;
    public IBox box;
    public P2EAM potential;
    public P1ImageHarmonic p1ImageHarmonic;

    public SimFeFCC(int numAtoms, double temperature, double density, double w, int offsetDim) {
        super(Space3D.getInstance());
        species = new SpeciesSpheresMono(space, Iron.INSTANCE);
        species.setIsDynamic(true);
        addSpecies(species);
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
        IVectorMutable l = space.makeVector();
        l.E(10);
        for (int i=0; i<=offsetDim; i++) {
            l.setX(i,20);
        }
        box.getBoundary().setBoxSize(l);

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();
    
        double n = 8.7932;
        double m = 8.14475;
        double eps = ElectronVolt.UNIT.toSim(0.0220225);
        double a = 3.48501;
        double C = 28.8474;
        double rc = 6;
        potential = new P2EAM(space, n, m, eps, a, C, rc, rc);
    
        potentialMaster = new PotentialMasterList(this, 1.2*rc, space);
        potentialMaster.setCellRange(2);
        double sigma = 1.0;
        integrator = new IntegratorVelocityVerlet(potentialMaster, random, 0.001, temperature, space);
        integrator.setTemperature(temperature);
        integrator.getEventManager().addListener(new IIntegratorListenerMD() {
            PotentialCalculationEnergySum pcEnergy = new PotentialCalculationEnergySum();
            IteratorDirective id = new IteratorDirective();
            @Override
            public void integratorForcePrecomputed(IIntegratorEvent e) {
                potential.disableEnergy();
                potential.reset();
                potentialMaster.calculate(box, id, pcEnergy);
                potential.enableEnergy();
            }
    
            @Override
            public void integratorForceComputed(IIntegratorEvent e) {}
    
            @Override
            public void integratorInitialized(IIntegratorEvent e) {}
    
            @Override
            public void integratorStepStarted(IIntegratorEvent e) {}
    
            @Override
            public void integratorStepFinished(IIntegratorEvent e) {}
        });
        integrator.setForceSum(new PotentialCalculationForceSum());

        ai = new ActivityIntegrate(integrator);
        getController().addAction(ai);

        IAtomType leafType = species.getLeafType();

        potentialMaster.addPotential(potential,new IAtomType[]{leafType,leafType});

        IVectorMutable offset = space.makeVector();
        offset.setX(offsetDim, box.getBoundary().getBoxSize().getX(offsetDim)*0.5);
        p1ImageHarmonic = new P1ImageHarmonic(space, offset, w);
        //potentialMaster.addPotential(p1ImageHarmonic, new IAtomType[]{leafType});

        integrator.setBox(box);
		
        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box);
    
        potentialMaster.getNeighborManager(box).reset();
    }
    
    public static void main(String[] args) {

        LjMC3DParams params = new LjMC3DParams();
        ParseArgs.doParseArgs(params, args);
        if (args.length==0) {
            params.graphics = true;
            params.numAtoms = 2048*2;
            params.steps = 1000000;
            params.density = 4.0/6.5;
            params.T = 1000;
            params.w = 0;
        }

        final int numAtoms = params.numAtoms;
        final double temperatureK = params.T;
        final double temperature = Kelvin.UNIT.toSim(temperatureK);
        final double density = params.density;
        long steps = params.steps;
        boolean graphics = params.graphics;
        double w = params.w;
        int offsetDim = params.offsetDim;

        if (!graphics) {
            System.out.println("Running LJ MC with N="+numAtoms+" at rho="+density+" T="+temperatureK);
            System.out.println(steps+" steps");
            System.out.println("w: "+w);
        }

        double L = Math.pow(numAtoms/density, 1.0/3.0);
        final SimFeFCC sim = new SimFeFCC(numAtoms, temperature, density, w, offsetDim);

        DataSourceEnergies dsEnergies = new DataSourceEnergies(sim.potentialMaster);
        dsEnergies.setPotentialCalculation(new DataSourceEnergies.PotentialCalculationEnergiesEAM(sim.potential));
        dsEnergies.setBox(sim.box);
        IData u = dsEnergies.getData();
        if (!graphics) System.out.println("Fe lattice energy: "+u.getValue(1)/numAtoms);

        if (graphics) {
            final String APP_NAME = "SimLJ";
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 3, sim.getSpace(), sim.getController());
            ColorScheme colorScheme = new ColorScheme() {
                @Override
                public Color getAtomColor(IAtom a) {
                    return a.getLeafIndex() < numAtoms/2 ? Color.RED : Color.BLUE;
                }
            };
            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);
            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

            simGraphic.makeAndDisplayFrame(APP_NAME);

            return;
        }

        long eqSteps = steps/10;
        sim.ai.setMaxSteps(eqSteps);
        sim.getController().actionPerformed();
        sim.getController().reset();
        sim.ai.setMaxSteps(steps);

        System.out.println("equilibration finished ("+eqSteps+" steps)");

        long t1 = System.currentTimeMillis();

        long blockSize = steps/numAtoms/100;
        if (blockSize==0) blockSize = 1;

        AccumulatorAverageFixed accEnergies = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pumpEnergies = new DataPumpListener(dsEnergies, accEnergies, numAtoms);
        sim.integrator.getEventManager().addListener(pumpEnergies);

        sim.getController().actionPerformed();

        IData avgEnergies = accEnergies.getData(accEnergies.AVERAGE);
        IData errEnergies = accEnergies.getData(accEnergies.ERROR);
        IData corEnergies = accEnergies.getData(accEnergies.BLOCK_CORRELATION);

        System.out.println("spring energy: "+avgEnergies.getValue(0)/numAtoms+"   error: "+errEnergies.getValue(0)/numAtoms+"  cor: "+corEnergies.getValue(0));
        System.out.println("LJ energy: "+avgEnergies.getValue(1)/numAtoms+"   error: "+errEnergies.getValue(1)/numAtoms+"  cor: "+corEnergies.getValue(1));

        long t2 = System.currentTimeMillis();
        System.out.println("time: "+(t2-t1)/1000.0+" seconds");
    }

    public static class LjMC3DParams extends ParameterBase {
        public int numAtoms = 500;
        public double T = 2.0;
        public double density = 0.3;
        public long steps = 100000;
        public double rc = 2.5;
        public boolean graphics = false;
        public double w = 1;
        public int offsetDim = 0;
    }

}
