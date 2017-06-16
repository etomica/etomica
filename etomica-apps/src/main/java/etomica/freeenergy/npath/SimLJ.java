/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.freeenergy.npath;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.CompoundUnit;
import etomica.units.Null;
import etomica.units.SimpleUnit;
import etomica.units.Unit;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
 
public class SimLJ extends Simulation {
    
    public final PotentialMasterCell potentialMasterCell;
    public final ActivityIntegrate ai;
    public IntegratorMC integrator;
    public SpeciesSpheresMono species;
    public Box box;
    public P2LennardJones potential;
    public MCMoveAtomNPath mcMoveAtom;
    public MCMoveAtomCoupled mcMoveAtomCoupled;
    public MCMoveAtomSwap mcMoveSwap;
    public P1ImageHarmonic p1ImageHarmonic;

    public SimLJ(int numAtoms, double temperature, double density, double rc, double w, int offsetDim) {
        super(Space3D.getInstance());
        species = new SpeciesSpheresMono(this, space);
        addSpecies(species);
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
        Vector l = space.makeVector();
        l.E(10);
        for (int i=0; i<=offsetDim; i++) {
            l.setX(i,20);
        }
        box.getBoundary().setBoxSize(l);

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();
        
        potentialMasterCell = new PotentialMasterCell(this, rc, space);
        potentialMasterCell.setCellRange(2);
        double sigma = 1.0;
        integrator = new IntegratorMC(this, potentialMasterCell);
        integrator.setTemperature(temperature);

        ai = new ActivityIntegrate(integrator);
        getController().addAction(ai);

        potential = new P2LennardJones(space, sigma, 1);
        AtomType leafType = species.getLeafType();
        P2SoftSphericalTruncated potentialTruncated = new P2SoftSphericalTruncated(space, potential, rc);

        potentialMasterCell.addPotential(potentialTruncated, new AtomType[]{leafType, leafType});

        Vector offset = space.makeVector();
        offset.setX(offsetDim, box.getBoundary().getBoxSize().getX(offsetDim)*0.5);
        p1ImageHarmonic = new P1ImageHarmonic(space, offset, w, true);
        potentialMasterCell.addPotential(p1ImageHarmonic, new AtomType[]{leafType});

        mcMoveAtom = new MCMoveAtomNPath(random, potentialMasterCell, space, p1ImageHarmonic);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        mcMoveAtomCoupled = new MCMoveAtomCoupled(random, potentialMasterCell, space);
        integrator.getMoveManager().addMCMove(mcMoveAtomCoupled);
//        ((MCMoveStepTracker)mcMoveAtom.getTracker()).setNoisyAdjustment(true);
//        ((MCMoveStepTracker)mcMoveAtomCoupled.getTracker()).setNoisyAdjustment(true);

        Vector boxLength = box.getBoundary().getBoxSize();
        double lMin = boxLength.getX(0);
        if (boxLength.getX(1) < lMin) lMin = boxLength.getX(1);
        if (boxLength.getX(2) < lMin) lMin = boxLength.getX(2);
        double ww = w / lMin;
        double swapDistance = 1.5*Math.sqrt(1.5*temperature/ww);
        if (swapDistance > lMin/4) swapDistance = lMin/4;
        if (swapDistance > rc) swapDistance = rc;
        if (swapDistance < 1) swapDistance = 1;
        mcMoveSwap = new MCMoveAtomSwap(random, potentialMasterCell, space, p1ImageHarmonic);
        mcMoveSwap.setNbrDistance(swapDistance);
        integrator.getMoveManager().addMCMove(mcMoveSwap);
        integrator.getMoveManager().setFrequency(mcMoveSwap,5);

        integrator.setBox(box);

        integrator.getMoveEventManager().addListener(potentialMasterCell.getNbrCellManager(box).makeMCMoveListener());
		
        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box);
        potentialMasterCell.getNbrCellManager(box).assignCellAll();
        p1ImageHarmonic.findNOffset(box);
    }
    
    public static void main(String[] args) {

        LjMC3DParams params = new LjMC3DParams();
        ParseArgs.doParseArgs(params, args);
        if (args.length==0) {
            params.graphics = true;
            params.numAtoms = 1000;
            params.steps = 1000000;
            params.density = 1.3;
            params.T = 2;
            params.w = 500;
        }

        final int numAtoms = params.numAtoms;
        final double temperature = params.T;
        final double density = params.density;
        long steps = params.steps;
        double rc = params.rc*Math.pow(params.density, -1.0/3.0);
        boolean graphics = params.graphics;
        double w = params.w;
        int offsetDim = params.offsetDim;

        if (!graphics) {
            System.out.println("Running LJ MC with N="+numAtoms+" at rho="+density+" T="+temperature);
            System.out.println(steps+" steps");
            System.out.println("cutoff: "+rc);
            System.out.println("w: "+w);
        }

        double L = Math.pow(numAtoms/density, 1.0/3.0);
        final SimLJ sim = new SimLJ(numAtoms, temperature, density, rc, w, offsetDim);

        DataSourceEnergies dsEnergies = new DataSourceEnergies(sim.potentialMasterCell);
        dsEnergies.setBox(sim.box);
        IData u = dsEnergies.getData();
        if (!graphics) System.out.println("LJ lattice energy: "+u.getValue(1)/numAtoms);

        if (graphics) {
            final String APP_NAME = "SimLJ";
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 3, sim.getSpace(), sim.getController());
            ColorScheme colorScheme = new ColorScheme() {
                @Override
                public Color getAtomColor(IAtom a) {
                    int idx0 = a.getLeafIndex();
                    int idx1 = sim.p1ImageHarmonic.getPartner(idx0);
                    return idx0 < idx1 ? Color.RED : Color.BLUE;
                }
            };
            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);
            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
            
            DataSourceCountSteps tSource = new DataSourceCountSteps(sim.integrator);
            AccumulatorHistory energyHist = new AccumulatorHistory(new HistoryCollapsingAverage());
            energyHist.setTimeDataSource(tSource);
            AccumulatorHistory springHist = new AccumulatorHistory(new HistoryCollapsingAverage());
            springHist.setTimeDataSource(tSource);
            DataSplitter splitter = new DataSplitter();
            DataPumpListener energyPump = new DataPumpListener(dsEnergies, splitter, numAtoms);
            sim.integrator.getEventManager().addListener(energyPump);
            splitter.setDataSink(0, springHist);
            splitter.setDataSink(1, energyHist);
            DisplayPlot energyPlot = new DisplayPlot();
            energyPlot.setLabel("LJ");
            energyPlot.setUnit(new CompoundUnit(new Unit[]{new SimpleUnit(Null.DIMENSION,1.0/numAtoms,"why do you want a name.  just use it.","per atom", false)},new double[]{-1}));
//            energyPlot.setUnit(new CompoundUnit(new Unit[]{ElectronVolt.UNIT, new SimpleUnit(Null.DIMENSION,1.0/numAtoms,"why do you want a name.  just use it.","per atom", false)},new double[]{1,-1}));
//            energyPlot.setUnit(ElectronVolt.UNIT);
            energyHist.addDataSink(energyPlot.getDataSet().makeDataSink());
            simGraphic.add(energyPlot);
            DisplayPlot springPlot = new DisplayPlot();
            springPlot.setLabel("spring");
            springPlot.setUnit(new CompoundUnit(new Unit[]{new SimpleUnit(Null.DIMENSION,1.0/numAtoms,"why do you want a name.  just use it.","per atom", false)},new double[]{-1}));
            springHist.addDataSink(springPlot.getDataSet().makeDataSink());
            simGraphic.add(springPlot);
    
            simGraphic.makeAndDisplayFrame(APP_NAME);

            return;
        }

        long eqSteps = steps/10;
        sim.ai.setMaxSteps(eqSteps);
        sim.getController().actionPerformed();
        sim.getController().reset();
        sim.integrator.resetStepCount();
        sim.integrator.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(steps);

        System.out.println("equilibration finished ("+eqSteps+" steps)");

        long t1 = System.currentTimeMillis();

        long blockSize = steps/numAtoms/100;
        if (blockSize==0) blockSize = 1;

        AccumulatorAverageCovariance accEnergies = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener pumpEnergies = new DataPumpListener(dsEnergies, accEnergies, numAtoms);
        sim.integrator.getEventManager().addListener(pumpEnergies);

        sim.getController().actionPerformed();

        IData avgEnergies = accEnergies.getData(AccumulatorAverage.AVERAGE);
        IData errEnergies = accEnergies.getData(AccumulatorAverage.ERROR);
        IData corEnergies = accEnergies.getData(AccumulatorAverage.BLOCK_CORRELATION);
        IData covEnergies = accEnergies.getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE);

        System.out.println("swap acceptance: "+sim.mcMoveSwap.getTracker().acceptanceProbability());
        System.out.println("simple move step size: " + sim.mcMoveAtom.getStepSize());
        System.out.println("coupled move step size: " + sim.mcMoveAtomCoupled.getStepSize());

        System.out.println("spring energy: "+avgEnergies.getValue(0)/numAtoms+"   error: "+errEnergies.getValue(0)/numAtoms+"  cor: "+corEnergies.getValue(0));
        System.out.println("LJ energy: "+avgEnergies.getValue(1)/numAtoms+"   error: "+errEnergies.getValue(1)/numAtoms+"  cor: "+corEnergies.getValue(1));
        System.out.println("du/dw: "+avgEnergies.getValue(2)/numAtoms+"   error: "+errEnergies.getValue(2)/numAtoms+"  cor: "+corEnergies.getValue(2));
        double var0 = covEnergies.getValue(0*3+0);
        double var1 = covEnergies.getValue(1*3+1);
        double var2 = covEnergies.getValue(2*3+2);
        double cor01 = 0;
        if (var0*var1>0) cor01 = covEnergies.getValue(0*3+1)/Math.sqrt(covEnergies.getValue(0*3+0)*covEnergies.getValue(1*3+1));
        double cor02 = 0;
        if (var0*var2>0) cor02 = covEnergies.getValue(0*3+2)/Math.sqrt(covEnergies.getValue(0*3+0)*covEnergies.getValue(2*3+2));
        System.out.println("spring correlation: 1 "+cor01+" "+cor02);
        double cor12 = covEnergies.getValue(1*3+2)/Math.sqrt(covEnergies.getValue(1*3+1)*covEnergies.getValue(2*3+2));
        System.out.println("Fe correlation: "+cor01+" 1 "+cor12);
        System.out.println("du/dw correlation: "+cor02+" "+cor12+" 1");
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
