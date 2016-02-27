/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import java.awt.Color;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtom;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.IData;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.liquidLJ.DataProcessorReweight;
import etomica.liquidLJ.DataProcessorReweightRatio;
import etomica.liquidLJ.ValueCache;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.Potential2SoftSpherical;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;



public class SimLJHTTISuper extends Simulation {

    public SimLJHTTISuper(Space _space, int numAtoms, double density, double temperature, double rc, boolean ss) {
        super(_space);
        
        potentialMaster = new PotentialMasterList(this, space);
                
        species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        // TARGET
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);

        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(box);
        atomMove = new MCMoveAtomCoupled(potentialMaster, meterPE, getRandom(), space);
        atomMove.setStepSize(0.1);
        atomMove.setStepSizeMax(0.5);
        atomMove.setDoExcludeNonNeighbors(true);
        integrator.getMoveManager().addMCMove(atomMove);
//        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);
        
        double L = Math.pow(4.0/density, 1.0/3.0);
        double nbrDistance = L / Math.sqrt(2);
        int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
        primitive = new PrimitiveCubic(space, n*L);
        
        nCells = new int[]{n,n,n};
        boundary = new BoundaryRectangularPeriodic(space, n * L);
        
        Basis basisFCC = new BasisCubicFcc();
        basis = new BasisBigCell(space, basisFCC, nCells);
    
        box.setBoundary(boundary);

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1,1,1});

        potential = ss ? new P2SoftSphere(space, 1.0, 4.0, 12) : new P2LennardJones(space, 1.0, 1.0);
        potential = new P2SoftSphericalTruncated(space, potential, rc);
        atomMove.setPotential(potential);
        IAtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(potential, new IAtomType[] {sphereType, sphereType });

        /*
         *  1-body Potential to Constraint the atom from moving too far
         *  	away from its lattice-site
         *  
         */

        P1ConstraintNbr p1Constraint = new P1ConstraintNbr(space, nbrDistance, this);
        p1Constraint.initBox(box);
        atomMove.setConstraint(p1Constraint);

        potentialMaster.lrcMaster().setEnabled(false);
    
        integrator.setBox(box);

        int cellRange = 7;
        potentialMaster.setRange(rc);
        potentialMaster.setCellRange(cellRange); // insanely high, this lets us have neighborRange close to dimensions/2
        // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
        potentialMaster.getNeighborManager(box).reset();
        int potentialCells = potentialMaster.getNbrCellManager(box).getLattice().getSize()[0];
        if (potentialCells < cellRange*2+1) {
            throw new RuntimeException("oops ("+potentialCells+" < "+(cellRange*2+1)+")");
        }
        
        activityIntegrate = new ActivityIntegrate(integrator);
        
        getController().addAction(activityIntegrate);
        
        // extend potential range, so that atoms that move outside the truncation range will still interact
        // atoms that move in will not interact since they won't be neighbors
        ((P2SoftSphericalTruncated)potential).setTruncationRadius(0.6*boundary.getBoxSize().getX(0));
    }
    
    public void initialize(long initSteps) {
        // equilibrate off the lattice to avoid anomolous contributions
        activityIntegrate.setMaxSteps(initSteps);
        getController().actionPerformed();
        getController().reset();
        integrator.getMoveManager().setEquilibrating(false);
    }
    
    /**
     * @param args filename containing simulation parameters
     * @see SimLJHTTI.SimOverlapParam
     */
    public static void main(String[] args) {
        //set up simulation parameters
        SimOverlapParam params = new SimOverlapParam();
        if (args.length == 0) {
            params.numAtoms = 864;
            params.numSteps = 10000000;
            params.temperature = 1;
            params.density = 1;
            params.rc = 3*Math.pow(params.density, -1.0/3.0);
            params.bpharm = new double[]{9.557960628737128,9.567450981563097};
        }
        else {
            ParseArgs.doParseArgs(params, args);
        }
        boolean ss = params.ss;
        double density = params.density;
        long numSteps = params.numSteps;
        final int numAtoms = params.numAtoms;
        double temperature = params.temperature;
        double rc = params.rc;
        double[] bpharm = params.bpharm;
        
        System.out.println("Running "+(ss?"soft-sphere":"Lennard-Jones")+" simulation");
        System.out.println(numAtoms+" atoms at density "+density+" and temperature "+temperature);
        System.out.println(numSteps+" steps");

        //instantiate simulation
        final SimLJHTTISuper sim = new SimLJHTTISuper(Space.getInstance(3), numAtoms, density, temperature, rc, ss);
        if (false) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, sim.space, sim.getController());
            simGraphic.setPaintInterval(sim.box, 1000);
            ColorScheme colorScheme = new ColorScheme() {
                public Color getAtomColor(IAtom a) {
                    if (allColors==null) {
                        allColors = new Color[768];
                        for (int i=0; i<256; i++) {
                            allColors[i] = new Color(255-i,i,0);
                        }
                        for (int i=0; i<256; i++) {
                            allColors[i+256] = new Color(0,255-i,i);
                        }
                        for (int i=0; i<256; i++) {
                            allColors[i+512] = new Color(i,0,255-i);
                        }
                    }
                    return allColors[(2*a.getLeafIndex()) % 768];
                }
                protected Color[] allColors;
            };
            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);
            
            DisplayTextBox timer = new DisplayTextBox();
            DataSourceCountSteps counter = new DataSourceCountSteps(sim.integrator);
            DataPumpListener counterPump = new DataPumpListener(counter, timer, 100);
            sim.integrator.getEventManager().addListener(counterPump);
            simGraphic.getPanel().controlPanel.add(timer.graphic());
            
            simGraphic.makeAndDisplayFrame("SS");
            return;
        }

        //start simulation

        double L1 = Math.pow(numAtoms, 1.0/3.0);
        double rc1 = rc*Math.pow(density, 1.0/3.0);
        double rcMax1 = 0.494*L1;
        double fac = 1;
        int nCutoffs = 1 + (int)((rcMax1-rc1)/fac);
        if (nCutoffs > bpharm.length) {
            throw new RuntimeException("need more beta P harmonic");
        }
        
        double[] cutoffs = new double[nCutoffs];
        cutoffs[0] = rc1;
        for (int i=1; i<nCutoffs; i++) {
            cutoffs[i] = cutoffs[i-1] + fac;
        }
        for (int i=0; i<nCutoffs; i++) {
            cutoffs[i] *= Math.pow(density, -1.0/3.0);
        }
        System.out.print("cutoffs: ");
        for (int i=0; i<nCutoffs; i++) {
            System.out.print(" "+cutoffs[i]);
        }
        System.out.println();

        PotentialMasterList potentialMasterData;
        potentialMasterData = new PotentialMasterList(sim, cutoffs[nCutoffs-1], sim.getSpace());
        Potential2SoftSpherical potential = ss ? new P2SoftSphere(sim.getSpace(), 1.0, 4.0, 12) : new P2LennardJones(sim.getSpace(), 1.0, 1.0);
        potential = new P2SoftSphericalTruncated(sim.getSpace(), potential, cutoffs[nCutoffs-1]-0.01);
        IAtomType sphereType = sim.species.getLeafType();
        potentialMasterData.addPotential(potential, new IAtomType[] {sphereType, sphereType });
        potentialMasterData.lrcMaster().setEnabled(false);

        int cellRange = 7;
        potentialMasterData.setCellRange(cellRange); // insanely high, this lets us have neighborRange close to dimensions/2
        // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
        potentialMasterData.getNeighborManager(sim.box).reset();
        int potentialCells = potentialMasterData.getNbrCellManager(sim.box).getLattice().getSize()[0];
        if (potentialCells < cellRange*2+1) {
            throw new RuntimeException("oops ("+potentialCells+" < "+(cellRange*2+1)+")");
        }
        
        // extend potential range, so that atoms that move outside the truncation range will still interact
        // atoms that move in will not interact since they won't be neighbors
        ((P2SoftSphericalTruncated)potential).setTruncationRadius(0.6*sim.box.getBoundary().getBoxSize().getX(0));
        
        // meter needs lattice energy, so make it now
        MeterSolidDACut meterSolid = new MeterSolidDACut(sim.getSpace(), potentialMasterData, sim.coordinateDefinition, cutoffs);
        meterSolid.setTemperature(temperature);
        meterSolid.setBPRes(bpharm);
        IData d = meterSolid.getData();
        
        MeterPotentialEnergy meterEnergyShort = new MeterPotentialEnergy(sim.potentialMaster);
        meterEnergyShort.setBox(sim.box);
        final double[] uFacCut = new double[cutoffs.length];
        for (int i=0; i<uFacCut.length; i++) {
            uFacCut[i] = d.getValue(5*i)*numAtoms - meterEnergyShort.getDataAsScalar();
        }
        System.out.print("bPharm ");
        for (int i=0; i<nCutoffs; i++) {
            System.out.print(" "+bpharm[i]);
        }
        System.out.println();
        
        if (args.length == 0) {
            // quick initialization
            sim.initialize(numSteps/10);
        }
        else {
            sim.initialize(numSteps/20 + 50*numAtoms + numAtoms*numAtoms*3);
        }

        int numBlocks = 100;
        int interval = numAtoms;
        long blockSize = numSteps/(numBlocks*interval);
        if (blockSize == 0) blockSize = 1;
        System.out.println("block size "+blockSize+" interval "+interval);

        final ValueCache energyFastCache = new ValueCache(meterEnergyShort, sim.integrator);

        DataProcessorReweight puReweight = new DataProcessorReweight(temperature, energyFastCache, uFacCut, sim.box, cutoffs.length);
        long bs = numSteps/(numAtoms*100);
        DataPumpListener pumpPU = new DataPumpListener(meterSolid, puReweight, numAtoms);
        sim.integrator.getEventManager().addListener(pumpPU);
        final AccumulatorAverageCovariance avgSolid = new AccumulatorAverageCovariance(bs == 0 ? 1 : bs);
        puReweight.setDataSink(avgSolid);

        DataProcessorReweightRatio puReweightRatio = new DataProcessorReweightRatio(cutoffs.length);
        avgSolid.setBlockDataSink(puReweightRatio);
        AccumulatorAverageCovariance accPUBlocks = new AccumulatorAverageCovariance(1, true);
        puReweightRatio.setDataSink(accPUBlocks);

        
        final long startTime = System.currentTimeMillis();
       
        sim.activityIntegrate.setMaxSteps(numSteps);

        sim.getController().actionPerformed();
        long endTime = System.currentTimeMillis();
        System.out.println();

        IData avgRawData = avgSolid.getData(avgSolid.AVERAGE);
        IData errRawData = avgSolid.getData(avgSolid.ERROR);
        IData corRawData = avgSolid.getData(avgSolid.BLOCK_CORRELATION);

        int j = 6;
        for (int i=1; i<cutoffs.length; i++) {
            double avgW = avgRawData.getValue(j+5);
            double errW = errRawData.getValue(j+5);
            double corW = corRawData.getValue(j+5);
            System.out.println(String.format("rc: %2d dbA:   % 21.15e  %10.4e  % 5.3f  %6.4f", i, -Math.log(avgW)/numAtoms, errW/avgW/numAtoms, corW, errW/avgW));
            j += 6;
        }
        System.out.println("\n");

        
        IData avgData = accPUBlocks.getData(avgSolid.AVERAGE);
        IData errData = accPUBlocks.getData(avgSolid.ERROR);
        IData corData = accPUBlocks.getData(avgSolid.BLOCK_CORRELATION);
        IData covData = accPUBlocks.getData(avgSolid.BLOCK_COVARIANCE);

        int n = avgData.getLength();
        j = 0;
        for  (int i=0; i<cutoffs.length; i++) {
            double avgU = avgData.getValue(j+0);
            double errU = errData.getValue(j+0);
            double corU = corData.getValue(j+0);
            double avgP = avgData.getValue(j+1);
            double errP = errData.getValue(j+1);
            double corP = corData.getValue(j+1);
            double avgBUc = avgData.getValue(j+2);
            double errBUc = errData.getValue(j+2);
            double corBUc = corData.getValue(j+2);
            double avgZc = avgData.getValue(j+3);
            double errZc = errData.getValue(j+3);
            double corZc = corData.getValue(j+3);
            // this is dbAc/drho at constant Y (for LJ)
            double avgDADv2 = avgData.getValue(j+4);
            double errDADv2 = errData.getValue(j+4);
            double corDADv2 = corData.getValue(j+4);
    
            double DADACor = covData.getValue(2*n+4)/Math.sqrt(covData.getValue(2*n+2)*covData.getValue(4*n+4));
            double ZcUcCor = covData.getValue(3*n+4)/Math.sqrt(covData.getValue(3*n+3)*covData.getValue(4*n+4));
            double facDADY = 4*density*density*density*density/temperature;
    
            System.out.print(String.format("rc: %2d DADY:  % 21.15e  %10.4e  % 5.3f\n", i, -facDADY*avgBUc, facDADY*errBUc, corBUc));
            if (!ss) System.out.print(String.format("rc: %2d DADv2: % 21.15e  %10.4e  % 5.3f  % 8.6f\n", i, avgDADv2, errDADv2, corDADv2, DADACor));
            System.out.print(String.format("rc: %2d Zc:    % 21.15e  %10.4e  % 5.3f\n", i, avgZc, errZc, corZc));
            System.out.print(String.format("rc: %2d bUc:   % 21.15e  %10.4e  % 5.3f  % 8.6f\n", i, avgBUc, errBUc, corBUc, ZcUcCor));
    
            double PUCor = covData.getValue(1*n+0)/Math.sqrt(covData.getValue(1*n+1)*covData.getValue(0*n+0));
            System.out.print(String.format("rc: %2d Uraw:  % 21.15e  %10.4e  % 5.3f\n", i, avgU, errU, corU));
            System.out.print(String.format("rc: %2d Praw:  % 21.15e  %10.4e  % 5.3f  % 8.6f\n", i, avgP, errP, corP, PUCor));
            System.out.println();
            j+=5;
        }

        System.out.println("time: " + (endTime - startTime)/1000.0);
    }

    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public IBox box;
    public Boundary boundary;
    public int[] nCells;
    public Basis basis;
    public Primitive primitive;
    public MCMoveAtomCoupled atomMove;
    public PotentialMasterList potentialMaster;
    public final CoordinateDefinitionLeaf coordinateDefinition;
    public Potential2SoftSpherical potential;
    public SpeciesSpheresMono species;
    
    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numAtoms = 256;
        public double density = 1.28;
        public long numSteps = 1000000;
        public double temperature = 0.1;
        public double rc = 2.7;
        public double[] bpharm = new double[0];
        public boolean ss = false;
    }
}
