/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.liquidLJ.*;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.list.BoxAgentSourceCellManagerList;
import etomica.nbr.list.NeighborListManagerSlanty;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Degree;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

import java.awt.*;
import java.util.Arrays;



public class SimLJHTTISuperHCP extends Simulation {

    public final CoordinateDefinitionLeaf coordinateDefinition;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public Box box;
    public BoundaryDeformablePeriodic boundary;
    public int[] nCells;
    public Basis basis;
    public Primitive primitive;
    public MCMoveAtomCoupled atomMove;
    public PotentialMasterList potentialMaster;
    public Potential2SoftSpherical potential;
    public SpeciesSpheresMono species;
    public SimLJHTTISuperHCP(Space _space, int numAtoms, double density, double coa, double temperature, double rc, boolean ss, int[] seeds) {
        super(_space);
        if (seeds != null) {
            setRandom(new RandomMersenneTwister(seeds));
        }

        species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        // TARGET
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);

        int n = (int)Math.round(Math.pow(numAtoms/8, 1.0/3.0));
        if (8*n*n*n != numAtoms) {
            throw new RuntimeException("Not compatible with HCP");
        }
        // V = nc
        // v = 2/density = (2^.5/rho)sqrt(8/3)*f
        //               = 4f/(rho sqrt(3))
        // f = sqrt(3)/2
        // v = sqrt(3)/2 a^3 coa
        // a = (2 v / (sqrt(3) coa))^(1/3)
        //   = (4 / (sqrt(3) rho coa))^(1/3)
        double a = Math.pow(4/(Math.sqrt(3)*density*coa), 1.0/3.0);
        double c = coa*a;  // sqrt(8/3)
        Vector[] boxDim = new Vector[3];
        boxDim[0] = space.makeVector(new double[]{2*n*a, 0, 0});
        boxDim[1] = space.makeVector(new double[]{-2*n*a*Math.cos(Degree.UNIT.toSim(60)), 2*n*a*Math.sin(Degree.UNIT.toSim(60)), 0});
        boxDim[2] = space.makeVector(new double[]{0, 0, n*c});

        primitive = new PrimitiveHexagonal(space, a, c);
        nCells = new int[]{2*n,2*n,n};
        boundary = new BoundaryDeformableLattice(primitive, nCells);
        boundary.setTruncationRadius(rc);
        basis = new BasisHcp();

        box.setBoundary(boundary);

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(nCells);

        if (rc > 0.494*c*n) {
            throw new RuntimeException("cutoff too big");
        }

        BoxAgentSourceCellManagerList boxAgentSource = new BoxAgentSourceCellManagerList(this, null, _space);
        BoxAgentManager<NeighborCellManager> boxAgentManager = new BoxAgentManager<NeighborCellManager>(boxAgentSource, this);
        potentialMaster = new PotentialMasterList(this, rc, boxAgentSource, boxAgentManager, new NeighborListManagerSlanty.NeighborListSlantyAgentSource(rc, space), space);

        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(box);
        atomMove = new MCMoveAtomCoupled(potentialMaster, meterPE, getRandom(), space);
        atomMove.setStepSize(0.1);
        atomMove.setStepSizeMax(0.5);
        atomMove.setDoExcludeNonNeighbors(true);
        integrator.getMoveManager().addMCMove(atomMove);
//        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);

        potential = ss ? new P2SoftSphere(space, 1.0, 4.0, 12) : new P2LennardJones(space, 1.0, 1.0);
        potential = new P2SoftSphericalTruncated(space, potential, rc);
        atomMove.setPotential(potential);
        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(potential, new AtomType[]{sphereType, sphereType});

        /*
         *  1-body Potential to Constraint the atom from moving too far
         *  	away from its lattice-site
         *
         */

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

    /**
     * @param args filename containing simulation parameters
     * @see SimLJHTTISuperHCP.SimOverlapParam
     */
    public static void main(String[] args) {
        //set up simulation parameters
        SimOverlapParam params = new SimOverlapParam();
        if (args.length == 0) {
            params.numAtoms = 4*4*4*8;
            params.numSteps = 10000000;
            params.temperature = 0.7716049382716044;
            params.density = 1.178511301977579;
            params.rcMax1 = 10;
            params.rcMax0 = 11;
            params.rc = 3;
            params.rc0 = 3;
            params.bpharm = new double[]{9.97423151884132,9.97721709086801,9.979113882088319,9.980514061272098,9.980920892526447,9.981107369390061,9.981164701696656};
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
        double rc0 = params.rc0;
        double rcMax0 = params.rcMax0;
        double rcMax1 = params.rcMax1;
        if (rcMax1 > rcMax0) rcMax1 = rcMax0;
        double[] bpharm = params.bpharm;
        double[] bpharmLJ = params.bpharmLJ;
        int[] seeds = params.randomSeeds;
        double alpha = params.alpha;

        System.out.println("Running "+(ss?"soft-sphere":"Lennard-Jones")+" simulation");
        System.out.println(numAtoms+" atoms at density "+density+" and temperature "+temperature);
        System.out.println(numSteps+" steps");

        //instantiate simulation
        final SimLJHTTISuperHCP sim = new SimLJHTTISuperHCP(Space.getInstance(3), numAtoms, density, alpha*Math.sqrt(8.0/3.0), temperature, rc*Math.pow(density, -1.0/3.0), ss, seeds);
        if (seeds == null) {
            seeds = ((RandomMersenneTwister)sim.getRandom()).getSeedArray();
        }
        System.out.println("Random seeds: "+Arrays.toString(seeds));
        if (false) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.setPaintInterval(sim.box, 1000);
            ColorScheme colorScheme = new ColorScheme() {
                protected Color[] allColors;

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


        Vector e2 = sim.box.getBoundary().getEdgeVector(2);
        double L = e2.getX(2)*Math.pow(density, 1.0/3.0);
        if (rcMax1 > 0.494*L) rcMax1 = 0.494*L;
        double delta = 0.5;
        int nCutoffs = 1;
        double c = rc0;
        for (nCutoffs=1; c<=rcMax1*1.0001; nCutoffs++) {
            c += delta;
            if (nCutoffs%2==0) delta += 0.5;
        }
        nCutoffs--;
        if (nCutoffs > bpharm.length) {
            throw new RuntimeException("need more beta P harmonic");
        }
        delta = 0.5;
        double[] cutoffs = new double[nCutoffs];
        cutoffs[0] = rc0;
        for (int i=1; i<nCutoffs; i++) {
            cutoffs[i] = cutoffs[i-1] + delta;
            if (i%2==0) delta += 0.5;
        }
        for (int i=0; i<nCutoffs; i++) {
            cutoffs[i] *= Math.pow(density, -1.0/3.0);
        }

        PotentialMasterList potentialMasterData;
        Potential2SoftSpherical potential = ss ? new P2SoftSphere(sim.getSpace(), 1.0, 4.0, 12) : new P2LennardJones(sim.getSpace(), 1.0, 1.0);
        {

            BoxAgentSourceCellManagerList boxAgentSource = new BoxAgentSourceCellManagerList(sim, null, sim.space);
            BoxAgentManager<NeighborCellManager> boxAgentManager = new BoxAgentManager<NeighborCellManager>(boxAgentSource, sim);
            potentialMasterData = new PotentialMasterList(sim, cutoffs[nCutoffs-1], boxAgentSource, boxAgentManager, new NeighborListManagerSlanty.NeighborListSlantyAgentSource(rc, sim.space), sim.space);

            // |potential| is our local potential used for data collection.
            P2SoftSphericalTruncated potentialT = new P2SoftSphericalTruncated(sim.getSpace(), potential, cutoffs[nCutoffs-1]-0.01);
            AtomType sphereType = sim.species.getLeafType();
            potentialMasterData.addPotential(potentialT, new AtomType[]{sphereType, sphereType});
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
            potentialT.setTruncationRadius(0.6*sim.box.getBoundary().getBoxSize().getX(0));
        }


        PotentialMasterList potentialMasterDataLJ = null;
        P2LennardJones p2LJ = null;
        Potential2SoftSpherical potentialLJ = null;
        if (ss) {
            if (true) throw new RuntimeException("still broken, need slanty nbrs");
            // |potential| is our local potential used for data collection.
            potentialMasterDataLJ = new PotentialMasterList(sim, cutoffs[nCutoffs-1], sim.getSpace());
            p2LJ = new P2LennardJones(sim.getSpace());
            System.out.println("Including 12+6 LJ");

            potentialLJ = new P2SoftSphericalTruncated(sim.getSpace(), p2LJ, cutoffs[nCutoffs-1]-0.01);
            AtomType sphereType = sim.species.getLeafType();
            potentialMasterDataLJ.addPotential(potentialLJ, new AtomType[]{sphereType, sphereType});
            potentialMasterDataLJ.lrcMaster().setEnabled(false);

            int cellRange = 7;
            potentialMasterDataLJ.setCellRange(cellRange); // insanely high, this lets us have neighborRange close to dimensions/2
            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
            potentialMasterDataLJ.getNeighborManager(sim.box).reset();
            int potentialCells = potentialMasterDataLJ.getNbrCellManager(sim.box).getLattice().getSize()[0];
            if (potentialCells < cellRange*2+1) {
                throw new RuntimeException("oops ("+potentialCells+" < "+(cellRange*2+1)+")");
            }

            // extend potential range, so that atoms that move outside the truncation range will still interact
            // atoms that move in will not interact since they won't be neighbors
            ((P2SoftSphericalTruncated)potentialLJ).setTruncationRadius(0.6*sim.box.getBoundary().getBoxSize().getX(0));
        }

        // meter needs lattice energy, so make it now
        MeterSolidDACut meterSolid = new MeterSolidDACut(sim.getSpace(), potentialMasterData, sim.coordinateDefinition, cutoffs);
        meterSolid.setTemperature(temperature);
        meterSolid.setBPRes(bpharm);
        IData d = meterSolid.getData();

        MeterPotentialEnergy meterEnergyShort = new MeterPotentialEnergy(sim.potentialMaster);
        meterEnergyShort.setBox(sim.box);
        final double[] uFacCut = new double[cutoffs.length];
        double uShort = meterEnergyShort.getDataAsScalar();
        for (int i=0; i<uFacCut.length; i++) {
            uFacCut[i] = d.getValue(6*i)*numAtoms - uShort;
        }


        if (ss) {
            if (bpharmLJ.length < cutoffs.length) {
                throw new RuntimeException("I need LJ harmonic pressures for all cutoffs");
            }
            meterSolid.setPotentialMasterDADv2(potentialMasterDataLJ, bpharmLJ);
        }

        double rcMaxLS = 3*0.494*L;
        if (rcMaxLS>rcMax0) rcMaxLS = rcMax0;
        if (rcMax1 >= rcMax0) rcMaxLS=0;

        delta = 0.5;
        int nCutoffsLS = 1;
        c = rc0;
        for (nCutoffsLS=1; c<rcMaxLS*1.0001; nCutoffsLS++) {
            c += delta;
            if (nCutoffsLS%2==0) delta += 0.5;
        }
        nCutoffsLS--;
        if (nCutoffsLS > bpharm.length) {
            throw new RuntimeException("need more beta P harmonic");
        }

        final double[] cutoffsLS = new double[nCutoffsLS];
        PotentialMasterMonatomic potentialMasterLS = new PotentialMasterMonatomic(sim);
        Potential2SoftSphericalLSMultiLatSlanty pLS = null;
        PotentialMasterMonatomic potentialMasterLJLS = null;
        Potential2SoftSphericalLSMultiLat pLJLS = null;
        final double[] uFacCutLS = new double[cutoffsLS.length];
        MeterSolidDACut meterSolidLS = null;
        if (nCutoffsLS>0) {
            cutoffsLS[0] = rc0;
            delta = 0.5;
            for (int i=1; i<cutoffsLS.length; i++) {
                cutoffsLS[i] = cutoffsLS[i-1] + delta;
                if (i%2==0) delta += 0.5;
            }
            for (int i=0; i<nCutoffsLS; i++) {
                cutoffsLS[i] *= Math.pow(density, -1.0/3.0);
            }
            pLS = new Potential2SoftSphericalLSMultiLatSlanty(sim.getSpace(), cutoffsLS, potential, sim.coordinateDefinition, new int[]{1,1,1});
            potentialMasterLS.addPotential(pLS, new AtomType[]{sim.species.getLeafType(), sim.species.getLeafType()});

            meterSolidLS = new MeterSolidDACut(sim.getSpace(), potentialMasterLS, sim.coordinateDefinition, cutoffsLS);
            meterSolidLS.setTemperature(temperature);
            meterSolidLS.setBPRes(bpharm);
            d = meterSolidLS.getData();

            if (params.ss) {
                potentialMasterLJLS = new PotentialMasterMonatomic(sim);
                pLJLS = new Potential2SoftSphericalLSMultiLat(sim.getSpace(), cutoffsLS, p2LJ, sim.coordinateDefinition);
                potentialMasterLJLS.addPotential(pLJLS, new AtomType[]{sim.species.getLeafType(), sim.species.getLeafType()});
                if (bpharmLJ.length < cutoffsLS.length) {
                    throw new RuntimeException("I need LJ harmonic pressures for all LS cutoffs");
                }
                meterSolidLS.setPotentialMasterDADv2(potentialMasterLJLS, bpharmLJ);
            }

            for (int i=0; i<uFacCut.length; i++) {
                uFacCutLS[i] = uFacCut[i];
            }
            for (int i=uFacCut.length; i<uFacCutLS.length; i++) {
                uFacCutLS[i] = d.getValue(6*i)*numAtoms - uShort;
            }
        }
        System.out.print("cutoffs: ");
        if (nCutoffsLS>0) {
            for (int i=0; i<nCutoffsLS; i++) {
                System.out.print(" "+cutoffsLS[i]);
            }
        }
        else {
            for (int i=0; i<nCutoffs; i++) {
                System.out.print(" "+cutoffs[i]);
            }
        }
        System.out.println();
        System.out.print("bPharm ");
        if (nCutoffsLS>0) {
            for (int i=0; i<nCutoffsLS; i++) {
                System.out.print(" "+bpharm[i]);
            }
        }
        else {
            for (int i=0; i<nCutoffs; i++) {
                System.out.print(" "+bpharm[i]);
            }
        }
        System.out.println();
        if (ss) {
            System.out.print("bPharmLJ ");
            if (nCutoffsLS>0) {
                for (int i=0; i<nCutoffsLS; i++) {
                    System.out.print(" "+bpharmLJ[i]);
                }
            }
            else {
                for (int i=0; i<nCutoffs; i++) {
                    System.out.print(" "+bpharmLJ[i]);
                }
            }
            System.out.println();
        }


        if (args.length == 0) {
            // quick initialization
            sim.initialize(numSteps/10);
        }
        else {
            long nSteps = numSteps/20 + 50*numAtoms + numAtoms*numAtoms*3;
            if (nSteps > numSteps/2) nSteps = numSteps/2;
            sim.initialize(nSteps);
        }

        int numBlocks = 100;
        int interval = 2*numAtoms;
        int intervalLS = 5*interval;
        long blockSize = numSteps/(numBlocks*interval);
        if (blockSize == 0) blockSize = 1;
        long blockSizeLS = numSteps/(numBlocks*intervalLS);
        if (blockSizeLS == 0) blockSizeLS = 1;
        int o=2;
        while (blockSize<numSteps/5 && (numSteps != numBlocks*intervalLS*blockSizeLS || numSteps != numBlocks*interval*blockSize)) {
            interval = 2*numAtoms+(o%2==0 ? (o/2) : -(o/2));
            if (interval < 1 || interval > numSteps/5) {
                throw new RuntimeException("oops interval "+interval);
            }
            // only need to enforce intervalLS if nCutoffsLS>0.  whatever.
            intervalLS = 5*interval;
            blockSize = numSteps/(numBlocks*interval);
            if (blockSize == 0) blockSize = 1;
            blockSizeLS = numSteps/(numBlocks*intervalLS);
            if (blockSizeLS == 0) blockSizeLS = 1;
            o++;
        }
        if (numSteps != numBlocks*intervalLS*blockSizeLS || numSteps != numBlocks*interval*blockSize) {
            throw new RuntimeException("unable to find appropriate intervals");
        }
        System.out.println("block size "+blockSize+" interval "+interval);

        final ValueCache energyFastCache = new ValueCache(meterEnergyShort, sim.integrator);

        DataProcessorReweight puReweight = new DataProcessorReweight(temperature, energyFastCache, uFacCut, sim.box, cutoffs.length);
        DataPumpListener pumpPU = new DataPumpListener(meterSolid, puReweight, interval);
        sim.integrator.getEventManager().addListener(pumpPU);
        final AccumulatorAverageCovariance avgSolid = new AccumulatorAverageCovariance(blockSize);
        puReweight.setDataSink(avgSolid);

        DataProcessorReweightRatio puReweightRatio = new DataProcessorReweightRatio(cutoffs.length);
        avgSolid.setBlockDataSink(puReweightRatio);
        AccumulatorAverageCovariance accPUBlocks = new AccumulatorAverageCovariance(1, true);
        puReweightRatio.setDataSink(accPUBlocks);

        AccumulatorAverageCovariance accPULSBlocks = null;
        AccumulatorAverageCovariance accPULS = null;
        if (nCutoffsLS>0) {
            DataProcessorReweight puLSReweight = new DataProcessorReweight(temperature, energyFastCache, uFacCutLS, sim.box, nCutoffsLS);
            DataPumpListener pumpPULS = new DataPumpListener(meterSolidLS, puLSReweight, intervalLS);
            sim.integrator.getEventManager().addListener(pumpPULS);
            accPULS = new AccumulatorAverageCovariance(blockSizeLS);
            puLSReweight.setDataSink(accPULS);

            DataProcessorReweightRatio puLSReweightRatio = new DataProcessorReweightRatio(nCutoffsLS, nCutoffs-1);
            accPULS.setBlockDataSink(puLSReweightRatio);

            accPULSBlocks = new AccumulatorAverageCovariance(1, true);
            puLSReweightRatio.setDataSink(accPULSBlocks);
        }


        final long startTime = System.currentTimeMillis();

        sim.activityIntegrate.setMaxSteps(numSteps);

        sim.getController().actionPerformed();
        long endTime = System.currentTimeMillis();
        System.out.println();

        IData avgRawData = avgSolid.getData(AccumulatorAverage.AVERAGE);
        IData errRawData = avgSolid.getData(AccumulatorAverage.ERROR);
        IData corRawData = avgSolid.getData(AccumulatorAverage.BLOCK_CORRELATION);
        IData covRawData = avgSolid.getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE);

        int j = 0;
        for (int i=0; i<cutoffs.length; i++) {
            double avgW = avgRawData.getValue(j+6);
            double errW = errRawData.getValue(j+6);
            double corW = corRawData.getValue(j+6);
            System.out.println(String.format("rc: %2d dbA:   % 21.15e  %10.4e  % 5.3f  % 6.4f", i, -Math.log(avgW)/numAtoms, errW/avgW/numAtoms, corW, errW/avgW));
            j += 7;
        }
        System.out.println("\n");

        if (nCutoffsLS>0) {
            avgRawData = accPULS.getData(AccumulatorAverage.AVERAGE);
            errRawData = accPULS.getData(AccumulatorAverage.ERROR);
            corRawData = accPULS.getData(AccumulatorAverage.BLOCK_CORRELATION);

            j = 0;
            for (int i=0; i<cutoffsLS.length; i++) {
                double avgW = avgRawData.getValue(j+6);
                double errW = errRawData.getValue(j+6);
                double corW = corRawData.getValue(j+6);
                System.out.println(String.format("rcLS: %2d dbA:   % 21.15e  %10.4e  % 5.3f  % 6.4f", i, -Math.log(avgW)/numAtoms, errW/avgW/numAtoms, corW, errW/avgW));
                j += 7;
            }
            System.out.println("\n");
        }

        IData avgData = accPUBlocks.getData(AccumulatorAverage.AVERAGE);
        IData errData = accPUBlocks.getData(AccumulatorAverage.ERROR);
        IData corData = accPUBlocks.getData(AccumulatorAverage.BLOCK_CORRELATION);
        IData covData = accPUBlocks.getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE);

        int n = errData.getLength();

        avgRawData = avgSolid.getData(AccumulatorAverage.AVERAGE);
        errRawData = avgSolid.getData(AccumulatorAverage.ERROR);
        covRawData = avgSolid.getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE);

        int jRaw = 0;
        j = 0;
        int nRaw = avgRawData.getLength();
        for  (int i=0; i<cutoffs.length; i++) {
            // compute bias based on Taylor series expansion
            // Xmeasured = Xtrue ( 1 + (ew/w)^2 + (cov_XW)/(X W) )
            double avgW = avgRawData.getValue(jRaw+6);
            double errW = errRawData.getValue(jRaw+6);
            double errWratio2 = errW*errW/(avgW*avgW);

            double avgU = avgRawData.getValue(jRaw+0);
            double errU = errRawData.getValue(jRaw+0);
            double corUW = covRawData.getValue((jRaw+6)*nRaw+(jRaw+0))/Math.sqrt(covRawData.getValue((jRaw+6)*nRaw+(jRaw+6))*covRawData.getValue((jRaw+0)*nRaw+(jRaw+0)));
            if (errW < 1e-7) corUW=0; // if this is sampled rc
            double biasU = errWratio2 - errU*errW/(avgU*avgW)*corUW;
            errU = Math.abs(avgU/avgW)*Math.sqrt(errU*errU/(avgU*avgU) + errWratio2 - 2*errU*errW/(avgU*avgW)*corUW);
            avgU /= avgW;
            double corU = corData.getValue(j+0);

            double avgP = avgRawData.getValue(jRaw+1);
            double errP = errRawData.getValue(jRaw+1);
            double corPW = covRawData.getValue((jRaw+6)*nRaw+(jRaw+1))/Math.sqrt(covRawData.getValue((jRaw+6)*nRaw+(jRaw+6))*covRawData.getValue((jRaw+1)*nRaw+(jRaw+1)));
            if (errW < 1e-7) corPW=0; // if this is sampled rc
            double biasP = errWratio2 - errP*errW/(avgP*avgW)*corPW;
            errP = Math.abs(avgP/avgW)*Math.sqrt(errP*errP/(avgP*avgP) + errWratio2 - 2*errP*errW/(avgP*avgW)*corPW);
            avgP /= avgW;
            double corP = corData.getValue(j+1);

            double avgBUc = avgRawData.getValue(jRaw+2);
            double errBUc = errRawData.getValue(jRaw+2);
            double corBUcW = covRawData.getValue((jRaw+6)*nRaw+(jRaw+2))/Math.sqrt(covRawData.getValue((jRaw+6)*nRaw+(jRaw+6))*covRawData.getValue((jRaw+2)*nRaw+(jRaw+2)));
            if (errW < 1e-7) corBUcW=0; // if this is sampled rc
            double biasBUc = errWratio2 - errBUc*errW/(avgBUc*avgW)*corBUcW;
            errBUc = Math.abs(avgBUc/avgW)*Math.sqrt(errBUc*errBUc/(avgBUc*avgBUc) + errWratio2 - 2*errBUc*errW/(avgBUc*avgW)*corBUcW);
            avgBUc /= avgW;
            double corBUc = corData.getValue(j+2);

            double avgZc = avgRawData.getValue(jRaw+3);
            double errZc = errRawData.getValue(jRaw+3);
            double corZcW = covRawData.getValue((jRaw+6)*nRaw+(jRaw+3))/Math.sqrt(covRawData.getValue((jRaw+6)*nRaw+(jRaw+6))*covRawData.getValue((jRaw+3)*nRaw+(jRaw+3)));
            if (errW < 1e-7) corZcW=0; // if this is sampled rc
            double biasZc = errWratio2 - errZc*errW/(avgZc*avgW)*corZcW;
            errZc = Math.abs(avgZc/avgW)*Math.sqrt(errZc*errZc/(avgZc*avgZc) + errWratio2 - 2*errZc*errW/(avgZc*avgW)*corZcW);
            avgZc /= avgW;
            double corZc = corData.getValue(j+3);

            // this is dbAc/dv2 at constant Y (for LJ)
            double avgDADv2 = avgRawData.getValue(jRaw+4);
            double errDADv2 = errRawData.getValue(jRaw+4);
            double corDADv2W = covRawData.getValue((jRaw+6)*nRaw+(jRaw+4))/Math.sqrt(covRawData.getValue((jRaw+6)*nRaw+(jRaw+6))*covRawData.getValue((jRaw+4)*nRaw+(jRaw+4)));
            if (errW < 1e-7) corDADv2W=0; // if this is sampled rc
            double biasDADv2 = errWratio2 - errDADv2*errW/(avgDADv2*avgW)*corDADv2W;
            errDADv2 = Math.abs(avgDADv2/avgW)*Math.sqrt(errDADv2*errDADv2/(avgDADv2*avgDADv2) + errWratio2 - 2*errDADv2*errW/(avgDADv2*avgW)*corDADv2W);
            avgDADv2 /= avgW;
            double corDADv2 = corData.getValue(j+4);

            double avgPcZ = avgRawData.getValue(jRaw+5);
            double errPcZ = errRawData.getValue(jRaw+5);
            double corPcZW = covRawData.getValue((jRaw+6)*nRaw+(jRaw+5))/Math.sqrt(covRawData.getValue((jRaw+6)*nRaw+(jRaw+6))*covRawData.getValue((jRaw+5)*nRaw+(jRaw+5)));
            if (errW < 1e-7) corPcZW=0; // if this is sampled rc
            double biasPcZ = errWratio2 - errPcZ*errW/(avgPcZ*avgW)*corPcZW;
            errPcZ = Math.abs(avgPcZ/avgW)*Math.sqrt(errPcZ*errPcZ/(avgPcZ*avgPcZ) + errWratio2 - 2*errPcZ*errW/(avgPcZ*avgW)*corPcZW);
            avgPcZ /= avgW;
            double corPcZ = corData.getValue(j+5);

            double facDADY = 4*density*density*density*density/temperature;
            double PUCor = covData.getValue((j+1)*n+(j+0))/Math.sqrt(covData.getValue((j+1)*n+(j+1))*covData.getValue((j+0)*n+(j+0)));
            double DADACor = -covData.getValue((j+2)*n+(j+4))/Math.sqrt(covData.getValue((j+2)*n+(j+2))*covData.getValue((j+4)*n+(j+4)));
            double ZcUcCor = covData.getValue((j+2)*n+(j+3))/Math.sqrt(covData.getValue((j+2)*n+(j+2))*covData.getValue((j+3)*n+(j+3)));

            System.out.print(String.format("rc: %2d DADY:  % 21.15e  %10.4e  % 10.4e  % 5.3f\n", i, -facDADY*avgBUc, facDADY*errBUc, -facDADY*avgBUc*biasBUc, corBUc));
            System.out.print(String.format("rc: %2d DADv2: % 21.15e  %10.4e  % 10.4e  % 5.3f  % 8.6f\n", i, avgDADv2, errDADv2, avgDADv2*biasDADv2, corDADv2, DADACor));
            System.out.print(String.format("rc: %2d Zc:    % 21.15e  %10.4e  % 10.4e  % 5.3f\n", i, avgZc, errZc, avgZc*biasZc, corZc));
            System.out.print(String.format("rc: %2d bUc:   % 21.15e  %10.4e  % 10.4e  % 5.3f  % 8.6f\n", i, avgBUc, errBUc, avgBUc*biasBUc, corBUc, ZcUcCor));
            System.out.print(String.format("rc: %2d Uraw:  % 21.15e  %10.4e  % 10.4e  % 5.3f\n", i, avgU, errU, avgU*biasU, corU));
            System.out.print(String.format("rc: %2d Praw:  % 21.15e  %10.4e  % 10.4e  % 5.3f  % 8.6f\n", i, avgP, errP, avgP*biasP, corP, PUCor));
            System.out.print(String.format("rc: %2d PcZ:   % 21.15e  %10.4e  % 10.4e  % 5.3f\n", i, avgPcZ, errPcZ, avgPcZ*biasPcZ, corPcZ));
            System.out.println();
            j+=6;
            jRaw+=7;
        }

        if (nCutoffsLS > 0) {

            avgRawData = accPULS.getData(AccumulatorAverage.AVERAGE);

            avgData = accPULSBlocks.getData(AccumulatorAverage.AVERAGE);
            errData = accPULSBlocks.getData(AccumulatorAverage.ERROR);
            covData = accPULSBlocks.getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE);
            corData = accPULSBlocks.getData(AccumulatorAverage.BLOCK_CORRELATION);

            n = errData.getLength();

            j = 0;
            jRaw = 0;
            double avgUref = 0, avgPref = 0, avgZCref = 0, avgBUCref = 0, avgDADv2ref = 0, avgPcZref = 0;
            for (int i=0; i<cutoffsLS.length; i++) {
                if (i<cutoffs.length-1) {
                    j+=6;
                    jRaw+=7;
                    continue;
                }
                // estimate bias based on difference between averages using all data and average of block data
                // that gives us bias in block data.  then divide by number of blocks (bias proportional to variance and covariance)
                double avgW = avgRawData.getValue(jRaw+6);
                double avgU = avgRawData.getValue(jRaw+0)/avgW;
                double avgU1 = avgData.getValue(j+0);
                double errU = errData.getValue(j+0);
                double corU = corData.getValue(j+0);
                double avgP = avgRawData.getValue(jRaw+1)/avgW;
                double avgP1 = avgData.getValue(j+1);
                double errP = errData.getValue(j+1);
                double corP = corData.getValue(j+1);
                double avgBUc = avgRawData.getValue(jRaw+2)/avgW;
                double avgBUc1 = avgData.getValue(j+2);
                double errBUc = errData.getValue(j+2);
                double corBUc = corData.getValue(j+2);
                double avgZc = avgRawData.getValue(jRaw+3)/avgW;
                double avgZc1 = avgData.getValue(j+3);
                double errZc = errData.getValue(j+3);
                double corZc = corData.getValue(j+3);
                // this is dbAc/drho at constant Y (for LJ)
                double avgDADv2 = avgRawData.getValue(jRaw+4)/avgW;
                double avgDADv21 = avgData.getValue(j+4);
                double errDADv2 = errData.getValue(j+4);
                double corDADv2 = corData.getValue(j+4);
                double avgPcZ = avgRawData.getValue(jRaw+5)/avgW;
                double avgPcZ1 = avgData.getValue(j+5);
                double errPcZ = errData.getValue(j+5);
                double corPcZ = corData.getValue(j+5);

                double DADACor = -covData.getValue((j+2)*n+(j+4))/Math.sqrt(covData.getValue((j+2)*n+(j+2))*covData.getValue((j+4)*n+(j+4)));
                double ZcUcCor = covData.getValue((j+2)*n+(j+3))/Math.sqrt(covData.getValue((j+2)*n+(j+2))*covData.getValue((j+3)*n+(j+3)));
                double facDADY = 4*density*density*density*density/temperature;
                double PUCor = covData.getValue((j+1)*n+(j+0))/Math.sqrt(covData.getValue((j+1)*n+(j+1))*covData.getValue((j+0)*n+(j+0)));

                if (i==cutoffs.length-1) {
                    avgUref = avgU;
                    avgPref = avgP;
                    avgZCref = avgZc;
                    avgBUCref = avgBUc;
                    avgDADv2ref = avgDADv2;
                    avgPcZref = avgPcZ;
                    j+=6;
                    jRaw+=7;
                    continue;
                }
                avgU -= avgUref;
                avgP -= avgPref;
                avgZc -= avgZCref;
                avgBUc -= avgBUCref;
                avgDADv2 -= avgDADv2ref;
                avgPcZ -= avgPcZref;

                System.out.print(String.format("rcLS: %2d DADY:  % 21.15e  %10.4e  % 11.4e  % 5.3f\n", i, -facDADY*avgBUc, facDADY*errBUc, -facDADY*(avgBUc1-avgBUc)/numBlocks, corBUc));
                System.out.print(String.format("rcLS: %2d DADv2: % 21.15e  %10.4e  % 11.4e  % 5.3f  % 8.6f\n", i, avgDADv2, errDADv2, (avgDADv21-avgDADv2)/numBlocks, corDADv2, DADACor));
                System.out.print(String.format("rcLS: %2d Zc:    % 21.15e  %10.4e  % 11.4e  % 5.3f\n", i, avgZc, errZc, (avgZc1-avgZc)/numBlocks, corZc));
                System.out.print(String.format("rcLS: %2d bUc:   % 21.15e  %10.4e  % 11.4e  % 5.3f  % 8.6f\n", i, avgBUc, errBUc, (avgBUc1-avgBUc)/numBlocks, corBUc, ZcUcCor));
                System.out.print(String.format("rcLS: %2d Uraw:  % 21.15e  %10.4e  % 11.4e  % 5.3f\n", i, avgU, errU, (avgU1-avgU)/numBlocks, corU));
                System.out.print(String.format("rcLS: %2d Praw:  % 21.15e  %10.4e  % 11.4e  % 5.3f  % 8.6f\n", i, avgP, errP, (avgP1-avgP)/numBlocks, corP, PUCor));
                System.out.print(String.format("rcLS: %2d PcZ:   % 21.15e  %10.4e  % 11.4e  % 5.3f\n", i, avgPcZ, errPcZ, (avgPcZ1-avgPcZ)/numBlocks, corPcZ));
                System.out.println();

                j+=6;
                jRaw+=7;
            }
        }

        System.out.println("time: " + (endTime - startTime)/1000.0);
    }

    public void initialize(long initSteps) {
        // equilibrate off the lattice to avoid anomalous contributions
        activityIntegrate.setMaxSteps(initSteps);
        getController().actionPerformed();
        getController().reset();
        integrator.getMoveManager().setEquilibrating(false);
    }
    
    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numAtoms = 256;
        public double density = 1.28;
        public long numSteps = 1000000;
        public double temperature = 0.1;
        public double rc = 2.5;
        public double rc0 = rc;
        public double rcMax1 = 100;
        public double rcMax0 = 100;
        public double[] bpharm = new double[0];
        public double[] bpharmLJ = new double[0];
        public boolean ss = false;
        public int[] randomSeeds = null;
        public double alpha = 1.0;
    }
}
