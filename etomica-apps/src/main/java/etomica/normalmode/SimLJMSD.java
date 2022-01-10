/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;

import java.awt.*;


public class SimLJMSD extends Simulation {

    public final CoordinateDefinitionLeaf coordinateDefinition;
    public IntegratorMC integrator;

    public Box box;
    public Boundary boundary;
    public int[] nCells;
    public Basis basis;
    public Primitive primitive;
    public MCMoveAtomCoupled atomMove;
    public PotentialMasterList potentialMaster;
    public Potential2SoftSpherical potential;
    public SpeciesGeneral species;
    MeterPotentialEnergy meterPE;

    public SimLJMSD(Space _space, int numAtoms, double density, double temperature, double rc, boolean ss) {
        super(_space);
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        potentialMaster = new PotentialMasterList(this, space);

        // TARGET
        double L = Math.pow(4.0 / density, 1.0 / 3.0);
        int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
        boundary = new BoundaryRectangularPeriodic(space, n * L);
        box = this.makeBox(boundary);
        box.setNMolecules(species, numAtoms);

        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature, box);
        meterPE = new MeterPotentialEnergy(potentialMaster, box);
        atomMove = new MCMoveAtomCoupled(potentialMaster, meterPE, getRandom(), space);
        atomMove.setStepSize(0.1);
        atomMove.setStepSizeMax(0.5);
        atomMove.setDoExcludeNonNeighbors(true);
        integrator.getMoveManager().addMCMove(atomMove);
//        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);

        primitive = new PrimitiveCubic(space, n * L);

        nCells = new int[]{n, n, n};
        Basis basisFCC = new BasisCubicFcc();
        basis = new BasisBigCell(space, basisFCC, nCells);

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});

        potential = ss ? new P2SoftSphere(space, 1.0, 4.0, 12) : new P2LennardJones(space, 1.0, 1.0);
        potential = new P2SoftSphericalTruncated(space, potential, rc);
        atomMove.setPotential(potential);
        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(potential, new AtomType[]{sphereType, sphereType});

        potentialMaster.lrcMaster().setEnabled(false);

        int cellRange = 2;
        potentialMaster.setRange(rc);
        potentialMaster.setCellRange(cellRange);
        potentialMaster.getNeighborManager(box).reset();
        int potentialCells = potentialMaster.getNbrCellManager(box).getLattice().getSize()[0];
        if (potentialCells < cellRange * 2 + 1) {
            throw new RuntimeException("oops (" + potentialCells + " < " + (cellRange * 2 + 1) + ")");
        }
        this.getController().addActivity(new ActivityIntegrate(integrator));
        ((P2SoftSphericalTruncated) potential).setTruncationRadius(0.6 * boundary.getBoxSize().getX(0));
    }

    /**
     * @param args filename containing simulation parameters
     */
    public static void main(String[] args) {
        //set up simulation parameters
        SimOverlapParam params = new SimOverlapParam();

        boolean ss = params.ss;
        double density = params.density;
        long numSteps = params.numSteps;
        final int numAtoms = params.numAtoms;
        double temperature = params.temperature;
        double rc = params.rc;

        System.out.println("Running "+(ss?"soft-sphere":"Lennard-Jones")+" simulation");
        System.out.println(numAtoms+" atoms at density "+density+" and temperature "+temperature);
        System.out.println(numSteps+" steps");

        //instantiate simulation
        final SimLJMSD sim = new SimLJMSD(Space.getInstance(3), numAtoms, density, temperature, rc, ss);
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

            simGraphic.makeAndDisplayFrame((ss?"SS":"LJ")+" FCC");

            return;
        }

        //U
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster, sim.box);
        double ULat = meterPE.getDataAsScalar();

        System.out.println(" Lattice Energy/N = " + (ULat/numAtoms));

        //Mapped
        MeterSolidPropsLJ meterSolid = new  MeterSolidPropsLJ(sim.space, new MeterPotentialEnergyFromIntegrator(sim.integrator),  sim.potentialMaster, sim.coordinateDefinition, temperature, new double[11]);

//Initialization
        System.out.flush();
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


        //U
        AccumulatorAverageFixed accumulatorPE = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPEPump = new DataPumpListener(meterPE, accumulatorPE, interval);
        sim.integrator.getEventManager().addListener(accumulatorPEPump);

        //Mapped
        AccumulatorAverageFixed accumulatorUP = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorUPPump = new DataPumpListener(meterSolid, accumulatorUP, interval);
        sim.integrator.getEventManager().addListener(accumulatorUPPump);

        final long startTime = System.currentTimeMillis();

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));

        System.out.println();
        //U
        IData dataPEAvg = accumulatorPE.getData(accumulatorPE.AVERAGE);
        IData dataPEErr = accumulatorPE.getData(accumulatorPE.ERROR);
        IData dataPECorrelation = accumulatorPE.getData(accumulatorPE.BLOCK_CORRELATION);
        double peAvg = dataPEAvg.getValue(0);
        double peErr = dataPEErr.getValue(0);
        double peCor = dataPECorrelation.getValue(0);
        System.out.println(" Udirect/N = "+peAvg/numAtoms+"  Err = "+peErr/numAtoms+"  cor: "+peCor);

        IData dataUPAvg = accumulatorUP.getData(accumulatorUP.AVERAGE);
        IData dataUPErr = accumulatorUP.getData(accumulatorUP.ERROR);
        IData dataUPCorr = accumulatorUP.getData(accumulatorUP.BLOCK_CORRELATION);

        double u2   = dataUPAvg.getValue(0);
        double u2Err = dataUPErr.getValue(0);
        double u2Corr = dataUPCorr.getValue(0);

        double u2_alpha   = dataUPAvg.getValue(1);
        double u2_alphaErr = dataUPErr.getValue(1);
        double u2_alphaCorr = dataUPCorr.getValue(1);

        System.out.println();
        System.out.println("************************************************************************");
        System.out.println("1 0 0");
        System.out.println("u2        : " + u2 + " err: " + u2Err + " corr: " + u2Corr);
        System.out.println("u2_alpha  : " + u2_alpha + " err: " + u2_alphaErr + " corr: " + u2_alphaCorr);


        long endTime = System.currentTimeMillis();
        System.out.println("time: " + (endTime - startTime)/1000.0);
    }



    public void initialize(long initSteps) {
        // equilibrate off the lattice to avoid anomalous contributions
        this.getController().runActivityBlocking(new ActivityIntegrate(this.integrator, initSteps));

        integrator.getMoveManager().setEquilibrating(false);
    }
    
    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numAtoms = 500;
        public double density = 1.2;
        public long numSteps = 10000000;
        public double temperature = 0.1;
        public double rc = 3.0;
        public boolean ss = false;
    }
}