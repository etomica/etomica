/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;


import etomica.action.activity.ActivityIntegrate2;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressure;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBoxSize;
import etomica.lattice.crystal.*;
import etomica.math.function.Function;
import etomica.nbr.list.NeighborListManagerSlanty;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;

/**
 * Simulation to run sampling with the hard sphere potential, but measuring
 * the harmonic potential based on normal mode data from a previous simulation.
 * 
 * The original Bennett's Overlapping Sampling Simulation
 * 	- used to check for the computation time
 * 
 * @author Tai Boon Tan
 */
public class SimOverlapSoftSphereTP extends Simulation {

    public IntegratorMC integrator;

    public Box box;
    public Boundary boundary;
    public CoordinateDefinitionLeaf coordinateDefinition;
    public int[] nCells;
    public Basis basis;
    public Primitive primitive;
    public AccumulatorAverageFixed accumulator;
    public DataPumpListener accumulatorPump;
    public MeterTargetTP meter;
    protected MCMoveAtomCoupled atomMove;
    protected PotentialMaster potentialMaster;
    protected double latticeEnergy;
    public SimOverlapSoftSphereTP(Space _space, int numAtoms, boolean slanty, boolean flex, double density, double temperature, double[] otherTemperatures, double[] alpha, int exponent, int numAlpha, double alphaSpan, long numSteps, double rc) {
        super(_space);

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        if (slanty) {
            potentialMaster = new PotentialMasterList(this, rc, new NeighborListManagerSlanty.NeighborListSlantyAgentSource(rc), space);
        } else {
            potentialMaster = new PotentialMasterList(this, space);
        }

        // TARGET
        double nbrDistance = 0;
        if (slanty) {
            int c = (int) Math.round(Math.pow(numAtoms, 1.0 / 3.0));
            nCells = new int[]{c, c, c};

            double L = Math.pow(Math.sqrt(2) / density, 1.0 / 3.0);
            nbrDistance = L;
            double angle = Math.PI / 3;

//            primitive = new PrimitiveFcc(space, L*c);
            primitive = new PrimitiveTriclinic(space, L * c, L * c, L * c, angle, angle, angle);

            boundary = new BoundaryDeformablePeriodic(space, primitive.vectors());
            ((BoundaryDeformablePeriodic) boundary).setTruncationRadius(rc);
            Basis basisSimple = new Basis(new Vector3D[]{new Vector3D(0.0, 0.0, 0.0)});
            basis = new BasisBigCell(space, basisSimple, nCells);
        } else {

            double L = Math.pow(4.0 / density, 1.0 / 3.0);
            nbrDistance = L / Math.sqrt(2);
            int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
            primitive = new PrimitiveCubic(space, n * L);

            nCells = new int[]{n, n, n};
            boundary = new BoundaryRectangularPeriodic(space, n * L);
            Basis basisFCC = new BasisCubicFcc();
            basis = new BasisBigCell(space, basisFCC, nCells);
        }
        box = this.makeBox(boundary);
        box.setNMolecules(species, numAtoms);

        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature, box);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster, box);
        atomMove = new MCMoveAtomCoupled(potentialMaster, meterPE, getRandom(), space);
        atomMove.setStepSize(0.1);
        atomMove.setStepSizeMax(0.5);
        atomMove.setDoExcludeNonNeighbors(true);
        integrator.getMoveManager().addMCMove(atomMove);
//        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});

        Potential2SoftSpherical potential = new P2SoftSphere(space, 1.0, 1.0, exponent);
        if (potentialMaster instanceof PotentialMasterList) {
            potential = new P2SoftSphericalTruncated(space, potential, rc);

        } else {
            potential = new P2SoftSphericalTruncatedShifted(space, potential, rc);

        }
        atomMove.setPotential(potential);
        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(potential, new AtomType[]{sphereType, sphereType});

        if (flex) {
            MCMoveBoxSize moveBoxSize = new MCMoveBoxSize(potentialMaster, random, space);
            integrator.getMoveManager().addMCMove(moveBoxSize);
        }

        /*
         *  1-body Potential to Constraint the atom from moving too far
         *  	away from its lattice-site
         *
         */

        P1ConstraintNbr p1Constraint = new P1ConstraintNbr(space, nbrDistance);
        p1Constraint.initBox(box);
        atomMove.setConstraint(p1Constraint);

        potentialMaster.lrcMaster().setEnabled(false);

        if (potentialMaster instanceof PotentialMasterList) {
            int cellRange = 7;
            ((PotentialMasterList) potentialMaster).setRange(rc);
            ((PotentialMasterList) potentialMaster).setCellRange(cellRange); // insanely high, this lets us have neighborRange close to dimensions/2
            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
            ((PotentialMasterList) potentialMaster).getNeighborManager(box).reset();
            int potentialCells = ((PotentialMasterList) potentialMaster).getNbrCellManager(box).getLattice().getSize()[0];
            if (potentialCells < cellRange * 2 + 1) {
                throw new RuntimeException("oops (" + potentialCells + " < " + (cellRange * 2 + 1) + ")");
            }
        }

        // this is the meterPE we gave to the MCMove but it hasn't called setTarget yet, so we're OK
        latticeEnergy = meterPE.getDataAsScalar();

        meter = new MeterTargetTP(potentialMaster, species, this, coordinateDefinition);
        meter.setLatticeEnergy(latticeEnergy);
        meter.setTemperature(temperature);
        meter.setOtherTemperatures(otherTemperatures);
        meter.setAlpha(alpha);
        meter.setAlphaSpan(alphaSpan);
        meter.setNumAlpha(numAlpha);
        meter.setConstraint(p1Constraint);
        int numBlocks = 100;
        int interval = numAtoms;
        long blockSize = numSteps / (numBlocks * interval);
        if (blockSize == 0) blockSize = 1;
        System.out.println("block size " + blockSize + " interval " + interval);
        if (otherTemperatures.length > 1) {
            accumulator = new AccumulatorAverageCovariance(blockSize);
        } else {
            accumulator = new AccumulatorAverageFixed(blockSize);
        }
        accumulatorPump = new DataPumpListener(meter, accumulator, interval);
        integrator.getEventManager().addListener(accumulatorPump);

        this.getController2().addActivity(new ActivityIntegrate2(integrator));

        if (potentialMaster instanceof PotentialMasterList) {
            // extend potential range, so that atoms that move outside the truncation range will still interact
            // atoms that move in will not interact since they won't be neighbors
            ((P2SoftSphericalTruncated) potential).setTruncationRadius(0.6 * boundary.getBoxSize().getX(0));
        }
    }

    /**
     * @param args filename containing simulation parameters
     * @see SimOverlapSoftSphereTP.SimOverlapParam
     */
    public static void main(String[] args) {
        //set up simulation parameters
        SimOverlapParam params = new SimOverlapParam();
        ParseArgs.doParseArgs(params, args);
        if (args.length==0) {
            params.alpha = new double[0];
            params.otherTemperatures = new double[0];
            params.numSteps = 10000000;
            params.temperature = 0.5;
            params.pEst = 9.345;
        }
        final double density = params.density;
        int exponentN = params.exponentN;
        long numSteps = params.numSteps;
        final int numAtoms = params.numAtoms;
        boolean slanty = params.slanty;
        boolean flex = params.flex;
        double temperature = params.temperature;
        double[] otherTemperatures = params.otherTemperatures;
        double[] alpha = params.alpha;
        int numAlpha = params.numAlpha;
        double alphaSpan = params.alphaSpan;
        double rc = params.rc;
        double pEst = params.pEst;

        System.out.println("Running soft sphere overlap simulation");
        System.out.println(numAtoms+" atoms at density "+density+" and temperature "+temperature);
        System.out.println("exponent N: "+ exponentN);
        System.out.println(numSteps+" steps");

        //instantiate simulation
        final SimOverlapSoftSphereTP sim = new SimOverlapSoftSphereTP(Space.getInstance(3), numAtoms, slanty, flex, density, temperature, otherTemperatures, alpha, exponentN, numAlpha, alphaSpan, numSteps, rc);
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

//        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster);
//        meterPE.setBox(sim.box);
//        meterPE.setIncludeLrc(false);
//        IVectorMutable l = sim.space.makeVector();
//        IVectorMutable l0 = sim.space.makeVector();
//        l0.E(sim.box.getBoundary().getBoxSize());
//        BoxInflate boxInflater = new BoxInflate(sim.space);
//        boxInflater.setBox(sim.box);
//        for (int i=-10; i<11; i++) {
//            double lntheta = i/100.0;
//            double theta = Math.exp(lntheta);
//            l.setX(0, Math.sqrt(theta));
//            l.setX(1, 1.0/Math.sqrt(theta));
//            l.setX(2, 1);
//            boxInflater.setVectorScale(l);
//            boxInflater.actionPerformed();
//            double u1 = (meterPE.getDataAsScalar()-sim.latticeEnergy);
//            boxInflater.undo();
//            l.setX(1, 1.0/Math.pow(theta, 0.25));
//            l.setX(2, l.getX(1));
//            boxInflater.setVectorScale(l);
//            boxInflater.actionPerformed();
//            double u2 = (meterPE.getDataAsScalar()-sim.latticeEnergy);
//            boxInflater.undo();
//            System.out.println(lntheta+" "+u1+" "+u2);
//        }
//        System.exit(0);

        //start simulation

        sim.initialize(numAtoms*50);
        System.out.flush();

        MeterPressure meterPressure = new MeterPressure(sim.getSpace());
        meterPressure.setIntegrator(sim.integrator);
        AccumulatorAverageFixed accP = new AccumulatorAverageFixed(sim.accumulator.getBlockSize());
        DataPumpListener pumpP = new DataPumpListener(meterPressure, accP, numAtoms);
        sim.integrator.getEventManager().addListener(pumpP);

        /*
        MeterSolidPressure meterSolidPressure = new MeterSolidPressure(sim.getSpace(), sim.potentialMaster, sim.coordinateDefinition);
        double pLattice = 12/(3*sim.box.getBoundary().volume())*sim.latticeEnergy;
        System.out.println("pLattice: "+pLattice);
        meterSolidPressure.setPressureResidual(pEst, pLattice);
        meterSolidPressure.setTemperature(temperature);
        AccumulatorAverageFixed accSolidP = new AccumulatorAverageFixed(sim.accumulator.getBlockSize());
        DataPumpListener pumpSolidP = new DataPumpListener(meterSolidPressure, accSolidP, numAtoms);
        sim.integrator.getEventManager().addListener(pumpSolidP);

        MeterSolidPressure meterSolidPressureLow = new MeterSolidPressure(sim.getSpace(), sim.potentialMaster, sim.coordinateDefinition);
        meterSolidPressureLow.setPressureEstimate(pEst-1, pLattice);
        meterSolidPressureLow.setTemperature(temperature);
        AccumulatorAverageFixed accSolidPLow = new AccumulatorAverageFixed(sim.accumulator.getBlockSize());
        DataPumpListener pumpSolidPLow = new DataPumpListener(meterSolidPressureLow, accSolidPLow, numAtoms);
        sim.integrator.getEventManager().addListener(pumpSolidPLow);

        MeterSolidPressure meterSolidPressureHigh = new MeterSolidPressure(sim.getSpace(), sim.potentialMaster, sim.coordinateDefinition);
        meterSolidPressureHigh.setPressureEstimate(pEst+1, pLattice);
        meterSolidPressureHigh.setTemperature(temperature);
        AccumulatorAverageFixed accSolidPHigh = new AccumulatorAverageFixed(sim.accumulator.getBlockSize());
        DataPumpListener pumpSolidPHigh = new DataPumpListener(meterSolidPressureHigh, accSolidPHigh, numAtoms);
        sim.integrator.getEventManager().addListener(pumpSolidPHigh);*/

        MeterDDP meterDP = new MeterDDP(sim.potentialMaster, sim.getSpecies(0), sim.getSpace(), sim);
        meterDP.setNumP(1);
        meterDP.setCoordinateDefinition(sim.coordinateDefinition);
        meterDP.setOtherDensities(new double[]{density-0.00001, density+0.00001});
        meterDP.setPCenter(new double[]{pEst, pEst});
        meterDP.setTemperature(temperature);
        meterDP.setULatFunction(new Function() {
            public double f(double x) {
                return Math.pow(x/density, 4)*sim.latticeEnergy/numAtoms;
            }
        });
        AccumulatorAverageFixed accDP = new AccumulatorAverageFixed(sim.accumulator.getBlockSize());
        DataPumpListener pumpDP = new DataPumpListener(meterDP, null, numAtoms);
        sim.integrator.getEventManager().addListener(pumpDP);


        final long startTime = System.currentTimeMillis();

        sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integrator), numSteps);

        //MeterTargetTP.openFW("x"+numMolecules+".dat");
        //MeterTargetTP.closeFW();

        double avgP = accP.getData(accP.AVERAGE).getValue(0);
        double errP = accP.getData(accP.ERROR).getValue(0);
        double corP = accP.getData(accP.BLOCK_CORRELATION).getValue(0);
        System.out.println(String.format("pressure: %20.15e %10.4e %3.2f\n", avgP, errP, corP));

        /*
        double avgSP = accSolidP.getData(accSolidP.AVERAGE).getValue(0);
        double errSP = accSolidP.getData(accSolidP.ERROR).getValue(0);
        double corSP = accSolidP.getData(accSolidP.BLOCK_CORRELATION).getValue(0);
        System.out.println(String.format("solid pressure: %20.15e %10.4e %3.2f\n", avgSP, errSP, corSP));

        double avgSPL = accSolidPLow.getData(accSolidPLow.AVERAGE).getValue(0);
        double errSPL = accSolidPLow.getData(accSolidPLow.ERROR).getValue(0);
        double corSPL = accSolidPLow.getData(accSolidPLow.BLOCK_CORRELATION).getValue(0);
        System.out.println(String.format("solid pressure (low): %20.15e %10.4e %3.2f\n", avgSPL, errSPL, corSPL));

        double avgSPH = accSolidPHigh.getData(accSolidPHigh.AVERAGE).getValue(0);
        double errSPH = accSolidPHigh.getData(accSolidPHigh.ERROR).getValue(0);
        double corSPH = accSolidPHigh.getData(accSolidPHigh.BLOCK_CORRELATION).getValue(0);
        System.out.println(String.format("solid pressure (high): %20.15e %10.4e %3.2f\n", avgSPH, errSPH, corSPH));
        */

        System.out.println("\nratio averages:\n");

        DataGroup data = (DataGroup)sim.accumulator.getData();
        IData dataErr = data.getData(sim.accumulator.ERROR.index);
        IData dataAvg = data.getData(sim.accumulator.AVERAGE.index);
        IData dataCorrelation = data.getData(sim.accumulator.BLOCK_CORRELATION.index);
        for (int i=0; i<otherTemperatures.length; i++) {
            System.out.println(otherTemperatures[i]);
            double[] iAlpha = sim.meter.getAlpha(i);
            for (int j=0; j<numAlpha; j++) {
                System.out.println("  "+iAlpha[j]+" "+dataAvg.getValue(i*numAlpha+j)
                        +" "+dataErr.getValue(i*numAlpha+j)
                        +" "+dataCorrelation.getValue(i*numAlpha+j));
            }
        }

        if (otherTemperatures.length == 2) {
            // we really kinda want each covariance for every possible pair of alphas,
            // but we're going to be interpolating anyway and the covariance is almost
            // completely insensitive to choice of alpha.  so just take the covariance for
            // the middle alphas.
            IData dataCov = data.getData(((AccumulatorAverageCovariance) sim.accumulator).BLOCK_COVARIANCE.index);
            System.out.print("covariance "+otherTemperatures[1]+" / "+otherTemperatures[0]+"   ");
            for (int i=0; i<numAlpha; i++) {
                i = (numAlpha-1)/2;
                double ivar = dataErr.getValue(i);
                ivar *= ivar;
                ivar *= (sim.accumulator.getBlockCount()-1);
                for (int j=0; j<numAlpha; j++) {
                    j = (numAlpha-1)/2;
                    double jvar = dataErr.getValue(numAlpha+j);
                    jvar *= jvar;
                    jvar *= (sim.accumulator.getBlockCount()-1);
                    System.out.print(dataCov.getValue(i*(2*numAlpha)+(numAlpha+j))/Math.sqrt(ivar*jvar)+" ");
                    break;
                }
                System.out.println();
                break;
            }
        }

        long endTime = System.currentTimeMillis();
        System.out.println("Time taken: " + (endTime - startTime)/1000.0);
    }

    public void initialize(long initSteps) {
        // equilibrate off the lattice to avoid anomolous contributions
        this.getController2().runActivityBlocking(new ActivityIntegrate2(this.integrator), initSteps);


        accumulator.reset();



    }
    
    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numAtoms = 256;
        public boolean slanty = false;
        public boolean flex = false;
        public double density = 1.0;
        public int exponentN = 12;
        public long numSteps = 1000000;
        public double temperature = 1.;
        public double[] alpha = new double[]{0.011};
        public int numAlpha = 11;
        public double alphaSpan = 1;
        public double[] otherTemperatures = new double[]{1.01};
        public double rc = 2.2;
        public double pEst = 0;
    }
}
