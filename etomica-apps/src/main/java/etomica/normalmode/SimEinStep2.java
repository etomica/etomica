/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.lattice.crystal.*;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeField;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;

/**
 * Simulation that samples a composite energy function (soft sphere and
 * Einstein crystal) and perturbs into overlap regions shared by systems with
 * more or less soft sphere contributions.
 * 
 * @author Andrew Schultz
 */
public class SimEinStep2 extends Simulation {

    public final PotentialMasterList potentialMaster;
    public final PotentialComputeField potentialMasterHarmonic;
    public IntegratorMC integrator;

    public Box box;
    public Boundary boundary;
    public int[] nCells;
    public Basis basis;
    public Primitive primitive;
    public MCMoveBoxStep atomMove;
    public SimEinStep2(Space _space, int numAtoms, double density, double temperature, double lambda, int exponent, double rc, boolean slanty) {
        super(_space);

        SpeciesGeneral species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        // TARGET
        if (slanty) {
            int c = (int) Math.round(Math.pow(numAtoms, 1.0 / 3.0));
            nCells = new int[]{c, c, c};

            double L = Math.pow(Math.sqrt(2) / density, 1.0 / 3.0);
            double angle = Math.PI / 3;

//            primitive = new PrimitiveFcc(space, L*c);
            primitive = new PrimitiveTriclinic(space, L * c, L * c, L * c, angle, angle, angle);

            boundary = new BoundaryDeformablePeriodic(space, primitive.vectors());
            ((BoundaryDeformablePeriodic) boundary).setTruncationRadius(rc);
            Basis basisSimple = new Basis(new Vector3D[]{new Vector3D(0.0, 0.0, 0.0)});
            basis = new BasisBigCell(space, basisSimple, nCells);
        } else {

            double L = Math.pow(4.0 / density, 1.0 / 3.0);
            int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
            primitive = new PrimitiveCubic(space, n * L);

            nCells = new int[]{n, n, n};
            boundary = new BoundaryRectangularPeriodic(space, n * L);
            Basis basisFCC = new BasisCubicFcc();
            basis = new BasisBigCell(space, basisFCC, nCells);
        }
        box = this.makeBox(boundary);
        box.setNMolecules(species, numAtoms);

        int cellRange = 7; // insanely high, this lets us have neighborRange close to dimensions/2
        potentialMaster = new PotentialMasterList(getSpeciesManager(), box, cellRange, rc, BondingInfo.noBonding());
        potentialMaster.doAllTruncationCorrection = false;
        potentialMaster.doOneTruncationCorrection = false;
        potentialMaster.setDoDownNbrs(true);

        // giving potentialMaster to integrator instead of potentialMasterComposite isn't quite right, but will work
        // fine and allows us to remove the listener (so that neighbors won't get updated)
        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature, box);
        integrator.getEventManager().removeListener(potentialMaster);

        CoordinateDefinitionLeaf coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});

        IPotential2 potential = null;
        if (exponent > 0) {
            potential = new P2SoftSphere(1.0, 1.0, exponent);
        } else {
            potential = new P2LennardJones(1.0, 1.0);
        }
        potential = new P2SoftSphericalTruncated(potential, rc);
        AtomType sphereType = species.getLeafType();
        potentialMaster.setPairPotential(sphereType, sphereType, potential);

        // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
        potentialMaster.init();

        P1HarmonicSite p1Harmonic = new P1HarmonicSite(space, coordinateDefinition.siteManager);
        p1Harmonic.setSpringConstant(lambda > 0 ? lambda : 1);
        potentialMasterHarmonic = new PotentialComputeField(getSpeciesManager(), box);
        potentialMasterHarmonic.setFieldPotential(sphereType, p1Harmonic);

        PotentialComputeAggregate potentialComputeComposite = new PotentialComputeAggregate(potentialMasterHarmonic, potentialMaster);

        atomMove = new MCMoveAtomCoupled(lambda == 0 ? potentialMaster : potentialComputeComposite, random, space);
        atomMove.setStepSize(0.1);
        atomMove.setStepSizeMax(0.5);
        integrator.getMoveManager().addMCMove(atomMove);

        this.getController().addActivity(new ActivityIntegrate(integrator));

        // extend potential range, so that atoms that move outside the truncation range will still interact
        // atoms that move in will not interact since they won't be neighbors
        //XXX we don't want to do this because our potential is shifted!
        ((P2SoftSphericalTruncated) potential).setTruncationRadius(0.6 * boundary.getBoxSize().getX(0));
    }

    /**
     * @param args filename containing simulation parameters
     * @see SimEinStep2.SimOverlapParam
     */
    public static void main(String[] args) {
        //set up simulation parameters
        SimOverlapParam params = new SimOverlapParam();
        if (args.length == 0) {
            params.numMolecules = 256;
            params.rc = 2.7;
            params.slanty = false;
            params.numSteps = 1000000;
            params.temperature = 1;
            params.exponentN = 0;
            params.f = 0.01;
        }
        else {
            ParseArgs.doParseArgs(params, args);
        }

        double density = params.density;
        int exponentN = params.exponentN;
        long numSteps = params.numSteps;
        final int numMolecules = params.numMolecules;
        double temperature = params.temperature;
        double rc = params.rc;
        double spring = params.spring;
        double f = params.f;
        double x0 = params.x0;
        boolean slanty = params.slanty;

        System.out.println("Running soft sphere overlap simulation");
        System.out.println(numMolecules+" atoms at density "+density+" and temperature "+temperature);
        System.out.println("exponent N: "+ exponentN);
        System.out.println(numSteps+" steps");
        System.out.println("spring: "+spring);
        System.out.println("f: "+f);

        final long startTime = System.currentTimeMillis();

        double c = Math.exp(x0);
        double xf = Math.log(spring + c);
        double x=x0+(xf-x0)*f;
        double lambda=(Math.exp(x)-c);
        System.out.println("lambda: "+lambda);
//        lambda = c*(Math.pow(spring/c+1, f) - 1)


        //instantiate simulation
        final SimEinStep2 sim = new SimEinStep2(Space.getInstance(3), numMolecules, density, temperature, lambda, exponentN, rc, slanty);
        final double latticeEnergy = sim.potentialMaster.computeAll(false);
        System.out.println("uLat="+latticeEnergy/numMolecules);

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

        sim.initialize(numMolecules*100);
        System.out.flush();

        final MeterPotentialEnergy meterPEHarmonic = new MeterPotentialEnergy(sim.potentialMasterHarmonic);
        int numBlocks = 100;
        int interval = numMolecules;
        long blockSize = numSteps/(numBlocks*interval);
        if (blockSize == 0) blockSize = 1;
        AccumulatorAverageFixed accumulator = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPump = new DataPumpListener(meterPEHarmonic, accumulator, interval);
        sim.integrator.getEventManager().addListener(accumulatorPump);

//        final MeterPotentialEnergyFromIntegrator meterPEInt = new MeterPotentialEnergyFromIntegrator(sim.integrator);
//        if (blockSize == 0) blockSize = 1;
//        AccumulatorAverageFixed accumulatorPEInt = new AccumulatorAverageFixed(blockSize);
//        DataPumpListener accumulatorPEIntPump = new DataPumpListener(meterPEInt, accumulatorPEInt, interval);
//        sim.integrator.getEventManager().addListener(accumulatorPEIntPump)

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));

        // potentialMasterHarmonic really just gives us sum[r^2]

        DataGroup data = (DataGroup)accumulator.getData();
        double l = lambda > 0 ? lambda : 1;
        double dataErr = data.getData(accumulator.ERROR.index).getValue(0) / l;
        double dataAvg = data.getData(accumulator.AVERAGE.index).getValue(0) / l;
        double dataCorrelation = data.getData(accumulator.BLOCK_CORRELATION.index).getValue(0);
        System.out.println("msd/T  "+dataAvg/(temperature*numMolecules)+" "+dataErr/(temperature*numMolecules)+" "+dataCorrelation);

//        DataGroup dataPEInt = (DataGroup)accumulatorPEInt.getData();
//        IData dataPEIntErr = dataPEInt.getData(accumulatorPEInt.ERROR.index);
//        IData dataPEIntAvg = dataPEInt.getData(accumulatorPEInt.AVERAGE.index);
//        IData dataPEIntCorrelation = dataPEInt.getData(accumulatorPEInt.BLOCK_CORRELATION.index);
//        System.out.println("U/T  "+(dataPEIntAvg.getValue(0))/(temperature*numMolecules)+" "+dataPEIntErr.getValue(0)/(temperature*numMolecules)+" "+dataPEIntCorrelation.getValue(0));

        long endTime = System.currentTimeMillis();
        System.out.println("time: " + (endTime - startTime)/1000.0);
    }

    public void initialize(long initSteps) {
        // equilibrate off the lattice to avoid anomolous contributions
        this.getController().runActivityBlocking(new ActivityIntegrate(this.integrator, initSteps));

    }

    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules = 256;
        public double density = 1.28;
        public int exponentN = 0;
        public long numSteps = 10000000;
        public double temperature = 2;
        public double rc = 2.7;
        public double spring = 2500*240*(3.405*3.405)/120;
        public double f = 1;
        public double x0 = 3.5;
        public boolean slanty = false;
    }
}
