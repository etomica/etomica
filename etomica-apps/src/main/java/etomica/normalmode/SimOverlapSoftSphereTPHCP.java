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
import etomica.data.types.DataGroup;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.nbr.list.NeighborListManagerSlanty;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesGeneral;
import etomica.units.Degree;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;

import java.awt.*;

/**
 * 
 * The original Bennett's Overlapping Sampling Simulation
 * 	- used to check for the computation time
 * 
 * @author Tai Boon Tan
 */
public class SimOverlapSoftSphereTPHCP extends Simulation {

    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;

    public Box box;
    public Boundary boundary;
    public int[] nCells;
    public Basis basis;
    public Primitive primitive;
    public AccumulatorAverageFixed accumulator;
    public DataPumpListener accumulatorPump;
    public MeterTargetTP meter;
    protected MCMoveAtomCoupled atomMove;
    protected PotentialMaster potentialMaster;
    protected double latticeEnergy;
    protected P1ConstraintNbrHcp p1Constraint;
    public SimOverlapSoftSphereTPHCP(Space _space, int numAtoms, double density, double temperature, double[] otherTemperatures, double[] alpha, int exponent, int numAlpha, double alphaSpan, long numSteps, double rc) {
        super(_space);

        SpeciesGeneral species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        potentialMaster = new PotentialMasterList(this, rc, new NeighborListManagerSlanty.NeighborListSlantyAgentSource(rc), space);

        // TARGET
        double a = Math.pow(Math.sqrt(2) / density, 1.0 / 3.0);
        double c = Math.sqrt(8.0 / 3.0) * a;
        int nC = (int) Math.ceil(Math.pow(numAtoms / 2, 1.0 / 3.0));
        System.out.println("nC: " + nC);
        Vector[] boxDim = new Vector[3];
        boxDim[0] = Vector.of(new double[]{nC * a, 0, 0});
        boxDim[1] = Vector.of(new double[]{-nC * a * Math.cos(Degree.UNIT.toSim(60)), nC * a * Math.sin(Degree.UNIT.toSim(60)), 0});
        boxDim[2] = Vector.of(new double[]{0, 0, nC * c});
        boundary = new BoundaryDeformablePeriodic(space,boxDim);
        box = this.makeBox(boundary);
        box.setNMolecules(species, numAtoms);

        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature, box);
        atomMove = new MCMoveAtomCoupled(potentialMaster, new MeterPotentialEnergy(potentialMaster), getRandom(), space);
        atomMove.setStepSize(0.1);
        atomMove.setStepSizeMax(0.5);
        atomMove.setDoExcludeNonNeighbors(true);
        integrator.getMoveManager().addMCMove(atomMove);
//        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);


        primitive = new PrimitiveHexagonal(space, nC * a, nC * c);

        nCells = new int[]{nC, nC, nC};
        Basis basisHCP = new BasisHcp();
        basis = new BasisBigCell(space, basisHCP, nCells);

        CoordinateDefinitionLeaf coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
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

        /*
         *  1-body Potential to Constraint the atom from moving too far
         *  	away from its lattice-site
         *
         */
        p1Constraint = new P1ConstraintNbrHcp(space, a, box);
        atomMove.setConstraint(p1Constraint);

        potentialMaster.lrcMaster().setEnabled(false);

        if (potentialMaster instanceof PotentialMasterList) {
            ((PotentialMasterList) potentialMaster).setRange(rc);
            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
            ((PotentialMasterList) potentialMaster).getNeighborManager(box).reset();

        }

        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster, box);
        latticeEnergy = meterPE.getDataAsScalar();
        System.out.println("lattice energy: " + latticeEnergy / numAtoms);

        meter = new MeterTargetTP(potentialMaster, species, this, coordinateDefinition);
        meter.setLatticeEnergy(latticeEnergy);
        meter.setTemperature(temperature);
        meter.setOtherTemperatures(otherTemperatures);
        meter.setAlpha(alpha);
        meter.setAlphaSpan(alphaSpan);
        meter.setNumAlpha(numAlpha);
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

        this.getController().addActivity(new ActivityIntegrate(integrator));

        if (potentialMaster instanceof PotentialMasterList) {
            // extend potential range, so that atoms that move outside the truncation range will still interact
            // atoms that move in will not interact since they won't be neighbors
            ((P2SoftSphericalTruncated) potential).setTruncationRadius(0.6 * boundary.getBoxSize().getX(0));
        }
    }

    /**
     * @param args filename containing simulation parameters
     * @see SimOverlapSoftSphereTPHCP.SimOverlapParam
     */
    public static void main(String[] args) {
        //set up simulation parameters
        SimOverlapParam params = new SimOverlapParam();
        String inputFilename = null;
        if (args.length > 0) {
            inputFilename = args[0];
        }
        if (inputFilename != null) {
            ReadParameters readParameters = new ReadParameters(inputFilename, params);
            readParameters.readParameters();
        }
        double density = params.density;
        int exponentN = params.exponentN;
        long numSteps = params.numSteps;
        final int numMolecules = params.numMolecules;
        double temperature = params.temperature;
        double[] otherTemperatures = params.otherTemperatures;
        double[] alpha = params.alpha;
        int numAlpha = params.numAlpha;
        double alphaSpan = params.alphaSpan;
        double rc = params.rc;

        System.out.println("Running HCP soft sphere overlap simulation");
        System.out.println(numMolecules+" atoms at density "+density+" and temperature "+temperature);
        System.out.println("exponent N: "+ exponentN);
        System.out.println(numSteps+" steps");

        //instantiate simulation
        final SimOverlapSoftSphereTPHCP sim = new SimOverlapSoftSphereTPHCP(Space.getInstance(3), numMolecules, density, temperature, otherTemperatures, alpha, exponentN, numAlpha, alphaSpan, numSteps, rc);
        if (true) {
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

        sim.initialize(numMolecules*50);
        System.out.flush();

        final long startTime = System.currentTimeMillis();

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));
        //MeterTargetTP.openFW("x"+numMolecules+".dat");
        //MeterTargetTP.closeFW();

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
        System.out.println("Time taken(s): " + (endTime - startTime)/1000.0);
        System.out.println("p1 Counter: " + sim.p1Constraint.getp1Counter());
    }

    public void initialize(long initSteps) {
        // equilibrate off the lattice to avoid anomolous contributions
        System.out.println("Equilibration Steps: " + initSteps);
        this.getController().runActivityBlocking(new ActivityIntegrate(this.integrator, initSteps));


        accumulator.reset();



    }
    
    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules = 1024;
        public double density = 1.1964;
        public int exponentN = 12;
        public long numSteps = 1000000;
        public double temperature = 1.1;
        public double[] alpha = new double[]{0.011, 0.011};
        public int numAlpha = 11;
        public double alphaSpan = 1;
        public double[] otherTemperatures = new double[]{1.09,1.11};
        public double rc = 2.2;
    }
}
