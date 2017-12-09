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
import etomica.data.types.DataGroup;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.*;
import etomica.molecule.IMolecule;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.list.BoxAgentSourceCellManagerList;
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
import etomica.units.dimensions.Energy;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.ReadParameters;

import java.awt.*;

/**
 * Simulation that samples a composite energy function (soft sphere and
 * Einstein crystal) and perturbs into overlap regions shared by systems with
 * more or less soft sphere contributions.
 * 
 * @author Andrew Schultz
 */
public class SimOverlapSoftSphereEin extends Simulation {

    private static final long serialVersionUID = 1L;
    public final PotentialMasterList potentialMaster;
    public final PotentialMasterMonatomic potentialMasterHarmonic;
    public final double latticeEnergy;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public Box box;
    public Boundary boundary;
    public int[] nCells;
    public Basis basis;
    public Primitive primitive;
    public AccumulatorAverageFixed accumulator;
    public DataPumpListener accumulatorPump;
    public MeterOverlapSwitch meter;
    public MCMoveAtomCoupled atomMove;
    public SimOverlapSoftSphereEin(Space _space, int numAtoms, double density, boolean slanty, double temperature, double spring, double frac, double[] otherFrac, double[] alpha, int exponent, int numAlpha, double alphaSpan, long numSteps, double rc) {
        super(_space);

        if (slanty) {
            BoxAgentSourceCellManagerList boxAgentSource = new BoxAgentSourceCellManagerList(this, null, space);
            BoxAgentManager<NeighborCellManager> boxAgentManager = new BoxAgentManager<NeighborCellManager>(boxAgentSource, NeighborCellManager.class, this);
            potentialMaster = new PotentialMasterList(this, rc, boxAgentSource, boxAgentManager, new NeighborListManagerSlanty.NeighborListSlantyAgentSource(rc, space), space);
        }
        else {
            potentialMaster = new PotentialMasterList(this, space);
        }

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        // TARGET
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);

        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature);
        double nbrDistance = 0;
        if (slanty) {
            int c = (int)Math.round(Math.pow(numAtoms, 1.0/3.0));
            nCells = new int[]{c,c,c};

            double L = Math.pow(Math.sqrt(2)/density, 1.0/3.0);
            nbrDistance = L;
            double angle = Math.PI/3;

//            primitive = new PrimitiveFcc(space, L*c);
            primitive = new PrimitiveTriclinic(space, L*c,L*c,L*c, angle,angle,angle);

            boundary = new BoundaryDeformablePeriodic(space, primitive.vectors());
            ((BoundaryDeformablePeriodic)boundary).setTruncationRadius(rc);
            Basis basisSimple = new Basis(new Vector3D[]{new Vector3D(0.0, 0.0, 0.0)});
            basis = new BasisBigCell(space, basisSimple, nCells);
        }
        else {

            double L = Math.pow(4.0/density, 1.0/3.0);
            nbrDistance = L / Math.sqrt(2);
            int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
            primitive = new PrimitiveCubic(space, n*L);

            nCells = new int[]{n,n,n};
            boundary = new BoundaryRectangularPeriodic(space, n * L);
            Basis basisFCC = new BasisCubicFcc();
            basis = new BasisBigCell(space, basisFCC, nCells);
        }
        System.out.println("nbr distance "+nbrDistance);

        box.setBoundary(boundary);

        CoordinateDefinitionLeaf coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1,1,1});

        Potential2SoftSpherical potential = new P2SoftSphere(space, 1.0, 1.0, exponent);
        potential = new P2SoftSphericalTruncated(space, potential, rc);
        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(potential, new AtomType[]{sphereType, sphereType});


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

        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(box);
        latticeEnergy = meterPE.getDataAsScalar();
        System.out.println("uLat "+latticeEnergy/numAtoms);

        P1HarmonicSite p1Harmonic = new P1HarmonicSite(space);
        p1Harmonic.setSpringConstant(spring);
        p1Harmonic.setAtomAgentManager(box,coordinateDefinition.siteManager);
        potentialMasterHarmonic = new PotentialMasterMonatomic(this);
        potentialMasterHarmonic.addPotential(p1Harmonic, new AtomType[]{sphereType, sphereType});

        MeterPotentialEnergyComposite meterPEComposite = new MeterPotentialEnergyComposite(potentialMasterHarmonic, potentialMaster, latticeEnergy);
        meterPEComposite.setBox(box);
        meterPEComposite.setFrac(frac);
        atomMove = new MCMoveAtomCoupled(potentialMaster, meterPEComposite, getRandom(), space);
        atomMove.setStepSize(0.1);
        atomMove.setStepSizeMax(0.5);
        atomMove.setDoExcludeNonNeighbors(true);
        atomMove.setPotential(potential);
        P1ConstraintNbr p1Constraint = new P1ConstraintNbr(space, nbrDistance, this);
//        atomMove.setConstraint(p1Constraint);
        integrator.getMoveManager().addMCMove(atomMove);
//      ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);


        meter = new MeterOverlapSwitch(potentialMasterHarmonic, potentialMaster);
        meter.setBox(box);
        meter.setLatticeEnergy(latticeEnergy);
        meter.setTemperature(temperature);
        meter.setSampledSwitchFrac(frac);
        meter.setOtherSwitchFrac(otherFrac);
        meter.setAlpha(alpha);
        meter.setAlphaSpan(alphaSpan);
        meter.setNumAlpha(numAlpha);
        int numBlocks = 100;
        int interval = numAtoms;
        long blockSize = numSteps/(numBlocks*interval);
        if (blockSize == 0) blockSize = 1;
        System.out.println("block size "+blockSize+" interval "+interval);
        if (otherFrac.length > 1) {
            accumulator = new AccumulatorAverageCovariance(blockSize);
        }
        else {
            accumulator = new AccumulatorAverageFixed(blockSize);
        }
        accumulatorPump = new DataPumpListener(meter, accumulator, interval);
        integrator.getEventManager().addListener(accumulatorPump);

        activityIntegrate = new ActivityIntegrate(integrator);

        getController().addAction(activityIntegrate);

        // extend potential range, so that atoms that move outside the truncation range will still interact
        // atoms that move in will not interact since they won't be neighbors
        ((P2SoftSphericalTruncated)potential).setTruncationRadius(0.6*boundary.getBoxSize().getX(0));
    }

    /**
     * @param args filename containing simulation parameters
     * @see SimOverlapSoftSphereEin.SimOverlapParam
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
        if (args.length > 1) {
            // we want to skip the first arg
            String[] otherArgs = new String[args.length-1];
            System.arraycopy(args, 1, otherArgs, 0, otherArgs.length);
            ParseArgs parseArgs = new ParseArgs(params);
            parseArgs.parseArgs(otherArgs);
        }
        double density = params.density;
        boolean slanty = params.slanty;
        int exponentN = params.exponentN;
        long numSteps = params.numSteps;
        final int numMolecules = params.numMolecules;
        double temperature = params.temperature;
        double frac = params.frac;
        double[] otherFrac = params.otherFrac;
        double[] alpha = params.alpha;
        int numAlpha = params.numAlpha;
        double alphaSpan = params.alphaSpan;
        double rc = params.rc;
        double spring = params.spring;

        System.out.println("Running soft sphere overlap simulation");
        System.out.println(numMolecules+" atoms at density "+density+" and temperature "+temperature);
        System.out.println("exponent N: "+ exponentN);
        System.out.println(numSteps+" steps");

        //instantiate simulation
        final SimOverlapSoftSphereEin sim = new SimOverlapSoftSphereEin(Space.getInstance(3), numMolecules, density, slanty, temperature, spring, frac, otherFrac, alpha, exponentN, numAlpha, alphaSpan, numSteps, rc);
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
        sim.meter.count = 0;
        sim.meter.targetSum = 0;
        sim.meter.refSum = 0;
        System.out.flush();

        final long startTime = System.currentTimeMillis();

        sim.activityIntegrate.setMaxSteps(numSteps);

        final MeterPotentialEnergy meterPEHarmonic = new MeterPotentialEnergy(sim.potentialMasterHarmonic);
        meterPEHarmonic.setBox(sim.box);
        final MeterPotentialEnergy meterPETarget = new MeterPotentialEnergy(sim.potentialMaster);
        meterPETarget.setBox(sim.box);
        meterPETarget.setIncludeLrc(false);
        DataSourceScalar meterPEdiff = new DataSourceScalar("PE diff", Energy.DIMENSION) {
            public double getDataAsScalar() {
                return meterPETarget.getDataAsScalar() - sim.latticeEnergy - meterPEHarmonic.getDataAsScalar();
            }
        };
        int numBlocks = 100;
        int interval = numMolecules;
        long blockSize = numSteps/(numBlocks*interval);
        if (blockSize == 0) blockSize = 1;
        AccumulatorAverageFixed accumulator = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPump = new DataPumpListener(meterPEdiff, accumulator, interval);
        sim.integrator.getEventManager().addListener(accumulatorPump);

        //MeterTargetTP.openFW("x"+numMolecules+".dat");
        sim.getController().actionPerformed();
        //MeterTargetTP.closeFW();

        System.out.println("average delta U " + accumulator.getData().getValue(AccumulatorAverage.AVERAGE.index) + " " + accumulator.getData().getValue(AccumulatorAverage.ERROR.index) + " " + accumulator.getData().getValue(AccumulatorAverage.BLOCK_CORRELATION.index));

        System.out.println("\nratio averages:\n");

        DataGroup data = (DataGroup)sim.accumulator.getData();
        IData dataErr = data.getData(AccumulatorAverage.ERROR.index);
        IData dataAvg = data.getData(AccumulatorAverage.AVERAGE.index);
        IData dataCorrelation = data.getData(AccumulatorAverage.BLOCK_CORRELATION.index);
        for (int i=0; i<otherFrac.length; i++) {
            System.out.println(otherFrac[i]);
            double[] iAlpha = sim.meter.getAlpha(i);
            for (int j=0; j<numAlpha; j++) {
                System.out.println("  "+iAlpha[j]+" "+dataAvg.getValue(i*numAlpha+j)
                        +" "+dataErr.getValue(i*numAlpha+j)
                        +" "+dataCorrelation.getValue(i*numAlpha+j));
            }
        }

        if (otherFrac.length == 2) {
            // we really kinda want each covariance for every possible pair of alphas,
            // but we're going to be interpolating anyway and the covariance is almost
            // completely insensitive to choice of alpha.  so just take the covariance for
            // the middle alphas.
            IData dataCov = data.getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE.index);
            System.out.print("covariance "+otherFrac[1]+" / "+otherFrac[0]+"   ");
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
        activityIntegrate.setMaxSteps(initSteps);
        getController().actionPerformed();
        getController().reset();

        accumulator.reset();

        getController().reset();

    }
    
    protected static class MeterPotentialEnergyComposite extends
            MeterPotentialEnergy {
        protected final MeterPotentialEnergy meterPE1, meterPE2;
        protected double frac, latticeEnergy;
        
        protected MeterPotentialEnergyComposite(PotentialMaster potentialMaster1, PotentialMaster potentialMaster2, double latticeEnergy) {
            super(null);
            meterPE1 = new MeterPotentialEnergy(potentialMaster1);
            meterPE2 = new MeterPotentialEnergy(potentialMaster2);
            this.latticeEnergy = latticeEnergy;
        }
        
        public double getFrac() {
            return frac;
        }

        public void setFrac(double newFrac) {
            frac = newFrac;
        }

        public Box getBox() {
            return meterPE1.getBox();
        }

        public void setBox(Box box) {
            meterPE1.setBox(box);
            meterPE2.setBox(box);
        }

        public boolean isIncludeLrc() {
            return meterPE1.isIncludeLrc();
        }

        public void setIncludeLrc(boolean b) {
            meterPE1.setIncludeLrc(b);
            meterPE2.setIncludeLrc(b);
        }

        public void setTarget(IAtom atom) {
            meterPE1.setTarget(atom);
            meterPE2.setTarget(atom);
        }

        public void setTarget(IMolecule mole) {
            meterPE1.setTarget(mole);
            meterPE2.setTarget(mole);
        }

        public double getDataAsScalar() {
            return (1-frac) * meterPE1.getDataAsScalar() + frac * (meterPE2.getDataAsScalar() - latticeEnergy);
        }
    }

    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules = 216;
        public boolean slanty = true;
        public double density = 1.1964;
        public int exponentN = 12;
        public long numSteps = 1000000;
        public double temperature = 0.01;
        public double[] alpha = new double[]{1};
        public int numAlpha = 11;
        public double alphaSpan = 10;
        public double frac = 0.1;
        public double[] otherFrac = new double[]{0};
        public double rc = 2.2;
        public double spring = 215;
    }
}
