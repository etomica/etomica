/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.data.AccumulatorRatioAverageCovariance;
import etomica.data.DataPump;
import etomica.data.IData;
import etomica.data.IDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorScheme;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.*;
import etomica.listener.IntegratorListenerAction;
import etomica.molecule.IMolecule;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.list.BoxAgentSourceCellManagerList;
import etomica.nbr.list.NeighborListManagerSlanty;
import etomica.nbr.list.PotentialMasterList;
import etomica.overlap.IntegratorOverlap;
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
import etomica.util.ReadParameters;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;

import java.awt.*;

/**
 * Overlap sampling simulation whose target system has a composite energy
 * function (soft sphere and Einstein crystal) and whose reference system is
 * defined by (harmonically coupled) normal modes and Einstein crystal.
 * 
 * @author Andrew Schultz
 */
public class SimOverlapSoftSphereEinHarm extends Simulation {

    private static final long serialVersionUID = 1L;
    public final PotentialMasterList potentialMaster;
    public final PotentialMasterMonatomic potentialMasterHarmonic;
    public final double latticeEnergy;
    public final NormalModes normalModes;
    public final IntegratorOverlap integratorOverlap;
    public IntegratorBox[] integrators;
    public DataSourceVirialOverlap dsvo;
    public ActivityIntegrate activityIntegrate;
    public Box box, boxRef;
    public Boundary boundary, boundaryRef;
    public int[] nCells;
    public Basis basis;
    public Primitive primitive;
    public MCMoveAtomCoupled atomMove;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public DataPump[] accumulatorPumps;
    public IDataSource[] meters;
    public double alphaCenter, alphaSpan;
    public int numAlpha;
    public SimOverlapSoftSphereEinHarm(Space _space, int numAtoms, double density, boolean slanty, double temperature, double spring, double frac, int exponent, double rc) {
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

        integrators = new IntegratorBox[2];
        accumulatorPumps = new DataPump[2];
        meters = new IDataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];

        integrators[1] = new IntegratorMC(potentialMaster, getRandom(), temperature);
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
            boundaryRef = new BoundaryDeformablePeriodic(space, primitive.vectors());
            ((BoundaryDeformablePeriodic)boundaryRef).setTruncationRadius(rc);
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
            boundaryRef = new BoundaryRectangularPeriodic(space, n * L);
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

        integrators[1].setBox(box);

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
        ((IntegratorMC)integrators[1]).getMoveManager().addMCMove(atomMove);
//      ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);


        // HARMONIC
        boxRef = new Box(boundaryRef, space);
        addBox(boxRef);
        boxRef.setNMolecules(species, numAtoms);

        IntegratorMC integratorRef = new IntegratorMC(null, random, 1.0); //null changed on 11/20/2009

        MCMoveHarmonicEin move = new MCMoveHarmonicEin(getRandom());
        move.setAlphaEin(spring);
        integratorRef.getMoveManager().addMCMove(move);
        integrators[0] = integratorRef;

        CoordinateDefinitionLeaf coordinateDefinitionRef = new CoordinateDefinitionLeaf(boxRef, primitive, basis, space);
        coordinateDefinitionRef.initializeCoordinates(new int[]{1,1,1});

        potentialMaster.getNeighborManager(boxRef).reset();

        String inFile = "inputSSDB"+numAtoms;
        normalModes = new NormalModesFromFile(inFile, space.D());
        /*
         * nuke this line if it is DB
         */
        //normalModes.setTemperature(temperature);

        WaveVectorFactory waveVectorFactory = normalModes.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(boxRef);
        move.setOmegaSquared(normalModes.getOmegaSquared());
        move.setEigenVectors(normalModes.getEigenvectors());
        move.setWaveVectors(waveVectorFactory.getWaveVectors());
        move.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
        move.setCoordinateDefinition(coordinateDefinitionRef);
        move.setTemperature(temperature);
        move.setFrac(frac);

        move.setBox(boxRef);

        integratorRef.setBox(boxRef);

        // OVERLAP
        integratorOverlap = new IntegratorOverlap(new IntegratorBox[]{integratorRef, integrators[1]});
        MeterHarmonicEnergy meterHarmonicEnergy = new MeterHarmonicEnergy(coordinateDefinition, normalModes);
        MeterBoltzmannTarget meterTarget = new MeterBoltzmannTarget(meterPE, meterHarmonicEnergy);
        meterTarget.setTemperature(temperature);
        meterTarget.setFrac(frac);
        meterTarget.setLatticeEnergy(latticeEnergy);
        meters[1] = meterTarget;

        MeterBoltzmannHarmonic meterHarmonic = new MeterBoltzmannHarmonic(move, potentialMaster);
        meterHarmonic.setFrac(frac);
        meterHarmonic.setTemperature(temperature);
        meterHarmonic.setLatticeEnergy(latticeEnergy);
        meters[0] = meterHarmonic;

        integratorOverlap.setRefStepFraction(0.5);
        integratorOverlap.setAdjustStepFraction(false);

        activityIntegrate = new ActivityIntegrate(integratorOverlap);

        getController().addAction(activityIntegrate);

        // extend potential range, so that atoms that move outside the truncation range will still interact
        // atoms that move in will not interact since they won't be neighbors
        ((P2SoftSphericalTruncated)potential).setTruncationRadius(0.6*boundary.getBoxSize().getX(0));
    }

    /**
     * @param args filename containing simulation parameters
     * @see SimOverlapSoftSphereEinHarm.SimOverlapParam
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
        int numAlpha = params.numAlpha;
        double alphaSpan = params.alphaSpan;
        double alphaCenter = params.alphaCenter;
        double rc = params.rc;
        double spring = params.spring;

        System.out.println("Running soft sphere overlap simulation");
        System.out.println(numMolecules+" atoms at density "+density+" and temperature "+temperature);
        System.out.println("exponent N: "+ exponentN);
        System.out.println(numSteps+" steps");

        //instantiate simulation
        final SimOverlapSoftSphereEinHarm sim = new SimOverlapSoftSphereEinHarm(Space.getInstance(3), numMolecules, density, slanty, temperature, spring, frac, exponentN, rc);
        sim.setNumAlpha(numAlpha);
        sim.setAlphaCenter(alphaCenter, alphaSpan);
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

            simGraphic.makeAndDisplayFrame("SS");
            return;
        }

        //start simulation

        System.out.flush();
        long steps = numMolecules/10;
        if (steps < numSteps/20) {
            steps = numSteps/20;
        }
        sim.equilibrate(numMolecules/5);

        final long startTime = System.currentTimeMillis();

        sim.activityIntegrate.setMaxSteps(numSteps);

        //MeterTargetTP.openFW("x"+numMolecules+".dat");
        sim.getController().actionPerformed();
        //MeterTargetTP.closeFW();

        System.out.println("\nratio averages:\n");
        double[] lnAlphaDiff = new double[numAlpha];
        double[] lnAlpha = new double[numAlpha];
        double[] err = new double[numAlpha];

        for (int i=0; i<2; i++) {
            System.out.println(i==0 ? "reference" : "target");
//            IData dataCorrelation = data.getData(AccumulatorRatioAverageCovariance.StatType.BLOCK_CORRELATION.index);
            for (int j=0; j<numAlpha; j++) {
                double jAlpha = sim.accumulators[i].getBennetBias(j);
                DataGroup data = (DataGroup)sim.accumulators[i].getData(j);
                IData dataErr = data.getData(AccumulatorRatioAverageCovariance.RATIO_ERROR.index);
                IData dataAvg = data.getData(AccumulatorRatioAverageCovariance.RATIO.index);
                System.out.println("  "+jAlpha+" "+dataAvg.getValue(1)
                        +" "+dataErr.getValue(1));
                lnAlphaDiff[j] += (1-i*2)*Math.log(dataAvg.getValue(1));
                double errRatio = dataErr.getValue(1) / dataAvg.getValue(1);
                err[j] += errRatio*errRatio;
            }
            System.out.println();
        }
        System.out.println("all alpha");
        for (int j=0; j<numAlpha; j++) {
            System.out.println(sim.accumulators[0].getBennetBias(j)+" "+sim.dsvo.getAverage(j)+" "+sim.dsvo.getError(j));
        }
        double[] avgerr = sim.dsvo.getOverlapAverageAndError();
        System.out.println("\nnew alpha "+avgerr[0]+" "+avgerr[1]);

        System.out.println("delta A "+(-temperature*Math.log(avgerr[0]))+" "+temperature*avgerr[1]/avgerr[0]);

        long endTime = System.currentTimeMillis();
        System.out.println("Time taken: " + (endTime - startTime)/1000.0);
    }

    public void setAlphaCenter(double newAlphaCenter, double span) {
        alphaCenter = newAlphaCenter;
        alphaSpan = span;
        if (accumulators[0] != null) {
            accumulators[0].setBennetParam(alphaCenter, span);
            accumulators[1].setBennetParam(alphaCenter, span);
        }
    }

    public void setAccumulator(AccumulatorVirialOverlapSingleAverage newAccumulator, int iBox) {

        accumulators[iBox] = newAccumulator;

        newAccumulator.setBlockSize(200); // setting the block size = 300

        if (accumulatorPumps[iBox] == null) {
            accumulatorPumps[iBox] = new DataPump(meters[iBox], newAccumulator);
            IntegratorListenerAction pumpListener = new IntegratorListenerAction(accumulatorPumps[iBox]);
            integrators[iBox].getEventManager().addListener(pumpListener);
            pumpListener.setInterval(box.getLeafList().getAtomCount());
        } else {
            accumulatorPumps[iBox].setDataSink(newAccumulator);
        }
        if (integratorOverlap != null && accumulators[0] != null && accumulators[1] != null) {
            dsvo = new DataSourceVirialOverlap(accumulators[0], accumulators[1]);
            integratorOverlap.setReferenceFracSource(dsvo);
        }
    }

    public void setNumAlpha(int newNumAlpha) {
        numAlpha = newNumAlpha;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, numAlpha, false), 1);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, numAlpha, true), 0);
        if (alphaSpan > 0 && alphaCenter > 0) {
            setAlphaCenter(alphaCenter, alphaSpan);
        }
    }

    public void equilibrate(long initSteps) {
        // run a short simulation to get reasonable MC Move step sizes and
        // (if needed) narrow in on a reference preference
        activityIntegrate.setMaxSteps(initSteps);

        for (int i = 0; i < 2; i++) {
            if (integrators[i] instanceof IntegratorMC)
                ((IntegratorMC) integrators[i]).getMoveManager().setEquilibrating(true);
        }
        getController().actionPerformed();
        getController().reset();
        for (int i = 0; i < 2; i++) {
            if (integrators[i] instanceof IntegratorMC)
                ((IntegratorMC) integrators[i]).getMoveManager().setEquilibrating(false);
        }

        dsvo.reset();
    }
    
    protected static class MeterPotentialEnergyComposite extends MeterPotentialEnergy {
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
        public int numMolecules = 108;
        public boolean slanty = false;
        public double density = 1.1964;
        public int exponentN = 12;
        public long numSteps = 100;
        public double temperature = 0.01;
        public int numAlpha = 3;
        public double alphaSpan = .5;
        public double alphaCenter = 1;
        public double frac = 0.01;
        public double rc = 2.2;
        public double spring = 200.195510966745843;
    }
}
