package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.IMolecule;
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
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.list.BoxAgentSourceCellManagerList;
import etomica.nbr.list.NeighborListManagerSlanty;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Degree;
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
public class SimEinStep2HCP extends Simulation {

    public final PotentialMasterList potentialMaster;
    public final PotentialMasterMonatomic potentialMasterHarmonic;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public Box box;
    public BoundaryDeformableLattice boundary;
    public int[] nCells;
    public Basis basis;
    public Primitive primitive;
    public MCMoveBoxStep atomMove;
    public SimEinStep2HCP(Space _space, int numAtoms, double density, double temperature, double lambda, int exponent, double rc, double coa) {
        super(_space);

        BoxAgentSourceCellManagerList boxAgentSource = new BoxAgentSourceCellManagerList(this, null, space);
        BoxAgentManager<NeighborCellManager> boxAgentManager = new BoxAgentManager<NeighborCellManager>(boxAgentSource, NeighborCellManager.class);
        potentialMaster = new PotentialMasterList(this, rc, boxAgentSource, boxAgentManager, new NeighborListManagerSlanty.NeighborListSlantyAgentSource(rc, space), space);

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        // TARGET
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);

        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature);

        int n = (int)Math.round(Math.pow(numAtoms/8, 1.0/3.0));
        if (8*n*n*n != numAtoms) {
            throw new RuntimeException("Not compatible with HCP");
        }

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


        CoordinateDefinitionLeaf coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(nCells);

        Potential2SoftSpherical potential = null;
        if (exponent > 0) {
            potential = new P2SoftSphere(space, 1.0, 1.0, exponent);
        }
        else {
            potential = new P2LennardJones(space);
        }
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

        final MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(box);
        double latticeEnergy = meterPE.getDataAsScalar();
//        System.out.println("uLat "+latticeEnergy/numAtoms);

        P1HarmonicSite p1Harmonic = new P1HarmonicSite(space);
        p1Harmonic.setSpringConstant(1);
        p1Harmonic.setAtomAgentManager(box,coordinateDefinition.siteManager);
        potentialMasterHarmonic = new PotentialMasterMonatomic(this);
        potentialMasterHarmonic.addPotential(p1Harmonic, new AtomType[]{sphereType, sphereType});

        MeterPotentialEnergyComposite meterPEComposite = new MeterPotentialEnergyComposite(potentialMasterHarmonic, potentialMaster, latticeEnergy, numAtoms);
        meterPEComposite.setBox(box);
        meterPEComposite.setLambda(lambda);
        if (false) {
            atomMove = new MCMoveAtom(null, meterPEComposite, random, space, 0.1, 1, false);
        }
        else {
            atomMove = new MCMoveAtomCoupled(potentialMaster, lambda==0 ? meterPE : meterPEComposite, getRandom(), space);
            atomMove.setStepSize(0.1);
            atomMove.setStepSizeMax(0.5);
            ((MCMoveAtomCoupled)atomMove).setPotential(potential);
            ((MCMoveAtomCoupled)atomMove).setDoExcludeNonNeighbors(true);
        }
        integrator.getMoveManager().addMCMove(atomMove);
//        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);
//        MeterPotentialEnergyComposite meterPEComposite2 = new MeterPotentialEnergyComposite(potentialMasterHarmonic, potentialMaster, latticeEnergy, numAtoms);
//        meterPEComposite2.setBox(box);
//        meterPEComposite2.setLambda(lambda);
//        integrator.setMeterPotentialEnergy(lambda==0 ? meterPE : meterPEComposite2);


        activityIntegrate = new ActivityIntegrate(integrator);

        getController().addAction(activityIntegrate);

        // extend potential range, so that atoms that move outside the truncation range will still interact
        // atoms that move in will not interact since they won't be neighbors
        //XXX we don't want to do this because our potential is shifted!
        ((P2SoftSphericalTruncated)potential).setTruncationRadius(0.6*boundary.getBoxSize().getX(0));
    }

    /**
     * @param args filename containing simulation parameters
     * @see SimEinStep2HCP.SimOverlapParam
     */
    public static void main(String[] args) {
        //set up simulation parameters
        SimOverlapParam params = new SimOverlapParam();
        if (args.length == 0) {
            params.numMolecules = 512;
            params.rc = 3;
            params.numSteps = 10000000;
            params.density = 1.4;
            params.temperature = 0.3;
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
        double rc = params.rc*Math.pow(density, -1.0/3.0);
        double spring = params.spring;
        double f = params.f;
        double x0 = params.x0;

        System.out.println("Running soft sphere overlap simulation");
        System.out.println(numMolecules+" atoms at density "+density+" and temperature "+temperature);
        System.out.println("exponent N: "+ exponentN);
        System.out.println("rc: "+rc);
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
        final SimEinStep2HCP sim = new SimEinStep2HCP(Space.getInstance(3), numMolecules, density, temperature, lambda, exponentN, rc, Math.sqrt(8.0/3.0));
        final MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster);
        meterPE.setBox(sim.box);
        final double latticeEnergy = meterPE.getDataAsScalar();
        System.out.println("uLat="+latticeEnergy/numMolecules);

        if (false) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, sim.space, sim.getController());
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

        sim.activityIntegrate.setMaxSteps(numSteps);

        // potentialMasterHarmonic really just gives us sum[r^2]
        final MeterPotentialEnergy meterPEHarmonic = new MeterPotentialEnergy(sim.potentialMasterHarmonic);
        meterPEHarmonic.setBox(sim.box);
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
//        sim.integrator.getEventManager().addListener(accumulatorPEIntPump);


        sim.getController().actionPerformed();

        DataGroup data = (DataGroup)accumulator.getData();
        IData dataErr = data.getData(AccumulatorAverage.ERROR.index);
        IData dataAvg = data.getData(AccumulatorAverage.AVERAGE.index);
        IData dataCorrelation = data.getData(AccumulatorAverage.BLOCK_CORRELATION.index);
        System.out.println("msd/T  "+dataAvg.getValue(0)/(temperature*numMolecules)+" "+dataErr.getValue(0)/(temperature*numMolecules)+" "+dataCorrelation.getValue(0));

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
        activityIntegrate.setMaxSteps(initSteps);
        getController().actionPerformed();
        getController().reset();
    }
    
    protected static class MeterPotentialEnergyComposite extends
            MeterPotentialEnergy {
        protected final MeterPotentialEnergy meterPE1, meterPE2;
        protected final int nMolecules;
        protected double lambda, latticeEnergy;
        protected boolean hasTarget = false;

        protected MeterPotentialEnergyComposite(PotentialMaster potentialMaster1, PotentialMaster potentialMaster2, double latticeEnergy, int nMolecules) {
            super(null);
            meterPE1 = new MeterPotentialEnergy(potentialMaster1);
            meterPE2 = new MeterPotentialEnergy(potentialMaster2);
            this.latticeEnergy = latticeEnergy;
            this.nMolecules = nMolecules;
        }

        public double getLambda() {
            return lambda;
        }

        public void setLambda(double newFrac) {
            lambda = newFrac;
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
            hasTarget = atom != null;
        }

        public void setTarget(IMolecule mole) {
            meterPE1.setTarget(mole);
            meterPE2.setTarget(mole);
            hasTarget = mole != null;
        }

        public double getDataAsScalar() {

//            System.out.println(hasTarget+" "+meterPE1.getDataAsScalar()+" "+meterPE2.getDataAsScalar()+" "+(latticeEnergy/(hasTarget ? nMolecules/2 : 1)));
//            System.out.println("hi "+(lambda * meterPE1.getDataAsScalar() + (meterPE2.getDataAsScalar() - latticeEnergy/(hasTarget ? nMolecules/2 : 1))));
            return lambda * meterPE1.getDataAsScalar() + (meterPE2.getDataAsScalar() - latticeEnergy/(hasTarget ? nMolecules/2 : 1));
        }
    }

    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules = 512;
        public double density = 1.28;
        public int exponentN = 0;
        public long numSteps = 10000000;
        public double temperature = 2;
        public double rc = 2.7;
        public double spring = 50000;
        public double f = 1;
        public double x0 = 3.5;
    }
}
