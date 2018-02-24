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
import etomica.math.function.Function;
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
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;

/**
 * Simulation to run density perturbation for soft spheres.
 * 
 * @author Andrew Schultz
 */
public class SimOverlapSoftSphereDP extends Simulation {

    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public Box box;
    public Boundary boundary;
    public int[] nCells;
    public Basis basis;
    public Primitive primitive;
    public AccumulatorAverageFixed accumulator;
    public DataPumpListener accumulatorPump;
    public MeterDP meter;
    protected MCMoveAtomCoupled atomMove;
    protected PotentialMaster potentialMaster;
    protected double latticeEnergy;
    public SimOverlapSoftSphereDP(Space _space, int numAtoms, boolean slanty, final double rho, double temperature, double[] otherRho, double[] P, final int exponent, int numP, double pSpan, long numSteps, double rc) {
        super(_space);

        // rc is the cutoff at unit density
        rc *= Math.pow(rho, -1.0 / 3.0);
        if (slanty) {
            BoxAgentSourceCellManagerList boxAgentSource = new BoxAgentSourceCellManagerList(this, null, space);
            BoxAgentManager<NeighborCellManager> boxAgentManager = new BoxAgentManager<NeighborCellManager>(boxAgentSource, this);
            potentialMaster = new PotentialMasterList(this, rc, boxAgentSource, boxAgentManager, new NeighborListManagerSlanty.NeighborListSlantyAgentSource(rc), space);
        } else {
            potentialMaster = new PotentialMasterList(this, space);
        }

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        // TARGET
        double nbrDistance = 0;
        if (slanty) {
            int c = (int) Math.round(Math.pow(numAtoms, 1.0 / 3.0));
            nCells = new int[]{c, c, c};

            double L = Math.pow(Math.sqrt(2) / rho, 1.0 / 3.0);
            nbrDistance = L;
            double angle = Math.PI / 3;

//            primitive = new PrimitiveFcc(space, L*c);
            primitive = new PrimitiveTriclinic(space, L * c, L * c, L * c, angle, angle, angle);

            boundary = new BoundaryDeformablePeriodic(space, primitive.vectors());
            ((BoundaryDeformablePeriodic) boundary).setTruncationRadius(rc);
            Basis basisSimple = new Basis(new Vector3D[]{new Vector3D(0.0, 0.0, 0.0)});
            basis = new BasisBigCell(space, basisSimple, nCells);
        } else {

            double L = Math.pow(4.0 / rho, 1.0 / 3.0);
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
//        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);
        integrator.getMoveManager().addMCMove(atomMove);



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

        meter = new MeterDP(potentialMaster, species, space, this);
        meter.setCoordinateDefinition(coordinateDefinition);
        meter.setTemperature(temperature);
        meter.setOtherDensities(otherRho);
        final double uLat = meterPE.getDataAsScalar() / numAtoms;
        double uLat1 = uLat * Math.pow(rho, -exponent / 3.0);
        System.out.println("uLat1 " + uLat1);
        for (int i = 0; i < P.length; i++) {
            if (P[i] == -1) {
                // no estimate for pressure.  use pressure from harmonic approximation
                double Y = Math.pow(rho, -exponent / 3.0);
                double A = uLat1 / Y - 1.5 * Math.log(Y) + Math.log(rho);
                double Yi = Math.pow(otherRho[i], -exponent / 3.0);
                double Ai = uLat1 / Yi - 1.5 * Math.log(Yi) + Math.log(otherRho[i]);
                double v = 1.0 / rho;
                double vi = 1.0 / otherRho[i];
                P[i] = -(Ai - A) / (vi - v);
                System.out.println("initializing P[" + i + "] to " + P[i]);
            }
        }
        meter.setPCenter(P);
        meter.setLinPSpan(pSpan);
        meter.setNumP(numP);
        meter.setConstraint(p1Constraint);
        meter.setULatFunction(new Function() {
            public double f(double xRho) {
                return uLat * Math.pow(xRho / rho, exponent / 3.0);
            }
        });
        int numBlocks = 100;
        int interval = numAtoms;
        long blockSize = numSteps / (numBlocks * interval);
        if (blockSize == 0) blockSize = 1;
        System.out.println("block size " + blockSize + " interval " + interval);
        if (otherRho.length > 1) {
            accumulator = new AccumulatorAverageCovariance(blockSize);
        } else {
            accumulator = new AccumulatorAverageFixed(blockSize);
        }
        accumulatorPump = new DataPumpListener(meter, accumulator, interval);
        integrator.getEventManager().addListener(accumulatorPump);

        activityIntegrate = new ActivityIntegrate(integrator);

        getController().addAction(activityIntegrate);

        if (potentialMaster instanceof PotentialMasterList) {
            // extend potential range, so that atoms that move outside the truncation range will still interact
            // atoms that move in will not interact since they won't be neighbors
            ((P2SoftSphericalTruncated) potential).setTruncationRadius(0.6 * boundary.getBoxSize().getX(0));
        }
    }

    /**
     * @param args filename containing simulation parameters
     * @see SimOverlapSoftSphereDP.SimOverlapParam
     */
    public static void main(String[] args) {
        //set up simulation parameters
        SimOverlapParam params = new SimOverlapParam();
        ParseArgs.doParseArgs(params, args);

        double rho = params.rho;
        int exponentN = params.exponentN;
        long numSteps = params.numSteps;
        final int numMolecules = params.numMolecules;
        boolean slanty = params.slanty;
        double temperature = params.temperature;
        double[] otherRho = params.otherRho;
        double[] P = params.P;
        int numAlpha = params.numP;
        double alphaSpan = params.pSpan;
        double rc = params.rc;

        System.out.println("Running soft sphere overlap simulation");
        System.out.println(numMolecules+" atoms at density "+rho+" and temperature "+temperature);
        System.out.println("exponent N: "+ exponentN);
        System.out.println(numSteps+" steps");

        //instantiate simulation
        final SimOverlapSoftSphereDP sim = new SimOverlapSoftSphereDP(Space.getInstance(3), numMolecules, slanty, rho, temperature, otherRho, P, exponentN, numAlpha, alphaSpan, numSteps, rc);
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

        sim.initialize(numMolecules*50);
        System.out.flush();

        final long startTime = System.currentTimeMillis();

        sim.activityIntegrate.setMaxSteps(numSteps);

        //MeterTargetTP.openFW("x"+numMolecules+".dat");
        sim.getController().actionPerformed();
        //MeterTargetTP.closeFW();

        System.out.println("\nratio averages:\n");

        DataGroup data = (DataGroup)sim.accumulator.getData();
        IData dataErr = data.getData(AccumulatorAverage.ERROR.index);
        IData dataAvg = data.getData(AccumulatorAverage.AVERAGE.index);
        IData dataCorrelation = data.getData(AccumulatorAverage.BLOCK_CORRELATION.index);
        for (int i=0; i<otherRho.length; i++) {
            System.out.println(otherRho[i]);
            double[] iAlpha = sim.meter.getP(i);
            for (int j=0; j<numAlpha; j++) {
                System.out.println("  "+iAlpha[j]+" "+dataAvg.getValue(i*numAlpha+j)
                        +" "+dataErr.getValue(i*numAlpha+j)
                        +" "+dataCorrelation.getValue(i*numAlpha+j));
            }
        }

        if (otherRho.length == 2) {
            // we really kinda want each covariance for every possible pair of alphas,
            // but we're going to be interpolating anyway and the covariance is almost
            // completely insensitive to choice of alpha.  so just take the covariance for
            // the middle alphas.
            IData dataCov = data.getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE.index);
            System.out.print("covariance "+otherRho[1]+" / "+otherRho[0]+"   ");
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
    
    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules = 108;
        public boolean slanty = false;
        public double rho = 1.03;
        public int exponentN = 100;
        public long numSteps = 1000000;
        public double temperature = 1;
        public double[] P = new double[]{-1,-1};
        public int numP = 3;
        public double pSpan = 0.2;
        public double[] otherRho = new double[]{1.02,1.04};
        public double rc = 1.2;
    }
}
