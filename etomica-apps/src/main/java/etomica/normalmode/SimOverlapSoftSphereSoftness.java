/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.IDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
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
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;
import etomica.virial.overlap.IntegratorOverlap;

/**
 * Conduct Bennett's Overlap Sampling to quantify the free-energy difference
 * between soft-sphere model with different softness
 * 
 * Soft-Sphere potential: U = epsilon*(sigma/r)^n
 * 
 * where n is the exponent, s = 1/n (where s is softness)
 * epsilon and sigma is set to 1 
 *  
 *   
 * @author Tai Boon Tan
 */
public class SimOverlapSoftSphereSoftness extends Simulation {

    private static final long serialVersionUID = 1L;
    public IntegratorOverlap integratorOverlap;
    public DataSourceVirialOverlap dsvo;
    public IntegratorBox[] integrators;
    public ActivityIntegrate activityIntegrate;
    public Box boxTarg, boxRef;
    public Boundary boundaryTarg, boundaryRef;
    public int[] nCells;
    public Basis basis;
    public Primitive primitive, primitiveUnitCell;
    public double refPref;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public DataPump[] accumulatorPumps;
    public IDataSource[] meters;
    public String fname;
    protected MCMoveAtomCoupled atomMoveTarg, atomMoveRef;
    protected PotentialMaster potentialMasterTarg, potentialMasterRef;
    protected double latticeEnergyTarg, latticeEnergyRef;

    public SimOverlapSoftSphereSoftness(Space _space, int numAtoms, double density, double temperature,
                                        double harmonicFudge, int[] exponent, double alpha, double alphaSpan, int numAlpha) {
        super(_space);

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        potentialMasterTarg = new PotentialMasterList(this, space);
        potentialMasterRef = new PotentialMasterList(this, space);

        integrators = new IntegratorBox[2];
        accumulatorPumps = new DataPump[2];
        meters = new IDataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];

        // TARGET
        double L = Math.pow(4.0 / density, 1.0 / 3.0);
        int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
        boundaryTarg = new BoundaryRectangularPeriodic(space, n * L);
        boxTarg = this.makeBox(boundaryTarg);
        boxTarg.setNMolecules(species, numAtoms);

        IntegratorMC integratorTarg = new IntegratorMC(potentialMasterTarg, getRandom(), temperature, boxTarg);
        atomMoveTarg = new MCMoveAtomCoupled(potentialMasterTarg, new MeterPotentialEnergy(potentialMasterTarg), getRandom(), space);
        atomMoveTarg.setStepSize(0.1);
        atomMoveTarg.setStepSizeMax(0.5);
        atomMoveTarg.setDoExcludeNonNeighbors(true);
        integratorTarg.getMoveManager().addMCMove(atomMoveTarg);
        ((MCMoveStepTracker) atomMoveTarg.getTracker()).setNoisyAdjustment(false);

        integrators[1] = integratorTarg;

        primitive = new PrimitiveCubic(space, n * L);
        primitiveUnitCell = new PrimitiveCubic(space, L);

        nCells = new int[]{n, n, n};
        Basis basisFCC = new BasisCubicFcc();
        basis = new BasisBigCell(space, basisFCC, nCells);


        CoordinateDefinitionLeaf coordinateDefinitionTarg = new CoordinateDefinitionLeaf(boxTarg, primitive, basis, space);
        coordinateDefinitionTarg.initializeCoordinates(new int[]{1, 1, 1});

        Potential2SoftSpherical potentialTarg = new P2SoftSphere(space, 1.0, 1.0, exponent[1]);
        double truncationRadius = 2.2;//boundaryTarg.getBoxSize().getX(0) * 0.495;

        System.out.println("Truncation Radius: " + truncationRadius);

        if (potentialMasterTarg instanceof PotentialMasterList) {
            potentialTarg = new P2SoftSphericalTruncated(space, potentialTarg, truncationRadius);

        } else {
            potentialTarg = new P2SoftSphericalTruncatedShifted(space, potentialTarg, truncationRadius);

        }
        atomMoveTarg.setPotential(potentialTarg);
        AtomType sphereType = species.getLeafType();
        potentialMasterTarg.addPotential(potentialTarg, new AtomType[]{sphereType, sphereType});

        /*
         *  1-body Potential to Constraint the atom from moving too far
         *  	away from its lattice-site
         *
         */

        P1ConstraintNbr p1ConstraintTarg = new P1ConstraintNbr(space, L / Math.sqrt(2.0));
        atomMoveTarg.setConstraint(p1ConstraintTarg);
        potentialMasterTarg.lrcMaster().setEnabled(false);

        if (potentialMasterTarg instanceof PotentialMasterList) {
            double neighborRange = truncationRadius;
            int cellRange = 7;
            ((PotentialMasterList) potentialMasterTarg).setRange(neighborRange);
            ((PotentialMasterList) potentialMasterTarg).setCellRange(cellRange); // insanely high, this lets us have neighborRange close to dimensions/2
            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
            ((PotentialMasterList) potentialMasterTarg).getNeighborManager(boxTarg).reset();
            int potentialCells = ((PotentialMasterList) potentialMasterTarg).getNbrCellManager(boxTarg).getLattice().getSize()[0];
            if (potentialCells < cellRange * 2 + 1) {
                throw new RuntimeException("oops (" + potentialCells + " < " + (cellRange * 2 + 1) + ")");
            }
            if (potentialCells > cellRange * 2 + 1) {
                System.out.println("could probably use a larger truncation radius (" + potentialCells + " > " + (cellRange * 2 + 1) + ")");
            }
            ((P2SoftSphericalTruncated) potentialTarg).setTruncationRadius(0.6 * boundaryTarg.getBoxSize().getX(0));
        }


        // Reference System
        boundaryRef = new BoundaryRectangularPeriodic(space);
        boxRef = this.makeBox(boundaryRef);
        boxRef.setNMolecules(species, numAtoms);

        IntegratorMC integratorRef = new IntegratorMC(potentialMasterRef, getRandom(), temperature, boxRef);
        atomMoveRef = new MCMoveAtomCoupled(potentialMasterRef, new MeterPotentialEnergy(potentialMasterRef), getRandom(), space);
        atomMoveRef.setStepSize(0.1);
        atomMoveRef.setStepSizeMax(0.5);
        atomMoveRef.setDoExcludeNonNeighbors(true);
        integratorRef.getMoveManager().addMCMove(atomMoveRef);
        ((MCMoveStepTracker) atomMoveRef.getTracker()).setNoisyAdjustment(false);

        integrators[0] = integratorRef;

        CoordinateDefinitionLeaf coordinateDefinitionRef = new CoordinateDefinitionLeaf(boxRef, primitive, basis, space);
        coordinateDefinitionRef.initializeCoordinates(new int[]{1, 1, 1});

        Potential2SoftSpherical potentialRef = new P2SoftSphere(space, 1.0, 1.0, exponent[0]);
        if (potentialMasterRef instanceof PotentialMasterList) {
            potentialRef = new P2SoftSphericalTruncated(space, potentialRef, truncationRadius);

        } else {
            potentialRef = new P2SoftSphericalTruncatedShifted(space, potentialRef, truncationRadius);

        }
        atomMoveRef.setPotential(potentialRef);
        potentialMasterRef.addPotential(potentialRef, new AtomType[]{sphereType, sphereType});

        /*
         *  1-body Potential to Constraint the atom from moving too far
         *  	away from its lattice-site
         *
         */

        P1ConstraintNbr p1ConstraintRef = new P1ConstraintNbr(space, L / Math.sqrt(2.0));
        //potentialMasterRef.addPotential(p1ConstraintRef, new IAtomType[] {sphereType});
        atomMoveRef.setConstraint(p1ConstraintRef);
        potentialMasterRef.lrcMaster().setEnabled(false);

        if (potentialMasterRef instanceof PotentialMasterList) {
            double neighborRange = truncationRadius;
            int cellRange = 7;
            ((PotentialMasterList) potentialMasterRef).setRange(neighborRange);
            ((PotentialMasterList) potentialMasterRef).setCellRange(cellRange); // insanely high, this lets us have neighborRange close to dimensions/2
            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
            ((PotentialMasterList) potentialMasterRef).getNeighborManager(boxRef).reset();
            int potentialCells = ((PotentialMasterList) potentialMasterRef).getNbrCellManager(boxRef).getLattice().getSize()[0];
            if (potentialCells < cellRange * 2 + 1) {
                throw new RuntimeException("oops (" + potentialCells + " < " + (cellRange * 2 + 1) + ")");
            }
            if (potentialCells > cellRange * 2 + 1) {
                System.out.println("could probably use a larger truncation radius (" + potentialCells + " > " + (cellRange * 2 + 1) + ")");
            }
            ((P2SoftSphericalTruncated) potentialRef).setTruncationRadius(0.6 * boundaryRef.getBoxSize().getX(0));
        }

        MeterPotentialEnergy meterPETarg = new MeterPotentialEnergy(potentialMasterRef, boxTarg);
        latticeEnergyTarg = meterPETarg.getDataAsScalar();
        System.out.println("lattice energy/N (targ n=" + exponent[1] + "): " + latticeEnergyTarg / numAtoms);

        MeterPotentialEnergy meterPERef = new MeterPotentialEnergy(potentialMasterTarg, boxRef);
        latticeEnergyRef = meterPERef.getDataAsScalar();
        System.out.println("lattice energy/N (ref n=" + exponent[0] + "): " + latticeEnergyRef / numAtoms);


        // OVERLAP
        integratorOverlap = new IntegratorOverlap(new IntegratorBox[]{integratorRef, integratorTarg});

        MeterBoltzmann meterTarg = new MeterBoltzmann(integratorTarg, meterPERef, latticeEnergyTarg, latticeEnergyRef);
        meters[1] = meterTarg;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, numAlpha, false), 1);

        MeterBoltzmann meterRef = new MeterBoltzmann(integratorRef, meterPETarg, latticeEnergyRef, latticeEnergyTarg);
        meters[0] = meterRef;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, numAlpha, true), 0);

        setRefPref(alpha, alphaSpan);

        activityIntegrate = new ActivityIntegrate(integratorOverlap);
        getController().addAction(activityIntegrate);
    }

    /**
     * @param args filename containing simulation parameters
     * @see SimOverlapSoftSphereSoftness.SimOverlapParam
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
        int[] exponentN = params.exponentN;
        long numSteps = params.numSteps;
        int numMolecules = params.numMolecules;
        double harmonicFudge = params.harmonicFudge;
        double temperature = params.temperature;
        double alpha = params.alpha;
        double alphaSpan = params.alphaSpan;
        int numAlpha = params.numAlpha;
        int D = params.D;

        numSteps *=2;

        System.out.println("Running "+(D==1 ? "1D" : (D==3 ? "FCC" : "2D hexagonal")) +" soft sphere overlap simulation");
        System.out.println(numMolecules+" atoms at density "+density+" and temperature "+temperature);
        System.out.println("Perturbation from n=" + exponentN[0] +" ("+(1.0/exponentN[0]) +") to n=" + exponentN[1]+" ("+(1.0/exponentN[1]) +")");
        System.out.println((numSteps/1000)+" total steps of 1000\n");


        SimOverlapSoftSphereSoftness sim = new SimOverlapSoftSphereSoftness(Space.getInstance(D), numMolecules, density, temperature,
                harmonicFudge, exponentN, alpha, alphaSpan, numAlpha);

        //start simulation
        sim.integratorOverlap.setNumSubSteps(1000);
        sim.integratorOverlap.setAdjustStepFreq(false);
        numSteps /= 1000;

        System.out.flush();

        sim.equilibrate(numSteps);
        System.out.println("equilibration finished");
        System.out.flush();

        final long startTime = System.currentTimeMillis();
        System.out.println("Start Time: " + startTime);

        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().actionPerformed();

        System.out.println("final reference optimal step frequency "+sim.integratorOverlap.getStepFreq0()
        		+" (actual: "+sim.integratorOverlap.getActualStepFreq0()+")");


        System.out.println("numPoint: " +sim.accumulators[0].getNBennetPoints()+ "\n");


        System.out.println("ratio averages: \n");
        System.out.println(exponentN[0] + " to " + exponentN[1]);
        //System.out.println("    i\t  alpha_set\t       alpha\t         ratio0\t            ratio0_err\t         ratio1\t          ratio1_err");

        for (int i=0; i<numAlpha; i++){

            double ratio = sim.dsvo.getAverage(i);
        	double ratio_err = sim.dsvo.getError(i);

            DataGroup dataRatio0 = (DataGroup)sim.accumulators[0].getData(i);
            double ratio0 = ((DataDoubleArray) dataRatio0.getData(AccumulatorAverage.AVERAGE.index)).getData()[1];
            double ratio0_err = ((DataDoubleArray) dataRatio0.getData(AccumulatorAverage.ERROR.index)).getData()[1];

        	DataGroup dataRatio1 = (DataGroup)sim.accumulators[1].getData(i);
            double ratio1 = ((DataDoubleArray) dataRatio1.getData(AccumulatorAverage.AVERAGE.index)).getData()[1];
            double ratio1_err = ((DataDoubleArray) dataRatio1.getData(AccumulatorAverage.ERROR.index)).getData()[1];

        	System.out.println("    "+sim.accumulators[0].getBennetBias(i)+
        			" "+ ratio  + " " + ratio_err +
                    " " + ratio0 + " " + ratio0_err +
                    " "+ ratio1 + " " + ratio1_err);

        }


        double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
        double ratio = ratioAndError[0];
        double error = ratioAndError[1];

        System.out.println("\nset_alpha: "+alpha);
        System.out.println("ratio_average: "+ratio+" ,error: "+error);
        //System.out.println("free energy difference: "+(-temperature*Math.log(ratio))+" ,error: "+temperature*(error/ratio));

        DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
        System.out.println("ref_ratio_average: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.AVERAGE.index)).getData()[1]
                + " error: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.ERROR.index)).getData()[1]);

        allYourBase = (DataGroup)sim.accumulators[1].getData(sim.dsvo.minDiffLocation());
        System.out.println("targ_ratio_average: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.AVERAGE.index)).getData()[1]
                + " error: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.ERROR.index)).getData()[1]);


        long endTime = System.currentTimeMillis();
        System.out.println("\nEnd Time: " + endTime);
        System.out.println("Time taken (s): " + (endTime - startTime)/1000);

    }

    public void setRefPref(double alpha, double alphaSpan) {
        refPref = alpha;
        accumulators[0].setBennetParam(alpha, alphaSpan);
        accumulators[1].setBennetParam(alpha, alphaSpan);

    }

    public void setAccumulator(AccumulatorVirialOverlapSingleAverage newAccumulator, int iBox) {

        accumulators[iBox] = newAccumulator;

        newAccumulator.setBlockSize(200); // setting the block size = 300

        if (accumulatorPumps[iBox] == null) {
            accumulatorPumps[iBox] = new DataPump(meters[iBox], newAccumulator);
            IntegratorListenerAction pumpListener = new IntegratorListenerAction(accumulatorPumps[iBox]);
            integrators[iBox].getEventManager().addListener(pumpListener);

            pumpListener.setInterval(boxTarg.getMoleculeList().size());
        } else {
            accumulatorPumps[iBox].setDataSink(newAccumulator);
        }

        if (integratorOverlap != null && accumulators[0] != null && accumulators[1] != null) {
            dsvo = new DataSourceVirialOverlap(accumulators[0], accumulators[1]);
            integratorOverlap.setDSVO(dsvo);
        }
    }

    public void equilibrate(long initSteps) {

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

        System.out.println("block size: " + accumulators[0].getBlockSize());
        dsvo.reset();

    }
    
    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules = 108;
        public double density = 1.1964;
        public int[] exponentN = new int[]{12, 14};
        public int D = 3;
        public double alpha =1.0000000000000002;
        public double alphaSpan = 1.0;
        public int numAlpha = 11;
        public long numSteps = 10000000;
        public double harmonicFudge = 1;
        public double temperature = 0.01;
    }
}
