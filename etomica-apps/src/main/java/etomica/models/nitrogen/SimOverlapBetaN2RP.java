/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;


import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.data.DataPump;
import etomica.data.IDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.normalmode.BasisBigCell;
import etomica.normalmode.MCMoveMoleculeCoupled;
import etomica.overlap.IntegratorOverlap;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.Degree;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;

/**
 * 
 *   
 * @author Tai Boon Tan
 */
public class SimOverlapBetaN2RP extends Simulation {

    public SimOverlapBetaN2RP(Space space, int numMolecules, double density, double temperature, 
    		double[] angle, double alpha, double alphaSpan, int numAlpha) {
        super(space);

        SpeciesN2 species = new SpeciesN2(space);
        addSpecies(species);

        potentialMasterTarg = new PotentialMaster();
        potentialMasterRef = new PotentialMaster();

        integrators = new IntegratorBox[2];
        accumulatorPumps = new DataPump[2];
        meters = new IDataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];

        // TARGET
        double ratio = 1.631;
        double aDim = Math.pow(4.0 / (Math.sqrt(3.0) * ratio * density), 1.0 / 3.0);
        double cDim = aDim * ratio;
        System.out.println("\n\naDim: " + aDim + " ;cDim: " + cDim);

        int nC = (int) Math.pow(numMolecules / 1.999999999, 1.0 / 3.0);
        Vector[] boxDim = new Vector[3];
        boxDim[0] = Vector.of(new double[]{nC * aDim, 0, 0});
        boxDim[1] = Vector.of(new double[]{-nC * aDim * Math.cos(Degree.UNIT.toSim(60)), nC * aDim * Math.sin(Degree.UNIT.toSim(60)), 0});
        boxDim[2] = Vector.of(new double[]{0, 0, nC * cDim});
        Boundary boundary = new BoundaryDeformablePeriodic(space, boxDim);
        boxTarg = this.makeBox(boundary);
        boxTarg.setNMolecules(species, numMolecules);


        Basis basisHCP = new BasisHcp();
        BasisBigCell basis = new BasisBigCell(space, basisHCP, new int[]{nC, nC, nC});


        int[] nCells = new int[]{1, 1, 1};
        Primitive primitive = new PrimitiveHexagonal(space, nC * aDim, nC * cDim);

        CoordinateDefinitionNitrogen coordinateDefTarg = new CoordinateDefinitionNitrogen(this, boxTarg, primitive, basis, space);
        coordinateDefTarg.setIsBeta();
        coordinateDefTarg.setOrientationVectorBeta(space);
        coordinateDefTarg.initializeCoordinates(nCells);

        double rCScale = 0.475;
        double rc = aDim * nC * rCScale;
        System.out.println("Truncation Radius (" + rCScale + " Box Length): " + rc);
        P2Nitrogen potentialTarg = new P2Nitrogen(space, rc);
        potentialTarg.setBox(boxTarg);

        PRotConstraint pRotConstraintTarg = new PRotConstraint(space, coordinateDefTarg, boxTarg);
        pRotConstraintTarg.setConstraintAngle(angle[0]);

        potentialMasterTarg.addPotential(potentialTarg, new ISpecies[]{species, species});
        potentialMasterTarg.addPotential(pRotConstraintTarg, new ISpecies[]{species});

        MCMoveMoleculeCoupled moveTarg = new MCMoveMoleculeCoupled(potentialMasterTarg, getRandom(), space);
        moveTarg.setBox(boxTarg);
        moveTarg.setPotential(potentialTarg);

        MCMoveRotateMolecule3D rotateTarg = new MCMoveRotateMolecule3D(potentialMasterTarg, getRandom(), space);
        rotateTarg.setBox(boxTarg);

        IntegratorMC integratorTarg = new IntegratorMC(potentialMasterTarg, getRandom(), temperature, boxTarg);
        integratorTarg.getMoveManager().addMCMove(moveTarg);
        integratorTarg.getMoveManager().addMCMove(rotateTarg);

        integrators[1] = integratorTarg;


        // Reference System
        boxRef = this.makeBox(boundary);
        boxRef.setNMolecules(species, numMolecules);

        CoordinateDefinitionNitrogen coordinateDefRef = new CoordinateDefinitionNitrogen(this, boxRef, primitive, basis, space);
        coordinateDefRef.setIsBeta();
        coordinateDefRef.setOrientationVectorBeta(space);
        coordinateDefRef.initializeCoordinates(nCells);

        P2Nitrogen potentialRef = new P2Nitrogen(space, rc);
        potentialRef.setBox(boxRef);

        PRotConstraint pRotConstraintRef = new PRotConstraint(space, coordinateDefRef, boxRef);
        pRotConstraintRef.setConstraintAngle(angle[1]);

        potentialMasterRef.addPotential(potentialRef, new ISpecies[]{species, species});
        potentialMasterRef.addPotential(pRotConstraintRef, new ISpecies[]{species});

        MCMoveMoleculeCoupled moveRef = new MCMoveMoleculeCoupled(potentialMasterRef, getRandom(), space);
        moveRef.setBox(boxRef);
        moveRef.setPotential(potentialRef);

        MCMoveRotateMolecule3D rotateRef = new MCMoveRotateMolecule3D(potentialMasterRef, getRandom(), space);
        rotateRef.setBox(boxRef);

        IntegratorMC integratorRef = new IntegratorMC(potentialMasterRef, getRandom(), temperature, boxRef);
        integratorRef.getMoveManager().addMCMove(moveRef);
        integratorRef.getMoveManager().addMCMove(rotateRef);

        integrators[0] = integratorRef;

        MeterPotentialEnergy meterPERef = new MeterPotentialEnergy(potentialMasterRef, boxTarg);
        //System.out.println("lattice energy (sim unit): " + meterPERef.getDataAsScalar());

        MeterPotentialEnergy meterPETarg = new MeterPotentialEnergy(potentialMasterTarg, boxRef);
        //System.out.println("lattice energy (sim unit): " + meterPETarg.getDataAsScalar());

        // OVERLAP
        integratorOverlap = new IntegratorOverlap(new IntegratorBox[]{integratorRef, integratorTarg});

        MeterBoltzmann meterTarg = new MeterBoltzmann(integratorTarg, meterPERef);
        meters[1] = meterTarg;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, numAlpha, false), 1);

        MeterBoltzmann meterRef = new MeterBoltzmann(integratorRef, meterPETarg);
        meters[0] = meterRef;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, numAlpha, true), 0);


//        AccumulatorAverageFixed accumulator = new AccumulatorAverageFixed(1);
//        DataPump dataPump = new DataPump(meterRef, accumulator);
//        IntegratorListenerAction listener = new IntegratorListenerAction(dataPump);
//        listener.setInterval(1);
//        integratorRef.getEventManager().addListener(listener);


        setRefPref(alpha, alphaSpan);

        this.getController().addActivity(new ActivityIntegrate(integratorOverlap));
    }

    public void setRefPref(double alpha, double alphaSpan) {
        refPref = alpha;
        accumulators[0].setBennetParam(alpha, alphaSpan);
        accumulators[1].setBennetParam(alpha, alphaSpan);

    }

    public void setAccumulator(AccumulatorVirialOverlapSingleAverage newAccumulator, int iBox) {

        accumulators[iBox] = newAccumulator;
    
        newAccumulator.setBlockSize(100); // setting the block size = 300
        
        if (accumulatorPumps[iBox] == null) {
            accumulatorPumps[iBox] = new DataPump(meters[iBox],newAccumulator);
            IntegratorListenerAction pumpListener = new IntegratorListenerAction(accumulatorPumps[iBox]);
            integrators[iBox].getEventManager().addListener(pumpListener);
           
            pumpListener.setInterval(100);//boxTarg.getMoleculeList().getMoleculeCount());            
        }
        else {
            accumulatorPumps[iBox].setDataSink(newAccumulator);
        }
        
        if (integratorOverlap != null && accumulators[0] != null && accumulators[1] != null) {
            dsvo = new DataSourceVirialOverlap(accumulators[0],accumulators[1]);
            integratorOverlap.setReferenceFracSource(dsvo);
        }
    }

    public void equilibrate(long initSteps) {

        for (int i=0; i<2; i++) {
            if (integrators[i] instanceof IntegratorMC) ((IntegratorMC)integrators[i]).getMoveManager().setEquilibrating(true);
        }
        this.getController().runActivityBlocking(new ActivityIntegrate(this.integratorOverlap), initSteps);
        for (int i=0; i<2; i++) {
            if (integrators[i] instanceof IntegratorMC) ((IntegratorMC)integrators[i]).getMoveManager().setEquilibrating(false);
        }

        System.out.println("block size: "+accumulators[0].getBlockSize());
        dsvo.reset();
        
    }

    /**
     * @param args filename containing simulation parameters
     * @see SimOverlapBetaN2RP.SimOverlapParam
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
        double[] angle = params.angle;
        long numSteps = params.numSteps;
        int numMolecules = params.numMolecules;
        double temperature = params.temperature;
        double alpha = params.alpha;
        double alphaSpan = params.alphaSpan;
        int numAlpha = params.numAlpha;
  
        System.out.println("Running beta-phase Nitrogen RP overlap simulation");
        System.out.println(numMolecules+" molecules at density "+density+" and temperature "+temperature + " K");
        System.out.println("perturbing from angle=" + angle[0] + " into " +angle[1]);
        System.out.println("with " + numSteps + " steps (1000 substeps)");
        
        SimOverlapBetaN2RP sim = new SimOverlapBetaN2RP(Space.getInstance(3), numMolecules, density, Kelvin.UNIT.toSim(temperature), 
        		angle, alpha, alphaSpan, numAlpha);
        
        //start simulation
        
//        sim.integratorOverlap.setNumSubSteps(1000);
//        sim.integratorOverlap.setAdjustStepFreq(false);
//        sim.activityIntegrate.setMaxSteps(100000000);
//        sim.getController().actionPerformed();
//        System.exit(1);
        
        sim.integratorOverlap.setNumSubSteps(1000);
        sim.integratorOverlap.setAdjustStepFraction(false);
        numSteps /= 1000;
        
        sim.equilibrate(numSteps);       
        System.out.println("equilibration finished");
        System.out.flush();
 
        final long startTime = System.currentTimeMillis();
        System.out.println("Start Time: " + startTime);
       
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorOverlap), numSteps);
         
        System.out.println("final reference optimal step frequency "+sim.integratorOverlap.getIdealRefStepFraction()
        		+" (actual: "+sim.integratorOverlap.getRefStepFraction()+")");
        System.out.println("numPoint: " +sim.accumulators[0].getNBennetPoints()+ "\n");
        
        System.out.println("ratio averages: \n");
        System.out.println(angle[0] + " to " + angle[1]);
      
        for (int i=0; i<numAlpha; i++){
        	
//        	double ratio = sim.dsvo.getAverage(i);
//        	double ratio_err = sim.dsvo.getError(i);
        	
        	DataGroup dataRatio0 = (DataGroup)sim.accumulators[0].getData(i);
        	double ratio0 = ((DataDoubleArray)dataRatio0.getData(sim.accumulators[0].AVERAGE.index)).getData()[1];
        	double ratio0_err = ((DataDoubleArray)dataRatio0.getData(sim.accumulators[0].ERROR.index)).getData()[1];
            
        	DataGroup dataRatio1 = (DataGroup)sim.accumulators[1].getData(i);
        	double ratio1 = ((DataDoubleArray)dataRatio1.getData(sim.accumulators[0].AVERAGE.index)).getData()[1];
        	double ratio1_err = ((DataDoubleArray)dataRatio1.getData(sim.accumulators[0].ERROR.index)).getData()[1];
        	
        	System.out.println("    "+sim.accumulators[0].getBennetBias(i)+
        			" "+ (ratio0/ratio1)  +//+ " " + ratio_err +
        			" "+ ratio0 + " " + ratio0_err + 
        			" "+ ratio1 + " " + ratio1_err);
      
        }
        
        
        double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
        double ratio = ratioAndError[0];
        double error = ratioAndError[1];
        
        System.out.println("\nset_alpha: "+alpha);
        System.out.println("ratio_average: "+ratio+" ,error: "+error);
        //System.out.println("free energy difference: "+(-temperature*Math.log(ratio))+" ,error: "+temperature*(error/ratio));
        
        DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
        System.out.println("ref_ratio_average: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].AVERAGE.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].ERROR.index)).getData()[1]);
        
        allYourBase = (DataGroup)sim.accumulators[1].getData(sim.dsvo.minDiffLocation());
        System.out.println("targ_ratio_average: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].AVERAGE.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].ERROR.index)).getData()[1]);
        
        
        long endTime = System.currentTimeMillis();
        System.out.println("\nEnd Time: " + endTime);
        System.out.println("Time taken (s): " + (endTime - startTime)/1000);
       
    }

    private static final long serialVersionUID = 1L;
    public IntegratorOverlap integratorOverlap;
    public DataSourceVirialOverlap dsvo;
    public IntegratorBox[] integrators;

    public Box boxTarg, boxRef;
    public int[] nCells;
    public Basis basis;
    public double refPref;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public DataPump[] accumulatorPumps;
    public IDataSource[] meters;
    protected PotentialMaster potentialMasterTarg, potentialMasterRef;
    
    /**
     * Inner class for parameters understood by the SimOverlapBetaN2RP constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules = 128;
        public double density = 0.025;
        public double[] angle = new double[]{160, 156};
        public int D = 3;
        public double alpha =1.0;
        public double alphaSpan = 1.0;
        public int numAlpha = 11;
        public long numSteps =100000;
        public double temperature = 40;
    }
}
