package etomica.normalmode;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.IEtomicaDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.listener.IntegratorListenerAction;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMaster;
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

    public SimOverlapSoftSphereSoftness(Space _space, int numAtoms, double density, double temperature, String filename, double harmonicFudge, int[] exponent) {
        super(_space);
        this.fname = filename;
        
        potentialMasterTarg = new PotentialMasterList(this, space);
        potentialMasterRef  = new PotentialMasterList(this, space);
        
        integrators = new IntegratorBox[2];
        accumulatorPumps = new DataPump[2];
        meters = new IEtomicaDataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        // TARGET
        boxTarg = new Box(space);
        addBox(boxTarg);
        boxTarg.setNMolecules(species, numAtoms);

        IntegratorMC integratorTarg = new IntegratorMC(potentialMasterTarg, getRandom(), temperature);
        atomMoveTarg = new MCMoveAtomCoupled(potentialMasterTarg, getRandom(), space);
        atomMoveTarg.setStepSize(0.1);
        atomMoveTarg.setStepSizeMax(0.5);
        atomMoveTarg.setDoExcludeNonNeighbors(true);
        integratorTarg.getMoveManager().addMCMove(atomMoveTarg);
        ((MCMoveStepTracker)atomMoveTarg.getTracker()).setNoisyAdjustment(false);
        
        integrators[1] = integratorTarg;

        double L = Math.pow(4.0/density, 1.0/3.0);
        int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
        primitive = new PrimitiveCubic(space, n*L);
        primitiveUnitCell = new PrimitiveCubic(space, L);
        
        nCells = new int[]{n,n,n};
        boundaryTarg = new BoundaryRectangularPeriodic(space, n * L);
        Basis basisFCC = new BasisCubicFcc();
        basis = new BasisBigCell(space, basisFCC, nCells);

        boxTarg.setBoundary(boundaryTarg);

        CoordinateDefinitionLeaf coordinateDefinitionTarg = new CoordinateDefinitionLeaf(boxTarg, primitive, basis, space);
        coordinateDefinitionTarg.initializeCoordinates(new int[]{1,1,1});

        Potential2SoftSpherical potentialTarg = new P2SoftSphere(space, 1.0, 1.0, exponent[1]);
        double truncationRadius = boundaryTarg.getBoxSize().getX(0) * 0.495;
     	if(potentialMasterTarg instanceof PotentialMasterList){
			potentialTarg = new P2SoftSphericalTruncated(space, potentialTarg, truncationRadius);
		
		} else {
			potentialTarg = new P2SoftSphericalTruncatedShifted(space, potentialTarg, truncationRadius);
			
		}
        atomMoveTarg.setPotential(potentialTarg);
        IAtomType sphereType = species.getLeafType();
        potentialMasterTarg.addPotential(potentialTarg, new IAtomType[] {sphereType, sphereType });
        
        /*
         *  1-body Potential to Constraint the atom from moving too far
         *  	away from its lattice-site
         *  
         */

        P1Constraint p1ConstraintTarg = new P1Constraint(space, primitiveUnitCell, boxTarg, coordinateDefinitionTarg);
        potentialMasterTarg.addPotential(p1ConstraintTarg, new IAtomType[] {sphereType});
        potentialMasterTarg.lrcMaster().setEnabled(false);
    
        integratorTarg.setBox(boxTarg);

		if (potentialMasterTarg instanceof PotentialMasterList) {
            double neighborRange = truncationRadius;
            int cellRange = 7;
            ((PotentialMasterList)potentialMasterTarg).setRange(neighborRange);
            ((PotentialMasterList)potentialMasterTarg).setCellRange(cellRange); // insanely high, this lets us have neighborRange close to dimensions/2
            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
            ((PotentialMasterList)potentialMasterTarg).getNeighborManager(boxTarg).reset();
            int potentialCells = ((PotentialMasterList)potentialMasterTarg).getNbrCellManager(boxTarg).getLattice().getSize()[0];
            if (potentialCells < cellRange*2+1) {
                throw new RuntimeException("oops ("+potentialCells+" < "+(cellRange*2+1)+")");
            }
            if (potentialCells > cellRange*2+1) {
                System.out.println("could probably use a larger truncation radius ("+potentialCells+" > "+(cellRange*2+1)+")");
            }
            //((P2SoftSphericalTruncated)potential).setTruncationRadius(0.6*boundaryTarget.getBoxSize().getX(0));
		}
        
        MeterPotentialEnergy meterPETarg = new MeterPotentialEnergy(potentialMasterTarg);
        meterPETarg.setBox(boxTarg);
        latticeEnergyTarg = meterPETarg.getDataAsScalar();
        System.out.println("lattice energy/N (targ n="+exponent[1]+"): " + latticeEnergyTarg/numAtoms);
       
        
        // Reference System
        boundaryRef = new BoundaryRectangularPeriodic(space);
        boxRef = new Box(boundaryRef, space);
        addBox(boxRef);
        boxRef.setNMolecules(species, numAtoms);
        
        IntegratorMC integratorRef = new IntegratorMC(potentialMasterRef, getRandom(), temperature);
        atomMoveRef = new MCMoveAtomCoupled(potentialMasterRef, getRandom(), space);
        atomMoveRef.setStepSize(0.1);
        atomMoveRef.setStepSizeMax(0.5);
        atomMoveRef.setDoExcludeNonNeighbors(true);
        integratorRef.getMoveManager().addMCMove(atomMoveRef);
        ((MCMoveStepTracker)atomMoveRef.getTracker()).setNoisyAdjustment(false);
        
        integrators[0] = integratorRef;
        
        boxRef.setBoundary(boundaryTarg);
        
        CoordinateDefinitionLeaf coordinateDefinitionRef = new CoordinateDefinitionLeaf(boxRef, primitive, basis, space);
        coordinateDefinitionRef.initializeCoordinates(new int[]{1,1,1});

        Potential2SoftSpherical potentialRef = new P2SoftSphere(space, 1.0, 1.0, exponent[0]);
     	if(potentialMasterRef instanceof PotentialMasterList){
			potentialRef = new P2SoftSphericalTruncated(space, potentialRef, truncationRadius);
		
		} else {
			potentialRef = new P2SoftSphericalTruncatedShifted(space, potentialRef, truncationRadius);
			
		}
        atomMoveRef.setPotential(potentialRef);
        potentialMasterRef.addPotential(potentialRef, new IAtomType[] {sphereType, sphereType });
        
        /*
         *  1-body Potential to Constraint the atom from moving too far
         *  	away from its lattice-site
         *  
         */

        P1Constraint p1ConstraintRef = new P1Constraint(space, primitiveUnitCell, boxRef, coordinateDefinitionRef);
        potentialMasterRef.addPotential(p1ConstraintRef, new IAtomType[] {sphereType});
        potentialMasterRef.lrcMaster().setEnabled(false);
    
        integratorRef.setBox(boxRef);

		if (potentialMasterRef instanceof PotentialMasterList) {
            double neighborRange = truncationRadius;
            int cellRange = 7;
            ((PotentialMasterList)potentialMasterRef).setRange(neighborRange);
            ((PotentialMasterList)potentialMasterRef).setCellRange(cellRange); // insanely high, this lets us have neighborRange close to dimensions/2
            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
            ((PotentialMasterList)potentialMasterRef).getNeighborManager(boxRef).reset();
            int potentialCells = ((PotentialMasterList)potentialMasterRef).getNbrCellManager(boxRef).getLattice().getSize()[0];
            if (potentialCells < cellRange*2+1) {
                throw new RuntimeException("oops ("+potentialCells+" < "+(cellRange*2+1)+")");
            }
            if (potentialCells > cellRange*2+1) {
                System.out.println("could probably use a larger truncation radius ("+potentialCells+" > "+(cellRange*2+1)+")");
            }
            //((P2SoftSphericalTruncated)potential).setTruncationRadius(0.6*boundaryTarget.getBoxSize().getX(0));
		}
        
	    MeterPotentialEnergy meterPERef = new MeterPotentialEnergy(potentialMasterRef);
        meterPERef.setBox(boxRef);
        latticeEnergyRef = meterPERef.getDataAsScalar();
        System.out.println("lattice energy/N (ref n="+exponent[0]+"): " + latticeEnergyRef/numAtoms);
       
        
        
        
        // OVERLAP
        integratorOverlap = new IntegratorOverlap(new IntegratorBox[]{integratorRef, integratorTarg});

        MeterBoltzmann meterTarg = new MeterBoltzmann(integratorTarg, meterPERef, latticeEnergyTarg, latticeEnergyRef);
        meters[1] = meterTarg;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, false), 1);
        
        MeterBoltzmann meterRef = new MeterBoltzmann(integratorRef, meterPETarg, latticeEnergyRef, latticeEnergyTarg);
        meters[0] = meterRef;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, true), 0);
        
        setRefPref(1.0, 2.0);
        
        activityIntegrate = new ActivityIntegrate(integratorOverlap);
        
        getController().addAction(activityIntegrate);
    }

    public void setRefPref(double refPrefCenter, double span) {
        refPref = refPrefCenter;
        accumulators[0].setBennetParam(refPrefCenter,span);
        accumulators[1].setBennetParam(refPrefCenter,span);
        // needed for Bennett sampling
//        if (accumulators[0].getNBennetPoints() == 1) {
//            ((MeterBoltzmannHarmonic)meters[0]).refPref = refPrefCenter;
//            ((MeterBoltzmannTarget)meters[1]).refPref = refPrefCenter;
//        }
    }

    public void setAccumulator(AccumulatorVirialOverlapSingleAverage newAccumulator, int iBox) {

        accumulators[iBox] = newAccumulator;
    
        newAccumulator.setBlockSize(200); // setting the block size = 300
        
        if (accumulatorPumps[iBox] == null) {
            accumulatorPumps[iBox] = new DataPump(meters[iBox],newAccumulator);
            IntegratorListenerAction pumpListener = new IntegratorListenerAction(accumulatorPumps[iBox]);
            integrators[iBox].getEventManager().addListener(pumpListener);
            if (iBox == 1) {
            	if (boxTarg.getMoleculeList().getMoleculeCount()==32){
            		
            	    pumpListener.setInterval(100);
            	
            	} else if (boxTarg.getMoleculeList().getMoleculeCount()==108){
            	    pumpListener.setInterval(300);
            	    
            	} else if (boxTarg.getMoleculeList().getMoleculeCount()==256){
            	    pumpListener.setInterval(500);
            	
            	} else 
            		
            	    pumpListener.setInterval(boxTarg.getMoleculeList().getMoleculeCount());
            }
        }
        else {
            accumulatorPumps[iBox].setDataSink(newAccumulator);
        }
        if (integratorOverlap != null && accumulators[0] != null && accumulators[1] != null) {
            dsvo = new DataSourceVirialOverlap(accumulators[0],accumulators[1]);
            integratorOverlap.setDSVO(dsvo);
        }
    }
    
    public void setRefPref(double newRefPref) {
        System.out.println("setting ref pref (explicitly) to "+newRefPref);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(11,true),0);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(11,false),1);
        setRefPref(newRefPref,1);
    }
    
    public void initRefPref(String fileName, long initSteps) {
        // refPref = -1 indicates we are searching for an appropriate value
        refPref = -1.0;
        if (fileName != null) {
            try { 
                FileReader fileReader = new FileReader(fileName);
                BufferedReader bufReader = new BufferedReader(fileReader);
                String refPrefString = bufReader.readLine();
                refPref = Double.parseDouble(refPrefString);
                bufReader.close();
                fileReader.close();
                System.out.println("setting ref pref (from file) to "+refPref);
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,true),0);
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,false),1);
                setRefPref(refPref,2.0);
            }
            catch (IOException e) {
                // file not there, which is ok.
            }
        }
        
        if (refPref == -1) {
            // equilibrate off the lattice to avoid anomolous contributions
            activityIntegrate.setMaxSteps(initSteps/2);
            getController().actionPerformed();
            getController().reset();
            System.out.println("target equilibration finished");

            setAccumulator(new AccumulatorVirialOverlapSingleAverage(41,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(41,false),1);
            setRefPref(1,200);
            activityIntegrate.setMaxSteps(initSteps);
            getController().actionPerformed();
            getController().reset();

            int newMinDiffLoc = dsvo.minDiffLocation();
            refPref = accumulators[0].getBennetAverage(newMinDiffLoc)
                /accumulators[1].getBennetAverage(newMinDiffLoc);
            if (Double.isNaN(refPref) || refPref == 0 || Double.isInfinite(refPref)) {
                throw new RuntimeException("Simulation failed to fiq" +
                		"nd a valid ref pref");
            }
            System.out.println("setting ref pref to "+refPref);
            
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(11,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(11,false),1);
            setRefPref(refPref,5);

            // set refPref back to -1 so that later on we know that we've been looking for
            // the appropriate value
            refPref = -1;
            getController().reset();
        }

    }
    
    public void equilibrate(String fileName, long initSteps) {
        // run a short simulation to get reasonable MC Move step sizes and
        // (if needed) narrow in on a reference preference
        activityIntegrate.setMaxSteps(initSteps);

        for (int i=0; i<2; i++) {
            if (integrators[i] instanceof IntegratorMC) ((IntegratorMC)integrators[i]).getMoveManager().setEquilibrating(true);
        }
        getController().actionPerformed();
        getController().reset();
        for (int i=0; i<2; i++) {
            if (integrators[i] instanceof IntegratorMC) ((IntegratorMC)integrators[i]).getMoveManager().setEquilibrating(false);
        }

        if (refPref == -1) {
            int newMinDiffLoc = dsvo.minDiffLocation();
            refPref = accumulators[0].getBennetAverage(newMinDiffLoc)
                /accumulators[1].getBennetAverage(newMinDiffLoc);
            System.out.println("setting ref pref to "+refPref+" ("+newMinDiffLoc+")");
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,true),0);
            
            System.out.println("block size (equilibrate) "+accumulators[0].getBlockSize());
            
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,false),1);
            setRefPref(refPref,2.0);
            if (fileName != null) {
                try {
                    FileWriter fileWriter = new FileWriter(fileName);
                    BufferedWriter bufWriter = new BufferedWriter(fileWriter);
                    bufWriter.write(String.valueOf(refPref)+"\n");
                    bufWriter.close();
                    fileWriter.close();
                }
                catch (IOException e) {
                    throw new RuntimeException("couldn't write to refpref file");
                }
            }
        }
        else {
            dsvo.reset();
        }
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
        final int numMolecules = params.numMolecules;
        double harmonicFudge = params.harmonicFudge;
        double temperature = params.temperature;
        int D = params.D;
        String filename = params.filename;
        if (filename.length() == 0) {
        	System.err.println("Need input files!!!");
            filename = "CB_FCC_n"+exponentN+"_T"+ (int)Math.round(temperature*10);
        }
        String refFileName = filename+"_ref";
        
        
        System.out.println("Running "+(D==1 ? "1D" : (D==3 ? "FCC" : "2D hexagonal")) +" soft sphere overlap simulation");
        System.out.println(numMolecules+" atoms at density "+density+" and temperature "+temperature);
        System.out.println("Perturbation from n=" + exponentN[0] +" to n=" + exponentN[1]);
        System.out.println((numSteps/1000)+" total steps of 1000");
        System.out.println("output data to "+filename);

        //instantiate simulation
        final SimOverlapSoftSphereSoftness sim = new SimOverlapSoftSphereSoftness(Space.getInstance(D), numMolecules, density, temperature, filename, harmonicFudge, exponentN);
        
        //start simulation
        sim.integratorOverlap.setNumSubSteps(1000);
        numSteps /= 1000;

        sim.initRefPref(refFileName, numSteps/20);
        if (Double.isNaN(sim.refPref) || sim.refPref == 0 || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("Simulation failed to find a valid ref pref");
        }
        System.out.flush();
        
        sim.equilibrate(refFileName, numSteps/10);
        if (Double.isNaN(sim.refPref) || sim.refPref == 0 || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("Simulation failed to find a valid ref pref");
        }
       
        System.out.println("equilibration finished");
        System.out.flush();
 
 
        final long startTime = System.currentTimeMillis();
        System.out.println("Start Time: " + startTime);
       
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().actionPerformed();
         
        System.out.println("final reference optimal step frequency "+sim.integratorOverlap.getStepFreq0()
        		+" (actual: "+sim.integratorOverlap.getActualStepFreq0()+")");
              
        double ratio = sim.dsvo.getDataAsScalar();
        double error = sim.dsvo.getError();
        
        System.out.println("\nratio average: "+ratio+" ,error: "+error);
        System.out.println("free energy difference: "+(-temperature*Math.log(ratio))+" ,error: "+temperature*(error/ratio));
        
        DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
        double betaFAW = -Math.log(((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]);
        System.out.println("ref ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1]);
        
        allYourBase = (DataGroup)sim.accumulators[1].getData(sim.dsvo.minDiffLocation());
        double betaFBW = -Math.log(((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]);
        System.out.println("targ ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1]);
        
        
        long endTime = System.currentTimeMillis();
        System.out.println("End Time: " + endTime);
        System.out.println("Time taken (s): " + (endTime - startTime)/1000);
       
    
    }

    private static final long serialVersionUID = 1L;
    public IntegratorOverlap integratorOverlap;
    public DataSourceVirialOverlap dsvo;
    public IntegratorBox[] integrators;
    public ActivityIntegrate activityIntegrate;
    public IBox boxTarg, boxRef;
    public Boundary boundaryTarg, boundaryRef;
    public int[] nCells;
    public Basis basis;
    public Primitive primitive, primitiveUnitCell;
    public double refPref;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public DataPump[] accumulatorPumps;
    public IEtomicaDataSource[] meters;
    public String fname;
    protected MCMoveAtomCoupled atomMoveTarg, atomMoveRef;
    protected PotentialMaster potentialMasterTarg, potentialMasterRef;
    protected double latticeEnergyTarg, latticeEnergyRef;
    
    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules =500;
        public double density = 1.1964;
        public int[] exponentN = new int[]{10, 12};
        public int D = 3;
        public long numSteps = 1000000;
        public double harmonicFudge = 1;
        public String filename = "inputSSDB32T0001";
        public double temperature =0.1;
    }
}
