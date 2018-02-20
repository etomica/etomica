/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPump;
import etomica.data.IDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.integrator.IntegratorListenerAction;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Simulation to sample Umbrella Overlap Sampling
 * 
 * The original Umbrella Sampling Simulation
 * 	- used to check for the computation time
 * 
 * @author Tai Boon Tan
 */
public class SimUmbrellaSoftSphere extends Simulation {

	private static final String APP_NAME = "Sim Umbrella's";
    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public Box box;
    public Boundary boundary;
    public Basis basis;
    public SpeciesSpheresMono species;
    public NormalModes normalModes;
    public int[] nCells;
    public CoordinateDefinition coordinateDefinition;
    public Primitive primitive, primitiveUnitCell;
    public PotentialMaster potentialMaster;
    public double latticeEnergy;
    public MeterHarmonicEnergy meterHarmonicEnergy;
    public MeterPotentialEnergy meterEnergy;
    public MCMoveAtomCoupledUmbrella move;
    public double refPref;
    public SimUmbrellaSoftSphere(Space _space, int numAtoms, double density, double temperature, String filename, int exponent) {
        super(_space);

        String refFileName = filename + "_ref";
        FileReader refFileReader;
        try {
            refFileReader = new FileReader(refFileName);
        } catch (IOException e) {
            throw new RuntimeException("Cannot find refPref file!! " + e.getMessage());
        }

        try {
            BufferedReader bufReader = new BufferedReader(refFileReader);
            String line = bufReader.readLine();

            refPref = Double.parseDouble(line);
            setRefPref(refPref);

        } catch (IOException e) {
            throw new RuntimeException(" Cannot read from file " + refFileName);
        }
        int D = space.D();

        potentialMaster = new PotentialMasterList(this, space);
        box = this.makeBox();
        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature, box);
        species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        //Target
        box.setNMolecules(species, numAtoms);

        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        double L = Math.pow(4.0 / density, 1.0 / 3.0);
        int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
        primitive = new PrimitiveCubic(space, n * L);
        primitiveUnitCell = new PrimitiveCubic(space, L);
        nCells = new int[]{n, n, n};
        boundary = new BoundaryRectangularPeriodic(space, n * L);
        Basis basisFCC = new BasisCubicFcc();
        basis = new BasisBigCell(space, basisFCC, nCells);

        box.setBoundary(boundary);

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        //String inFile = "inputSSDB"+numAtoms;
        normalModes = new NormalModesFromFile(filename, D);
        /*
         * nuke this line when it is derivative-based
         */
        normalModes.setTemperature(temperature);
        coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});

        Potential2SoftSpherical potential = new P2SoftSphere(space, 1.0, 1.0, exponent);
        double truncationRadius = boundary.getBoxSize().getX(0) * 0.495;

        if (potentialMaster instanceof PotentialMasterList) {
            potential = new P2SoftSphericalTruncated(space, potential, truncationRadius);

        } else {
            potential = new P2SoftSphericalTruncatedShifted(space, potential, truncationRadius);

        }

        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(potential, new AtomType[]{sphereType, sphereType});
        potentialMaster.lrcMaster().setEnabled(false);

        /*
         *  1-body Potential to Constraint the atom from moving too far
         *  	away from its lattice-site
         *
         */

        P1Constraint p1Constraint = new P1Constraint(space, primitiveUnitCell.getSize()[0], box, coordinateDefinition);
        potentialMaster.addPotential(p1Constraint, new AtomType[]{sphereType});
        potentialMaster.lrcMaster().setEnabled(false);

        if (potentialMaster instanceof PotentialMasterList) {
            double neighborRange = truncationRadius;
            int cellRange = 7;
            ((PotentialMasterList) potentialMaster).setRange(neighborRange);
            ((PotentialMasterList) potentialMaster).setCellRange(cellRange); // insanely high, this lets us have neighborRange close to dimensions/2
            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
            ((PotentialMasterList) potentialMaster).getNeighborManager(box).reset();
            int potentialCells = ((PotentialMasterList) potentialMaster).getNbrCellManager(box).getLattice().getSize()[0];
            if (potentialCells < cellRange * 2 + 1) {
                throw new RuntimeException("oops (" + potentialCells + " < " + (cellRange * 2 + 1) + ")");
            }
            if (potentialCells > cellRange * 2 + 1) {
                System.out.println("could probably use a larger truncation radius (" + potentialCells + " > " + (cellRange * 2 + 1) + ")");
            }
            //((P2SoftSphericalTruncated)potential).setTruncationRadius(0.6*boundary.getBoxSize().getX(0));
        }

        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster, box);
        latticeEnergy = meterPE.getDataAsScalar();

        move = new MCMoveAtomCoupledUmbrella(potentialMaster, getRandom(),
                coordinateDefinition, normalModes, getRefPref(), space);
        move.setTemperature(temperature);
        move.setLatticeEnergy(latticeEnergy);
        integrator.getMoveManager().addMCMove(move);
        ((MCMoveStepTracker) move.getTracker()).setNoisyAdjustment(true);

        meterHarmonicEnergy = new MeterHarmonicEnergy(coordinateDefinition, normalModes);

        meterEnergy = new MeterPotentialEnergy(potentialMaster, box);

    }

    /**
     * @param args
     */
    public static void main(String[] args) {

        SimBennetParam params = new SimBennetParam();
        String inputFilename = null;
        if (args.length > 0) {
            inputFilename = args[0];
        }
        if (inputFilename != null) {
            ReadParameters readParameters = new ReadParameters(inputFilename, params);
            readParameters.readParameters();
        }
        double density = params.density/10000;
        int exponentN = params.exponentN;
        long numSteps = params.numSteps;
        final int numAtoms = params.numMolecules;
        double temperature = params.temperature;
        int D = params.D;
        String filename = params.filename;
        if (filename.length() == 0) {
        	System.err.println("Need input files!!!");
            filename = "FCC_SoftSphere_n"+exponentN+"_T"+ (int)Math.round(temperature*10);
        }

        long startTime = System.currentTimeMillis();
        System.out.println("Start Time: " + startTime);
        System.out.println("Running "+(D==1 ? "1D" : (D==3 ? "FCC" : "2D hexagonal")) +" soft sphere Umbrella's Sampling simulation");
        System.out.println(numAtoms+" atoms at density "+density+" and temperature "+temperature);
        System.out.println("exponent N: "+ exponentN );
        System.out.println("total steps: "+ numSteps);
        System.out.println("output data to "+filename);

        //construct simulation
        final SimUmbrellaSoftSphere sim = new SimUmbrellaSoftSphere(Space.getInstance(D), numAtoms, density, temperature, filename, exponentN);

        IDataSource[] samplingMeters = new IDataSource[2];

        sim.activityIntegrate.setMaxSteps(numSteps/10);
        sim.getController().actionPerformed();
        System.out.println("System Equilibrated!");

        sim.getController().reset();

        /*
         *
         */

        // Harmonic Sampling
        final MeterSamplingHarmonic meterSamplingHarmonic = new MeterSamplingHarmonic(sim.integrator, sim.move);
        meterSamplingHarmonic.setRefPref(sim.refPref);
        samplingMeters[0] = meterSamplingHarmonic;

        final AccumulatorAverageFixed dataAverageSamplingHarmonic = new AccumulatorAverageFixed();

        DataPump pumpSamplingHarmonic = new DataPump(samplingMeters[0], dataAverageSamplingHarmonic);
        dataAverageSamplingHarmonic.setBlockSize(200);
        IntegratorListenerAction pumpSamplingHarmonicListener = new IntegratorListenerAction(pumpSamplingHarmonic);
        sim.integrator.getEventManager().addListener(pumpSamplingHarmonicListener);


        // Target Sampling
        final MeterSamplingTarget meterSamplingTarget = new MeterSamplingTarget(sim.integrator, sim.move);
        meterSamplingTarget.setRefPref(sim.refPref);
        samplingMeters[1] = meterSamplingTarget;

        final AccumulatorAverageFixed dataAverageSamplingTarget = new AccumulatorAverageFixed();

        DataPump pumpSamplingTarget = new DataPump(samplingMeters[1], dataAverageSamplingTarget);
        dataAverageSamplingTarget.setBlockSize(200);
        IntegratorListenerAction pumpSamplingTargetListener = new IntegratorListenerAction(pumpSamplingTarget);
        sim.integrator.getEventManager().addListener(pumpSamplingTargetListener);

        if (numAtoms == 32){
            pumpSamplingHarmonicListener.setInterval(100);
            pumpSamplingTargetListener.setInterval(100);
        } else if (numAtoms == 108){
        	pumpSamplingHarmonicListener.setInterval(300);
            pumpSamplingTargetListener.setInterval(300);

        } else if (numAtoms >= 256){
        	pumpSamplingHarmonicListener.setInterval(500);
            pumpSamplingTargetListener.setInterval(500);

        } else {
        	pumpSamplingHarmonicListener.setInterval(1);
            pumpSamplingTargetListener.setInterval(1);
        }

        int totalCells = 1;
        for (int i=0; i<D; i++) {
            totalCells *= sim.nCells[i];
        }

        double  AHarmonic = CalcHarmonicA.doit(sim.normalModes, D, temperature, numAtoms);
        System.out.println("Harmonic-reference free energy, A: "+AHarmonic + " " + AHarmonic/numAtoms);

        FileWriter fileWriter1, fileWriter2;

        try{
        	fileWriter1 = new FileWriter(filename+"_UmblnQharm");
        	fileWriter2 = new FileWriter(filename+"_UmblnQtarg");

        } catch (IOException e){
        	fileWriter1 = null;
        	fileWriter2 = null;
        }
        final FileWriter fileWriterHarm = fileWriter1;
        final FileWriter fileWriterTarg = fileWriter2;


        IAction outputAction = new IAction(){
        	public void actionPerformed(){
        		long idStep = sim.integrator.getStepCount();

                double Qharmonic = meterSamplingHarmonic.getData().getValue(0);
        		double Qtarget = meterSamplingTarget.getData().getValue(0);

                try {
        			fileWriterHarm.write(idStep + " " + Qharmonic +"\n");
        			fileWriterTarg.write(idStep + " " + Qtarget +"\n");

                } catch (IOException e){

                }

            }
        };

        IntegratorListenerAction outputActionListener = new IntegratorListenerAction(outputAction);
        outputActionListener.setInterval(10000);
        sim.integrator.getEventManager().addListener(outputActionListener);

        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().actionPerformed();

        try{
        	fileWriterHarm.close();
        	fileWriterTarg.close();

        } catch (IOException e){

        }

        double Qharmonic = dataAverageSamplingHarmonic.getData().getValue(AccumulatorAverage.AVERAGE.index);
        double Qtarget = dataAverageSamplingTarget.getData().getValue(AccumulatorAverage.AVERAGE.index);

        double eQharmonic = dataAverageSamplingHarmonic.getData().getValue(AccumulatorAverage.ERROR.index);
        double eQtarget = dataAverageSamplingTarget.getData().getValue(AccumulatorAverage.ERROR.index);

	    /*
	     * deltaFE_harmonic: beta*(FE_harmonic - FE_umbrella) = - ln(Qharmonic)
	     *  deltaFE_target : beta*(FE_target - FE_umbrella) = - ln(Qtarget)
	     */
	    double deltaFE_harmonic = - Math.log(Qharmonic);
	    double deltaFE_target = - Math.log(Qtarget);
	    double deltaFE = temperature*deltaFE_target - temperature*deltaFE_harmonic;

        System.out.println("Qharmonic: "+ Qharmonic+ " ;error: " +eQharmonic);
	    System.out.println("Qtarget: "+ Qtarget+ " ;error: " +eQtarget);
	    System.out.println("Q: " +(Qtarget/Qharmonic) + " ;error: "+ (Qtarget/Qharmonic)*Math.sqrt((eQharmonic/Qharmonic)*(eQharmonic/Qharmonic) + (eQtarget/Qtarget)*(eQtarget/Qtarget)));
	    System.out.println("\nfree energy difference (umbrella --> harmonic): "+ (temperature*deltaFE_harmonic) + "; error: " + temperature*eQharmonic/Qharmonic);
	    System.out.println("free energy difference (umbrella --> target): "+ (temperature*deltaFE_target) + "; error: " + temperature*eQtarget/Qtarget);

        System.out.println(" beta*deltaFE (harmonic): " + deltaFE_harmonic);
	    System.out.println(" beta*deltaFE (target): " + deltaFE_target);

        System.out.println("\nHarmonic-reference free energy: " + AHarmonic);
        System.out.println("free energy difference (harmonic --> target): " + deltaFE
                + ", error: " + temperature*Math.sqrt( (eQharmonic/Qharmonic)*(eQharmonic/Qharmonic) + (eQtarget/Qtarget)*(eQtarget/Qtarget) ));
		System.out.println("target free energy: " + (AHarmonic+ deltaFE));
		System.out.println("target free energy per particle: " + (AHarmonic+ deltaFE)/numAtoms);
		long endTime = System.currentTimeMillis();
		System.out.println("End Time: " + endTime);
		System.out.println("Total time taken: " + (endTime - startTime));

    }

    public double getRefPref() {
        return refPref;
    }

    public void setRefPref(double refPref) {
        this.refPref = refPref;
    }
    
    public static class SimBennetParam extends ParameterBase {
    	public int numMolecules = 256;
    	public double density = 12560;
    	public int exponentN = 12;
    	public int D = 3;
    	public long numSteps = 10000;
    	public double harmonicFudge =1;
    	public String filename = "inputSSDB256T01";
    	public double temperature = 0.1;
    }

}
