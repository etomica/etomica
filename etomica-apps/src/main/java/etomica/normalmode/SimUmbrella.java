/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.IAction;

import etomica.action.activity.ActivityIntegrate2;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.histogram.HistogramExpanding;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.math.DoubleRange;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMasterMonatomic;
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
 * Simulation to sample Umbrella Sampling
 */
public class SimUmbrella extends Simulation {

	private static final String APP_NAME = "Sim Umbrella's";
    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;

    public Box box;
    public Boundary boundary;
    public Basis basis;
    public SpeciesSpheresMono species;
    public NormalModes normalModes;
    public int[] nCells;
    public CoordinateDefinition coordinateDefinition;
    public Primitive primitive;
    public PotentialMasterMonatomic potentialMasterMonatomic;
    public double latticeEnergy;
    public MeterHarmonicEnergy meterHarmonicEnergy;
    public MeterPotentialEnergy meterEnergy;
    public double refPref;
    public MCMoveAtomCoupledUmbrella move;
    public SimUmbrella(Space _space, int numAtoms, double density, double temperature, String filename, int exponent) {
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
        //System.out.println("refPref is: "+ refPref);

        species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        int D = space.D();

        potentialMasterMonatomic = new PotentialMasterMonatomic(this);
        double L = Math.pow(4.0 / density, 1.0 / 3.0);
        primitive = new PrimitiveCubic(space, L);
        int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
        boundary = new BoundaryRectangularPeriodic(space,n * L);
        box = this.makeBox(boundary);
        integrator = new IntegratorMC(potentialMasterMonatomic, getRandom(), temperature, box);

        //Target
        box.setNMolecules(species, numAtoms);

        this.getController2().addActivity(new ActivityIntegrate2(integrator));

        nCells = new int[]{n, n, n};
        basis = new BasisCubicFcc();

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        normalModes = new NormalModesFromFile(filename, D);
        /*
         * nuke this line when it is derivative-based
         */
        normalModes.setTemperature(temperature);
        coordinateDefinition.initializeCoordinates(nCells);

        Potential2SoftSpherical potential = new P2SoftSphere(space);
        double truncationRadius = boundary.getBoxSize().getX(0) * 0.495;
        P2SoftSphericalTruncatedShifted pTruncated = new P2SoftSphericalTruncatedShifted(space, potential, truncationRadius);
        AtomType sphereType = species.getLeafType();
        potentialMasterMonatomic.addPotential(pTruncated, new AtomType[]{sphereType, sphereType});

        /*
         *  1-body Potential to Constraint the atom from moving too far
         *  	away from its lattice-site
         *
         */

        P1Constraint p1Constraint = new P1Constraint(space, primitive.getSize()[0], box, coordinateDefinition);
        potentialMasterMonatomic.addPotential(p1Constraint, new AtomType[]{sphereType});


        potentialMasterMonatomic.lrcMaster().setEnabled(false);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMasterMonatomic, box);
        latticeEnergy = meterPE.getDataAsScalar();

        move = new MCMoveAtomCoupledUmbrella(potentialMasterMonatomic, getRandom(),
                coordinateDefinition, normalModes, getRefPref(), space);
        move.setTemperature(temperature);
        move.setLatticeEnergy(latticeEnergy);
        integrator.getMoveManager().addMCMove(move);

        meterHarmonicEnergy = new MeterHarmonicEnergy(coordinateDefinition, normalModes);

        meterEnergy = new MeterPotentialEnergy(potentialMasterMonatomic, box);

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
        double density = params.density/1000;
        int exponentN = params.exponentN;
        long numSteps = params.numSteps;
        int numAtoms = params.numMolecules;
        double temperature = params.temperature;
        double harmonicFudge = params.harmonicFudge;
        int D = params.D;
        String filename = params.filename;
        if (filename.length() == 0) {
        	System.err.println("Need input files!!!");
            filename = "FCC_SoftSphere_n"+exponentN+"_T"+ (int)Math.round(temperature*10);
        }


        System.out.println("Running "+(D==1 ? "1D" : (D==3 ? "FCC" : "2D hexagonal")) +" soft sphere Umbrella's Sampling perturbation simulation");
        System.out.println(numAtoms+" atoms at density "+density+" and temperature "+temperature);
        System.out.println("exponent N: "+ exponentN );
        System.out.println("total steps: "+ numSteps);
        System.out.println("output data to "+filename);

        //construct simulation
        final SimUmbrella sim = new SimUmbrella(Space.getInstance(D), numAtoms, density, temperature, filename, exponentN);

        IDataSource[] workMeters = new IDataSource[2];
        IDataSource[] samplingMeters = new IDataSource[2];

      // Work Umbrella Sampling ---> Harmonic
        final MeterWorkUmbrellaHarmonic meterWorkUmbrellaHarmonic = new MeterWorkUmbrellaHarmonic(sim.integrator, sim.move);
        meterWorkUmbrellaHarmonic.setRefPref(sim.refPref);
        workMeters[0] = meterWorkUmbrellaHarmonic;

        DataFork dataForkHarmonic = new DataFork();
        DataPump pumpHarmonic = new DataPump(workMeters[0], dataForkHarmonic);

        final AccumulatorAverageFixed dataAverageHarmonic = new AccumulatorAverageFixed();
        dataAverageHarmonic.setBlockSize(100);
        dataForkHarmonic.addDataSink(dataAverageHarmonic);
        IntegratorListenerAction pumpHarmonicListener = new IntegratorListenerAction(pumpHarmonic);
        sim.integrator.getEventManager().addListener(pumpHarmonicListener);

        //Histogram Harmonic
        final AccumulatorHistogram histogramHarmonic = new AccumulatorHistogram(new HistogramExpanding(1, new DoubleRange(0,1)));
        dataForkHarmonic.addDataSink(histogramHarmonic);

        // Work Umbrella Sampling ---> Target
        final MeterWorkUmbrellaTarget meterWorkUmbrellaTarget = new MeterWorkUmbrellaTarget(sim.integrator, sim.move);
        meterWorkUmbrellaTarget.setRefPref(sim.refPref);
        workMeters[1] = meterWorkUmbrellaTarget;

        DataFork dataForkTarget = new DataFork();
        DataPump pumpTarget = new DataPump(workMeters[1], dataForkTarget);

        final AccumulatorAverageFixed dataAverageTarget = new AccumulatorAverageFixed();
        dataAverageTarget.setBlockSize(100);
        dataForkTarget.addDataSink(dataAverageTarget);
        IntegratorListenerAction pumpTargetListener = new IntegratorListenerAction(pumpTarget);
        sim.integrator.getEventManager().addListener(pumpTargetListener);


        //Histogram Target
        final AccumulatorHistogram histogramTarget = new AccumulatorHistogram(new HistogramExpanding(1,new DoubleRange(0,1)));
        dataForkTarget.addDataSink(histogramTarget);

        /*
         *
         */

        // Umbrella Sampling --> Harmonic System

        final MeterSamplingHarmonic meterSamplingHarmonic = new MeterSamplingHarmonic(sim.integrator, sim.move);
        meterSamplingHarmonic.setRefPref(sim.refPref);
        samplingMeters[0] = meterSamplingHarmonic;

        final AccumulatorAverageFixed dataAverageSamplingHarmonic = new AccumulatorAverageFixed();
        dataAverageSamplingHarmonic.setBlockSize(100);
        DataPump pumpSamplingHarmonic = new DataPump(samplingMeters[0], dataAverageSamplingHarmonic);
        IntegratorListenerAction pumpSamplingHarmonicListener = new IntegratorListenerAction(pumpSamplingHarmonic);
        sim.integrator.getEventManager().addListener(pumpSamplingHarmonicListener);

        // Umbrella Sampling --> Target System

        final MeterSamplingTarget meterSamplingTarget = new MeterSamplingTarget(sim.integrator, sim.move);
        meterSamplingTarget.setRefPref(sim.refPref);
        samplingMeters[1] = meterSamplingTarget;

        final AccumulatorAverageFixed dataAverageSamplingTarget = new AccumulatorAverageFixed();
        dataAverageSamplingTarget.setBlockSize(100);
        DataPump pumpSamplingTarget = new DataPump(samplingMeters[1], dataAverageSamplingTarget);
        IntegratorListenerAction pumpSamplingTargetListener = new IntegratorListenerAction(pumpSamplingTarget);
        sim.integrator.getEventManager().addListener(pumpSamplingTargetListener);


        if (numAtoms == 32){
            pumpHarmonicListener.setInterval(100);
            pumpTargetListener.setInterval(100);
            pumpSamplingHarmonicListener.setInterval(100);
            pumpSamplingTargetListener.setInterval(100);
        } else if (numAtoms == 108){
            pumpHarmonicListener.setInterval(300);
            pumpTargetListener.setInterval(300);
            pumpSamplingHarmonicListener.setInterval(300);
            pumpSamplingTargetListener.setInterval(300);
            if (temperature >= 1.1){
            	dataAverageHarmonic.setBlockSize(200);
            	dataAverageTarget.setBlockSize(200);
            	dataAverageSamplingHarmonic.setBlockSize(200);
            	dataAverageSamplingTarget.setBlockSize(200);
            }
        } else {
            pumpHarmonicListener.setInterval(1);
            pumpTargetListener.setInterval(1);
            pumpSamplingHarmonicListener.setInterval(1);
            pumpSamplingTargetListener.setInterval(1);
        }


        double[][] omega2 = sim.normalModes.getOmegaSquared();
		double[] coeffs = sim.normalModes.getWaveVectorFactory().getCoefficients();
		double AHarmonic = 0;
		for(int i=0; i<omega2.length; i++) {
		      for(int j=0; j<omega2[0].length; j++) {
		            if (!Double.isInfinite(omega2[i][j])) {
		               AHarmonic += coeffs[i]*Math.log(omega2[i][j]*coeffs[i]/(temperature*Math.PI));
		            }
		      }
		}

        int totalCells = 1;
        for (int i=0; i<D; i++) {
        		totalCells *= sim.nCells[i];
        }
		int basisSize = sim.basis.getScaledCoordinates().length;
		double fac = 1;
		if (totalCells % 2 == 0) {
			fac = Math.pow(2,D);
		}
		AHarmonic -= Math.log(Math.pow(2.0, basisSize*D*(totalCells - fac)/2.0) / Math.pow(totalCells,0.5*D));
		System.out.println("Harmonic-reference free energy: "+AHarmonic*temperature);

        FileWriter fileWriter1, fileWriter2;

        try{
        	fileWriter1 = new FileWriter(filename+"_SimUmb");
        	fileWriter2 = new FileWriter(filename+"_SimUmbPI");
        } catch(IOException e){
        	fileWriter1 = null;
        	fileWriter2 = null;
        }

        final String outFileName = filename;
        final double temp = temperature;
        final long steps = numSteps;
        final FileWriter fileWriterSimUmb = fileWriter1;
        final FileWriter fileWriterSimUmbPI = fileWriter2;


        IAction outputAction = new IAction(){
        	public void actionPerformed(){
        		long idStep = sim.integrator.getStepCount();


                double wHarmonic = dataAverageHarmonic.getData().getValue(dataAverageHarmonic.AVERAGE.index);
                double wTarget = dataAverageTarget.getData().getValue(dataAverageTarget.AVERAGE.index);

		        /*
                 * Qharmonic = < e0 / [sqrt(e1^2 + alpha^2 * e0^2)]>umbrella
		         *  Qtarget  = < e1 / [sqrt(e1^2 + alpha^2 * e0^2)]>umbrella
		         *
		         */
                double Qharmonic = dataAverageSamplingHarmonic.getData().getValue(dataAverageSamplingHarmonic.AVERAGE.index);
                double Qtarget = dataAverageSamplingTarget.getData().getValue(dataAverageSamplingTarget.AVERAGE.index);

                double eQharmonic = dataAverageSamplingHarmonic.getData().getValue(dataAverageSamplingHarmonic.ERROR.index);
                double eQtarget = dataAverageSamplingTarget.getData().getValue(dataAverageSamplingTarget.ERROR.index);

		        /*
		         * deltaFE_harmonic: beta*(FE_harmonic - FE_umbrella) = - ln(Qharmonic)
		         *  deltaFE_target : beta*(FE_target - FE_umbrella) = - ln(Qtarget)
		         */
		        double beta_deltaFE_harmonic = - Math.log(Qharmonic);
		        double beta_deltaFE_target = - Math.log(Qtarget);
		        double beta_deltaFE = beta_deltaFE_target - beta_deltaFE_harmonic;

                double sHarmonic = wHarmonic - beta_deltaFE_harmonic;
		        double sTarget = wTarget - beta_deltaFE_target;

                double reweightedwHarmonic = meterWorkUmbrellaHarmonic.getDataReweighted();
		        double reweightedwTarget = meterWorkUmbrellaTarget.getDataReweighted();

                double sRHarmonic = - reweightedwHarmonic + beta_deltaFE_harmonic;
		        double sRTarget = - reweightedwTarget + beta_deltaFE_target;

                double M_harm = Math.sqrt(4*Math.PI*sRHarmonic*Math.exp(2*sRHarmonic))+1;
		        double M_targ = Math.sqrt(4*Math.PI*sRTarget*Math.exp(2*sRTarget))+1;
		        try {
		        	fileWriterSimUmb.write( idStep + " " + beta_deltaFE_harmonic + " " + (eQharmonic/Qharmonic)+ " "
		        									  + beta_deltaFE_target + " " + (eQtarget/Qtarget) + " "
		        									  + beta_deltaFE + " "
		        			                          + wHarmonic + " " + wTarget + " "
		        			                          + sHarmonic + " " + sTarget +  "\n");

                    fileWriterSimUmbPI.write(idStep + " " + reweightedwHarmonic + " " + reweightedwTarget + " "
		        									+ " " + beta_deltaFE_harmonic + " " + beta_deltaFE_target + " "
		        									+ " " + sRHarmonic + " " + sRTarget + " "
		        									+ " " + M_harm + " " + M_targ + "\n");

                } catch (IOException e){

                }

            }
        };

        IntegratorListenerAction outputActionListener = new IntegratorListenerAction(outputAction);
        outputActionListener.setInterval((int)numSteps/100);
        sim.integrator.getEventManager().addListener(outputActionListener);

        sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integrator), numSteps);

        try{
	        fileWriterSimUmb.close();
	        fileWriterSimUmbPI.close();
        } catch (IOException e){

        }


        double wHarmonic = dataAverageHarmonic.getData().getValue(dataAverageHarmonic.AVERAGE.index);
        double wTarget = dataAverageTarget.getData().getValue(dataAverageTarget.AVERAGE.index);

        double eHarmonic = dataAverageHarmonic.getData().getValue(dataAverageHarmonic.ERROR.index);
        double eTarget = dataAverageTarget.getData().getValue(dataAverageTarget.ERROR.index);
        /*
         * Qharmonic = < e0 / [sqrt(e1^2 + alpha^2 * e0^2)]>umbrella
         *  Qtarget  = < e1 / [sqrt(e1^2 + alpha^2 * e0^2)]>umbrella
         *
         */
        double Qharmonic = dataAverageSamplingHarmonic.getData().getValue(dataAverageSamplingHarmonic.AVERAGE.index);
        double Qtarget = dataAverageSamplingTarget.getData().getValue(dataAverageSamplingTarget.AVERAGE.index);

        double eQharmonic = dataAverageSamplingHarmonic.getData().getValue(dataAverageSamplingHarmonic.ERROR.index);
        double eQtarget = dataAverageSamplingTarget.getData().getValue(dataAverageSamplingTarget.ERROR.index);

        /*
         * deltaFE_harmonic: beta*(FE_harmonic - FE_umbrella) = - ln(Qharmonic)
         *  deltaFE_target : beta*(FE_target - FE_umbrella) = - ln(Qtarget)
         */
        double deltaFE_harmonic = - Math.log(Qharmonic);
        double deltaFE_target = - Math.log(Qtarget);
        double deltaFE = temp*deltaFE_target - temp*deltaFE_harmonic;

        double sHarmonic = wHarmonic - deltaFE_harmonic;
        double sTarget = wTarget - deltaFE_target;

        double er_sHarmonic = Math.sqrt(eHarmonic*eHarmonic + (eQharmonic/Qharmonic)*(eQharmonic/Qharmonic));
        double er_sTarget   = Math.sqrt(eTarget*eTarget + (eQtarget/Qtarget)*(eQtarget/Qtarget));

        double reweightedwHarmonic = meterWorkUmbrellaHarmonic.getDataReweighted();
        double reweightedwTarget = meterWorkUmbrellaTarget.getDataReweighted();

        System.out.println("\nQharmonic: "+ Qharmonic+ " ;error: " +eQharmonic);
	    System.out.println("Qtarget: "+ Qtarget+ " ;error: " +eQtarget);
	    System.out.println("Q: " +(Qtarget/Qharmonic) + " ;error: "+ (Qtarget/Qharmonic)*Math.sqrt((eQharmonic/Qharmonic)*(eQharmonic/Qharmonic) + (eQtarget/Qtarget)*(eQtarget/Qtarget)));

        System.out.println("\nfree energy difference (umbrella --> harmonic): "+ (temp*deltaFE_harmonic) + "; error: " + temp*eQharmonic/Qharmonic);
        System.out.println("free energy difference (umbrella --> target): "+ (temp*deltaFE_target) + "; error: " + temp*eQtarget/Qtarget);
        System.out.println("free energy difference (harmonic --> target): " + deltaFE
                + ", error: " + temp*Math.sqrt( (eQharmonic/Qharmonic)*(eQharmonic/Qharmonic) + (eQtarget/Qtarget)*(eQtarget/Qtarget) ));
        System.out.println("target free energy: " + (temperature * AHarmonic + deltaFE) +
                " ,error: " + temp*Math.sqrt( (eQharmonic/Qharmonic)*(eQharmonic/Qharmonic) + (eQtarget/Qtarget)*(eQtarget/Qtarget) ));
        System.out.println("target free energy per particle: " + (temperature*AHarmonic+ deltaFE)/numAtoms);


        System.out.println("\nwUmbrellaHarmonic: "+ wHarmonic + ", error: "+ eHarmonic);
        System.out.println(" wUmbrellaTarget : "  + wTarget   + ", error: "+ eTarget);

        System.out.println("\nwHarmonicUmbrella: "+ reweightedwHarmonic );
        System.out.println(" wTargetUmbrella : "  + reweightedwTarget   );

        System.out.println(" beta*deltaFE (harmonic): " + deltaFE_harmonic);
        System.out.println(" beta*deltaFE (target): " + deltaFE_target);

        System.out.println("\nUmbrella (perturbed into Harmonic) entropy: " + sHarmonic + ", error: "+ er_sHarmonic);
        System.out.println("Umbrella (perturbed into Target) entropy: " + sTarget + ", error: "+ er_sTarget);

        /*
         * Histogram
         */
        //Target
		DataLogger dataLogger = new DataLogger();
		DataTableWriter dataTableWriter = new DataTableWriter();
		dataLogger.setFileName(outFileName + "_hist_UmbTarg");
		dataLogger.setDataSink(dataTableWriter);
		dataTableWriter.setIncludeHeader(false);
		dataLogger.putDataInfo(histogramTarget.getDataInfo());

        dataLogger.setWriteInterval(1);
		dataLogger.setAppending(false); //overwrite the file 8/5/08
		dataLogger.putData(histogramTarget.getData());
		dataLogger.closeFile();

        //Harmonic
        dataLogger.setFileName(outFileName + "_hist_UmbHarm");
        dataLogger.setAppending(false); // overwrite the file 8/5/08
        dataLogger.putData(histogramHarmonic.getData());
        dataLogger.closeFile();


    }

    public double getRefPref() {
        return refPref;
    }

    public void setRefPref(double refPref) {
        this.refPref = refPref;
    }
    
    public static class SimBennetParam extends ParameterBase {
    	public int numMolecules = 32;
    	public double density = 1256;
    	public int exponentN = 12;
    	public int D = 3;
    	public long numSteps = 1000000;
    	public double harmonicFudge =1;
    	public String filename = "CB_FCC_n12_T14";
    	public double temperature = 1.4;
    }

}
