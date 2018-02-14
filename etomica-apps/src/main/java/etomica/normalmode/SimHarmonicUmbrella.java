/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.histogram.HistogramSimple;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.integrator.IntegratorListenerAction;
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
 * Simulation to sample harmonic system and perturbing into
 * 	the umbrella sampling region
 */
public class SimHarmonicUmbrella extends Simulation {

	private static final String APP_NAME = "Sim Harmonic-Umbrella";
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
    public Primitive primitive;
    public PotentialMasterMonatomic potentialMasterMonatomic;
    public double latticeEnergy;
    public MeterHarmonicEnergy meterHarmonicEnergy;
    public double refPref;
    public SimHarmonicUmbrella(Space _space, int numAtoms, double density, double temperature, String filename, int exponent) {
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


        int D = space.D();

        potentialMasterMonatomic = new PotentialMasterMonatomic(this);
        box = new Box(space);
        integrator = new IntegratorMC(potentialMasterMonatomic, getRandom(), temperature, box);
        species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        //Target
        addBox(box);
        box.setNMolecules(species, numAtoms);

        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        double L = Math.pow(4.0 / density, 1.0 / 3.0);
        primitive = new PrimitiveCubic(space, L);
        int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
        nCells = new int[]{n, n, n};
        boundary = new BoundaryRectangularPeriodic(space, n * L);
        basis = new BasisCubicFcc();

        box.setBoundary(boundary);

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        normalModes = new NormalModesFromFile(filename, D);
        /*
         * nuke this line when it is derivative-based
         */
        //normalModes.setTemperature(temperature);
        coordinateDefinition.initializeCoordinates(nCells);

        Potential2SoftSpherical potential = new P2SoftSphere(space, 1.0, 1.0, exponent);
        double truncationRadius = boundary.getBoxSize().getX(0) * 0.495;
        P2SoftSphericalTruncatedShifted pTruncated = new P2SoftSphericalTruncatedShifted(space, potential, truncationRadius);
        AtomType sphereType = species.getLeafType();
        potentialMasterMonatomic.addPotential(pTruncated, new AtomType[]{sphereType, sphereType});

        potentialMasterMonatomic.lrcMaster().setEnabled(false);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMasterMonatomic);
        meterPE.setBox(box);
        latticeEnergy = meterPE.getDataAsScalar();

        MCMoveHarmonic move = new MCMoveHarmonic(getRandom());
        integrator.getMoveManager().addMCMove(move);

        WaveVectorFactory waveVectorFactory = normalModes.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);
        move.setOmegaSquared(normalModes.getOmegaSquared());
        move.setEigenVectors(normalModes.getEigenvectors());
        move.setWaveVectors(waveVectorFactory.getWaveVectors());
        move.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
        move.setCoordinateDefinition(coordinateDefinition);
        move.setTemperature(temperature);

        move.setBox(box);

        meterHarmonicEnergy = new MeterHarmonicEnergy(coordinateDefinition, normalModes);
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


        System.out.println("Running "+(D==1 ? "1D" : (D==3 ? "FCC" : "2D hexagonal")) +" soft sphere Harmonic-Umbrella perturbation simulation");
        System.out.println(numAtoms+" atoms at density "+density+" and temperature "+temperature);
        System.out.println("exponent N: "+ exponentN );
        System.out.println("total steps: "+ numSteps);
        System.out.println("output data to "+filename);

        //construct simulation
        final SimHarmonicUmbrella sim = new SimHarmonicUmbrella(Space.getInstance(D), numAtoms, density, temperature, filename, exponentN);

        IDataSource[] workMeters = new IDataSource[1];

        // Harmonic
        MeterWorkHarmonicUmbrella meterWorkHarmonicUmbrella = new MeterWorkHarmonicUmbrella(sim.integrator, sim.meterHarmonicEnergy);
        meterWorkHarmonicUmbrella.setRefPref(sim.refPref);
        meterWorkHarmonicUmbrella.setLatticeEnergy(sim.latticeEnergy);
        workMeters[0] = meterWorkHarmonicUmbrella;

        DataFork dataFork = new DataFork();
        DataPump pumpHarmonic = new DataPump(workMeters[0], dataFork);

        final AccumulatorAverageFixed dataAverageHarmonic = new AccumulatorAverageFixed();
        dataFork.addDataSink(dataAverageHarmonic);
        IntegratorListenerAction pumpHarmonicListener = new IntegratorListenerAction(pumpHarmonic);
        pumpHarmonicListener.setInterval(1);
        sim.integrator.getEventManager().addListener(pumpHarmonicListener);

        //Histogram Harmonic
        final AccumulatorHistogram histogramHarmonic = new AccumulatorHistogram(new HistogramSimple(600, new DoubleRange(-150,450)));
        dataFork.addDataSink(histogramHarmonic);


       FileWriter fileWriter;

        try{
        	fileWriter = new FileWriter(filename + "_HarmUmb");
        } catch (IOException e){
        	fileWriter = null;
        }

        final String outFileName = filename;
        final FileWriter fileWriterHarmUmb = fileWriter;

        IAction outputAction = new IAction(){
        	public void actionPerformed(){
        		long idStep = sim.integrator.getStepCount();
		        /*
		         * Histogram
		         */
		        //Harmonic
				DataLogger dataLogger = new DataLogger();
				DataTableWriter dataTableWriter = new DataTableWriter();
				dataLogger.setFileName(outFileName + "_hist_HarmUmb");
				dataLogger.setDataSink(dataTableWriter);
				dataTableWriter.setIncludeHeader(false);
				dataLogger.putDataInfo(histogramHarmonic.getDataInfo());

                dataLogger.setWriteInterval(1);
				dataLogger.setAppending(false); //overwrite the file 8/5/08
				dataLogger.putData(histogramHarmonic.getData());
				dataLogger.closeFile();

                System.out.println("\n*****************************************************************");
		        System.out.println("********** Harmonic-to-Umbrella Sampling "+ idStep + "   *************");
		        System.out.println("*****************************************************************");

                double wHarmonic = dataAverageHarmonic.getData().getValue(AccumulatorAverage.AVERAGE.index);
                double eHarmonic = dataAverageHarmonic.getData().getValue(AccumulatorAverage.ERROR.index);
                System.out.println("\n wHarmonicUmbrella: "  + wHarmonic   + " ,error: "+ eHarmonic);

                try {
					fileWriterHarmUmb.write(idStep + " " + wHarmonic + " " + eHarmonic+ "\n");

                } catch(IOException e){

                }
        	}
        };

        IntegratorListenerAction outputActionListener = new IntegratorListenerAction(outputAction);
        outputActionListener.setInterval(20000);
        sim.integrator.getEventManager().addListener(outputActionListener);

        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().actionPerformed();

        try {
	        fileWriterHarmUmb.close();

        } catch(IOException e){

        }
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
    	public String filename = "CB_FCC_n12_T01";
    	public double temperature = 0.1;
    }

}
