/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.histogram.HistogramSimple;
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
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * Simulation to sample Bennet's Overlap Region
 */
public class SimBennet extends Simulation {

	private static final String APP_NAME = "Sim Bennet's";
    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;

    public Box box;
    public Boundary boundary;
    public Basis basis;
    public SpeciesGeneral species;
    public NormalModes normalModes;
    public int[] nCells;
    public CoordinateDefinition coordinateDefinition;
    public Primitive primitive;
    public PotentialMasterMonatomic potentialMasterMonatomic;
    public double latticeEnergy;
    public MeterHarmonicEnergy meterHarmonicEnergy;
    public double refPref;
    public SimBennet(Space _space, int numAtoms, double density, double temperature, String filename, int exponent) {
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
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        potentialMasterMonatomic = new PotentialMasterMonatomic(getSpeciesManager());
        double L = Math.pow(4.0 / density, 1.0 / 3.0);
        primitive = new PrimitiveCubic(space, L);
        int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
        boundary = new BoundaryRectangularPeriodic(space, n * L);
        box = this.makeBox(boundary);
        integrator = new IntegratorMC(potentialMasterMonatomic, getRandom(), temperature, box);

        //Target
        box.setNMolecules(species, numAtoms);

        this.getController().addActivity(new ActivityIntegrate(integrator));

        nCells = new int[]{n, n, n};
        basis = new BasisCubicFcc();

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        normalModes = new NormalModesFromFile(filename, D);
        /*
         * nuke this line when it is derivative-based
         */
        normalModes.setTemperature(temperature);
        coordinateDefinition.initializeCoordinates(nCells);

        Potential2SoftSpherical potential = new P2SoftSphere(space, 1.0, 1.0, exponent);
        double truncationRadius = boundary.getBoxSize().getX(0) * 0.495;
        P2SoftSphericalTruncatedShifted pTruncated = new P2SoftSphericalTruncatedShifted(space, potential, truncationRadius);
        AtomType sphereType = species.getLeafType();
        potentialMasterMonatomic.addPotential(pTruncated, new AtomType[]{sphereType, sphereType});

        potentialMasterMonatomic.lrcMaster().setEnabled(false);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMasterMonatomic, box);
        latticeEnergy = meterPE.getDataAsScalar();

        meterHarmonicEnergy = new MeterHarmonicEnergy(coordinateDefinition, normalModes);

        MCMoveAtomCoupledBennet move = new MCMoveAtomCoupledBennet(potentialMasterMonatomic, getRandom(),
                coordinateDefinition, normalModes, getRefPref(), space);
        move.setTemperature(temperature);
        move.setLatticeEnergy(latticeEnergy);
        integrator.getMoveManager().addMCMove(move);


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
            filename = "CB_FCC_n"+exponentN+"_T"+ (int)Math.round(temperature*10);
        }


        System.out.println("Running "+(D==1 ? "1D" : (D==3 ? "FCC" : "2D hexagonal")) +" soft sphere Bennet's overlap perturbation simulation");
        System.out.println(numAtoms+" atoms at density "+density+" and temperature "+temperature);
        System.out.println("exponent N: "+ exponentN );
        System.out.println("total steps: "+ numSteps);
        System.out.println("output data to "+filename);

        //construct simulation
        SimBennet sim = new SimBennet(Space.getInstance(D), numAtoms, density, temperature, filename, exponentN);

        IDataSource[] workMeters = new IDataSource[2];

        // Bennet's Overlap ---> Harmonic
        MeterWorkBennetHarmonic meterWorkBennetHarmonic = new MeterWorkBennetHarmonic(sim.integrator, sim.meterHarmonicEnergy);
        meterWorkBennetHarmonic.setRefPref(sim.refPref);
        meterWorkBennetHarmonic.setLatticeEnergy(sim.latticeEnergy);
        workMeters[0] = meterWorkBennetHarmonic;

        DataFork dataForkHarmonic = new DataFork();
        DataPump pumpHarmonic = new DataPump(workMeters[0], dataForkHarmonic);

        AccumulatorAverageFixed dataAverageHarmonic = new AccumulatorAverageFixed();
        dataForkHarmonic.addDataSink(dataAverageHarmonic);
        AccumulatorHistogram histogramHarmonic = new AccumulatorHistogram(new HistogramSimple(2500, new DoubleRange(-150,100)));
        dataForkHarmonic.addDataSink(histogramHarmonic);

        IntegratorListenerAction pumpHarmonicListener = new IntegratorListenerAction(pumpHarmonic);
        pumpHarmonicListener.setInterval(1);
        sim.integrator.getEventManager().addListener(pumpHarmonicListener);

        //Bennet's Overlap ---> Target
        MeterWorkBennetTarget meterWorkBennetTarget = new MeterWorkBennetTarget(sim.integrator, sim.meterHarmonicEnergy);
        meterWorkBennetTarget.setRefPref(sim.refPref);
        meterWorkBennetTarget.setLatticeEnergy(sim.latticeEnergy);
        workMeters[1] = meterWorkBennetTarget;

        DataFork dataForkTarget = new DataFork();
        DataPump pumpTarget = new DataPump(workMeters[1], dataForkTarget);

        AccumulatorAverageFixed dataAverageTarget = new AccumulatorAverageFixed();
        dataForkTarget.addDataSink(dataAverageTarget);
        AccumulatorHistogram histogramTarget = new AccumulatorHistogram(new HistogramSimple(2500, new DoubleRange(-150,100)));
        dataForkTarget.addDataSink(histogramTarget);

        IntegratorListenerAction pumpTargetListener = new IntegratorListenerAction(pumpTarget);
        pumpTargetListener.setInterval(1000);
        sim.integrator.getEventManager().addListener(pumpTargetListener);

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));

        /*
         *
         */

        double wHarmonic = dataAverageHarmonic.getData().getValue(dataAverageHarmonic.AVERAGE.index);
        double wTarget = dataAverageTarget.getData().getValue(dataAverageTarget.AVERAGE.index);

        double eHarmonic = dataAverageHarmonic.getData().getValue(dataAverageHarmonic.ERROR.index);
        double eTarget = dataAverageTarget.getData().getValue(dataAverageTarget.ERROR.index);

        System.out.println("wBennetHarmonic: "  + wHarmonic  + " ,error: "+ eHarmonic);
        System.out.println("wBennetTarget: " + wTarget + " ,error: "+ eTarget);


        /*
         * To get the entropy (overlap --> target or harmonic system),
         * all you need to do is to take:
         *
         *  sHarmonic = wHarmonic + ln(ratioHarmonicAverage)
         *   sTarget  =  wTarget  + ln( ratioTargetAverage )
         *
         *  both ratioHarmonicAverage and ratioTargetAverage are from SimPhaseOverlapSoftSphere class
         *
         */


        //Target
		DataLogger dataLogger1 = new DataLogger();
		DataTableWriter dataTableWriter1 = new DataTableWriter();
		dataLogger1.setFileName(filename + "_hist_BennTargSim");
		dataLogger1.setDataSink(dataTableWriter1);
		dataTableWriter1.setIncludeHeader(false);
		dataLogger1.putDataInfo(histogramTarget.getDataInfo());

        dataLogger1.setWriteInterval(1);
		dataLogger1.setAppending(false);
		dataLogger1.putData(histogramTarget.getData());
		dataLogger1.closeFile();


        //Harmonic
		DataLogger dataLogger2 = new DataLogger();
		DataTableWriter dataTableWriter2 = new DataTableWriter();
		dataLogger2.setFileName(filename + "_hist_BennHarmSim");
		dataLogger2.setDataSink(dataTableWriter2);
		dataTableWriter2.setIncludeHeader(false);
		dataLogger2.putDataInfo(histogramHarmonic.getDataInfo());

        dataLogger2.setWriteInterval(1);
		dataLogger2.setAppending(false);
		dataLogger2.putData(histogramHarmonic.getData());
		dataLogger2.closeFile();
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
    	public long numSteps = 10000000;
    	public double harmonicFudge =1;
    	public String filename = "CB_FCC_n12_T10";
    	public double temperature = 1.0;
    }

}
