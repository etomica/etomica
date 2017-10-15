/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.gaussianwork;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.IController;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.histogram.HistogramSimple;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorMC;
import etomica.listener.IntegratorListenerAction;
import etomica.math.DoubleRange;
import etomica.overlap.IntegratorOverlap;
import etomica.potential.P1Harmonic;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space1d.Space1D;
import etomica.space1d.Vector1D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;



/**
 * MultiHarmonic Overlap Simulation
 * 
 * 	a. To study the probability distribution the reverse-direction perturbation
 * 		by reweighting the ensemble average
 *   	
 * 
 * @author taitan
 *
 */
public class SimOverlapMultiHarmonic extends Simulation{

    private static final long serialVersionUID = 1L;
    protected SpeciesSpheresMono species;
    protected Box boxA, boxB;
    protected IController controller;
    protected ActivityIntegrate activityIntegrate;
    protected IntegratorMC integratorA, integratorB;
    protected P1Harmonic potentialA, potentialB;
    protected PotentialMaster potentialMasterA, potentialMasterB;
    protected IntegratorOverlap integratorOverlap;
    protected DataSourceVirialOverlap dsvo;
    protected IntegratorBox[] integrators;
    protected AccumulatorVirialOverlapSingleAverage[] accumulators;
    protected IDataSource[] meters;
    protected DataPump[] accumulatorPumps;
    protected double refPref;

    public SimOverlapMultiHarmonic(int numAtoms, double wA, double wB, double temperature){
		super(Space1D.getInstance());

		refPref =(wA/wB); // the analytic solution for free energy difference

		potentialMasterA = new PotentialMasterMonatomic(this);
		potentialMasterB = new PotentialMasterMonatomic(this);

		integrators = new IntegratorBox[2];
		accumulatorPumps = new DataPump[2];
        meters = new IDataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];

		species = new SpeciesSpheresMono(this, space);
		addSpecies(species);

		//System A
		boxA = new Box(new BoundaryRectangularNonperiodic(space), space);
		addBox(boxA);
		boxA.getBoundary().setBoxSize(new Vector1D(3.0));
		boxA.setNMolecules(species, numAtoms);

		integratorA = new IntegratorMC(this, potentialMasterA);
		integratorA.setBox(boxA);
		integratorA.setTemperature(temperature);
		integrators[0] = integratorA;

		potentialA = new P1Harmonic(space);
		potentialA.setSpringConstant(wA);
		integratorA.getMoveManager().addMCMove(new MCMoveMultiHarmonic(potentialA, random));
        potentialMasterA.addPotential(potentialA, new AtomType[]{species.getLeafType()});

		MeterPotentialEnergy meterPEA = new MeterPotentialEnergy(potentialMasterA);
		meterPEA.setBox(boxA);

		//System B
		boxB = new Box(new BoundaryRectangularNonperiodic(space), space);
		addBox(boxB);
		boxB.setNMolecules(species, numAtoms);

		integratorB = new IntegratorMC(this, potentialMasterB);
		integratorB.setBox(boxB);
		integratorB.setTemperature(temperature);
		integrators[1] = integratorB;

		potentialB = new P1Harmonic(space);
		potentialB.setSpringConstant(wB);
		potentialB.setX0(new Vector1D(1.0));
		integratorB.getMoveManager().addMCMove(new MCMoveMultiHarmonic(potentialB, random));
        potentialMasterB.addPotential(potentialB, new AtomType[]{species.getLeafType()});

		MeterPotentialEnergy meterPEB = new MeterPotentialEnergy(potentialMasterB);
		meterPEB.setBox(boxB);

		//Overlap
		integratorOverlap = new IntegratorOverlap(new IntegratorBox[] {integratorA, integratorB});
		MeterBoltzmannA meterA = new MeterBoltzmannA(integratorA, potentialMasterB);
		meterA.setTemperature(temperature);
		meters[0] = meterA;
		setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, true), 0);

		MeterBoltzmannB meterB = new MeterBoltzmannB(integratorB, potentialMasterA);
		meterB.setTemperature(temperature);
		meters[1] = meterB;
		setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, false), 1);

		activityIntegrate = new ActivityIntegrate(integratorOverlap);
		getController().addAction(activityIntegrate);

	}

	public static void main(String[] args){

		SimOverlapParam params = new SimOverlapParam();
		String inputFilename = null;
		if (args.length > 0){
			inputFilename = args[0];
		}
		if (inputFilename != null){
			ReadParameters readParameters = new ReadParameters(inputFilename);
			readParameters.readParameters();
		}

		int numAtoms = params.numAtoms;
		double wA = params.wA;
		double wB = params.wB;
		long numSteps = params.numSteps;
		String filename = params.filename;
		if (filename.length() == 0){
			System.err.println("Need Input file");
			filename = "Mharm_Default";
		}
		double temperature = 1.0;

        System.out.println("Running 1D uncoupled harmonic oscillator simulation");
		System.out.println(numAtoms + " atoms at temperature " + temperature);
		System.out.println((numSteps/1000) + " total steps of 1000");

        SimOverlapMultiHarmonic sim = new SimOverlapMultiHarmonic(numAtoms, wA, wB, temperature);
		sim.integratorOverlap.setNumSubSteps(1000);
		numSteps /= 1000;

		/*
		 * Equilibration Step
		 */
		sim.setRefPref(sim.refPref);
		sim.activityIntegrate.setMaxSteps(numSteps/10);
		for (int i=0; i<2; i++) {
			if (sim.integrators[i] instanceof IntegratorMC) ((IntegratorMC)sim.integrators[i]).getMoveManager().setEquilibrating(true);
		}
	    sim.getController().actionPerformed();
	    sim.getController().reset();

        System.out.println("Equilibration finished");
	    System.out.flush();
	    /*
         *
	     */

        // System A to System B Sampling
        MeterWorkABandReweighting meterWorkA =
                new MeterWorkABandReweighting(sim.integratorA, sim.potentialMasterB, sim.refPref);
        meterWorkA.setTemperature(temperature);

        DataFork dataForkA = new DataFork();
        DataPumpListener dataPumpA = new DataPumpListener(meterWorkA, dataForkA,1);

        AccumulatorAverageFixed dataAverageA = new AccumulatorAverageFixed(1);
        AccumulatorHistogram histogramAB = new AccumulatorHistogram(new HistogramSimple(4000, new DoubleRange(-200, 200)));
        dataForkA.addDataSink(dataAverageA);
        dataForkA.addDataSink(histogramAB);

        sim.integrators[0].getEventManager().addListener(dataPumpA);

        // System B to System A Sampling
        MeterWorkBAandReweighting meterWorkB =
                new MeterWorkBAandReweighting(sim.integratorB, sim.potentialMasterA, sim.refPref);
        meterWorkB.setTemperature(temperature);

        DataFork dataForkB = new DataFork();
        DataPumpListener dataPumpB = new DataPumpListener(meterWorkB, dataForkB, 1);

        AccumulatorAverageFixed dataAverageB = new AccumulatorAverageFixed(1);
        AccumulatorHistogram histogramBA = new AccumulatorHistogram(new HistogramSimple(4000, new DoubleRange(-200, 200)));
        dataForkB.addDataSink(dataAverageB);
        dataForkB.addDataSink(histogramBA);

        sim.integrators[1].getEventManager().addListener(dataPumpB);

        sim.activityIntegrate.setMaxSteps(numSteps);
	    sim.getController().actionPerformed();

        System.out.println(" ");
	    System.out.println("final reference optimal step frequency "+sim.integratorOverlap.getIdealRefStepFraction()
        		+" (actual: "+sim.integratorOverlap.getRefStepFraction()+")");
        double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
        double ratio = ratioAndError[0];
        double error = ratioAndError[1];
        double betaFAB = -Math.log(ratio);

        System.out.println("\nratio average: "+ratio+" ,error: "+error);
        System.out.println("free energy difference: "+(-temperature*Math.log(ratio))+" ,error: "+temperature*(error/ratio));

        DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
        double betaFAW = -Math.log(((DataDoubleArray) allYourBase.getData(AccumulatorAverage.AVERAGE.index)).getData()[1]);
        System.out.println("System-A ratio average: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.AVERAGE.index)).getData()[1]
                + " stdev: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[1]
                + " error: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.ERROR.index)).getData()[1]);

        allYourBase = (DataGroup)sim.accumulators[1].getData(sim.dsvo.minDiffLocation());
        double betaFBW = -Math.log(((DataDoubleArray) allYourBase.getData(AccumulatorAverage.AVERAGE.index)).getData()[1]);
        System.out.println("System-B ratio average: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.AVERAGE.index)).getData()[1]
                + " stdev: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[1]
                + " error: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.ERROR.index)).getData()[1]);

        /*
         * Refer Wu & Kofke JCp 123,054103(2003) Eq (6)
         */
        double betaUAB = dataAverageA.getData().getValue(AccumulatorAverage.AVERAGE.index);
        double betaUBA = dataAverageB.getData().getValue(AccumulatorAverage.AVERAGE.index);

        double SAB = betaUAB - betaFAB;
        double SBA = betaUBA + betaFAB;

        double betaUAWf = meterWorkA.getBetaUAWf();  // < beta*U_AW>A
        double betaUAWr = meterWorkA.getBetaUAWr();  // < beta*U_AW>W

        double betaUBWf = meterWorkB.getBetaUBWf();	// < beta*U_BW>B
        double betaUBWr = meterWorkB.getBetaUBWr();	// < beta*U_BW>W

        double SAW =  betaUAWf - betaFAW;
        double SWA = - betaUAWr + betaFAW;

        double SBW = betaUBWf - betaFBW;
        double SWB = - betaUBWr + betaFBW;

        System.out.println("");
        System.out.println("SAB: "+ SAB + " ; betaUAB: " + betaUAB + " ;betaFAB: " + betaFAB);
        System.out.println("SBA: "+ SBA + " ; betaUBA: " + betaUBA + " ;betaFAB: " + betaFAB);

        System.out.println("");
        System.out.println("SAW: "+ SAW + " ; betaUAWf: " + betaUAWf + " ;betaFAW: " + betaFAW);
        System.out.println("SWA: "+ SWA + " ; betaUAWr: " + betaUAWr + " ;betaFAW: " + betaFAW);

        System.out.println("");
        System.out.println("SBW: "+ SBW + " ; betaUBWf: " + betaUBWf + " ;betaFBW: " + betaFBW);
        System.out.println("SWB: "+ SWB + " ; betaUBWr: " + betaUBWr + " ;betaFBW: " + betaFBW);

        /*
         * System A -----> B
         */
        // AB
        DataLogger dataLogger = new DataLogger();
        DataTableWriter dataTableWriter = new DataTableWriter();
        dataLogger.setFileName(filename + "_hist_AB");
        dataTableWriter.setIncludeHeader(false);

        dataLogger.setDataSink(dataTableWriter);
        dataLogger.putDataInfo(histogramAB.getDataInfo());
        dataLogger.putData(histogramAB.getData());
        dataLogger.closeFile();

        // AW
        dataLogger = new DataLogger();
        dataTableWriter = new DataTableWriter();
        dataLogger.setFileName(filename + "_hist_AW");
        dataTableWriter.setIncludeHeader(false);

        dataLogger.setDataSink(dataTableWriter);
        dataLogger.putDataInfo(meterWorkA.getDataInfoHistogramBetaUAWf());
        dataLogger.putData(meterWorkA.getDataHistogramBetaUAWf());
        dataLogger.closeFile();

        // WA
        dataLogger = new DataLogger();
        dataTableWriter = new DataTableWriter();
        dataLogger.setFileName(filename + "_hist_WA");
        dataTableWriter.setIncludeHeader(false);

        dataLogger.setDataSink(dataTableWriter);
        dataLogger.putDataInfo(meterWorkA.getDataInfoHistogramBetaUAWr());
        dataLogger.putData(meterWorkA.getDataHistogramBetaUAWr());
        dataLogger.closeFile();

        /*
         * System B -----> A
         */
        //BA
        dataLogger = new DataLogger();
        dataTableWriter = new DataTableWriter();
        dataLogger.setFileName(filename + "_hist_BA");
        dataTableWriter.setIncludeHeader(false);

        dataLogger.setDataSink(dataTableWriter);
        dataLogger.putDataInfo(histogramBA.getDataInfo());
        dataLogger.putData(histogramBA.getData());
        dataLogger.closeFile();

        //BW
        dataLogger = new DataLogger();
        dataTableWriter = new DataTableWriter();
        dataLogger.setFileName(filename + "_hist_BA");
        dataTableWriter.setIncludeHeader(false);

        dataLogger.setDataSink(dataTableWriter);
        dataLogger.putDataInfo(meterWorkB.getDataInfoHistogramBetaUBWf());
        dataLogger.putData(meterWorkB.getDataHistogramBetaUBWf());
        dataLogger.closeFile();

        // WB
        dataLogger = new DataLogger();
        dataTableWriter = new DataTableWriter();
        dataLogger.setFileName(filename + "_hist_WB");
        dataTableWriter.setIncludeHeader(false);

        dataLogger.setDataSink(dataTableWriter);
        dataLogger.putDataInfo(meterWorkB.getDataInfoHistogramBetaUBWr());
        dataLogger.putData(meterWorkB.getDataHistogramBetaUBWr());
        dataLogger.closeFile();


    }

    public void setAccumulator(AccumulatorVirialOverlapSingleAverage newAccumulator, int iBox) {

        accumulators[iBox] = newAccumulator;

        newAccumulator.setBlockSize(1); // setting the block size = 300

        if (accumulatorPumps[iBox] == null) {
            accumulatorPumps[iBox] = new DataPump(meters[iBox], newAccumulator);
            IntegratorListenerAction pumpListener = new IntegratorListenerAction(accumulatorPumps[iBox]);
            integrators[iBox].getEventManager().addListener(pumpListener);
            if (iBox == 1) {
                pumpListener.setInterval(boxB.getMoleculeList().getMoleculeCount());
            }
        } else {
            accumulatorPumps[iBox].setDataSink(newAccumulator);
        }
        if (integratorOverlap != null && accumulators[0] != null && accumulators[1] != null) {
            dsvo = new DataSourceVirialOverlap(accumulators[0], accumulators[1]);
            integratorOverlap.setReferenceFracSource(dsvo);
        }
    }

    public void setRefPref(double refPrefCenter, double span) {
        refPref = refPrefCenter;
        accumulators[0].setBennetParam(refPrefCenter, span);
        accumulators[1].setBennetParam(refPrefCenter, span);

    }

    public void setRefPref(double newRefPref) {
        System.out.println("setting ref pref to " + newRefPref);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1, true), 0);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1, false), 1);
        setRefPref(newRefPref, 1);
    }

    public static class SimOverlapParam extends ParameterBase {
		public int numAtoms = 10;
		public double wA = 1.0;
		public double wB = 0.8;
		public long numSteps = 1000000;
		public String filename = "Mharm_01_01_01";
	}
	
}
