package etomica.normalmode;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.DataPump;
import etomica.data.DataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorPhase;
import etomica.integrator.IntervalActionAdapter;
import etomica.lattice.BravaisLattice;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.lattice.crystal.PrimitiveFcc;
import etomica.math.SpecialFunctions;
import etomica.nbr.list.PotentialMasterList;
import etomica.phase.Phase;
import etomica.potential.P1HardPeriodic;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
import etomica.potential.Potential2HardSpherical;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;
import etomica.virial.overlap.IntegratorOverlap;

/**
 * Simulation to run sampling with the hard sphere potential, but measuring
 * the harmonic potential based on normal mode data from a previous simulation.
 * 
 * @author Andrew Schultz
 */
public class SimOverlap extends Simulation {

    public SimOverlap(Space space, int numAtoms, double density, String filename, double harmonicFudge) {
        super(space, true, (space.D() == 1 ? new PotentialMasterList(space) : new PotentialMaster(space)));

        defaults.makeLJDefaults();
        defaults.atomSize = 1.0;
        integrators = new IntegratorPhase[2];
        accumulatorAAs = new IntervalActionAdapter[2];
        accumulatorPumps = new DataPump[2];
        meters = new DataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];

        SpeciesSpheresMono species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);

        // TARGET
        
        phaseTarget = new Phase(this);
        addPhase(phaseTarget);
        phaseTarget.getAgent(species).setNMolecules(numAtoms);

        IntegratorHard integratorTarget = new IntegratorHard(potentialMaster, getRandom(), 4, 1.0);

        integratorTarget.setIsothermal(false);
        integrators[1] = integratorTarget;

        Potential2 p2 = new P2HardSphere(space, 1.0, true);
        if (space.D() == 1) {
            p2 = new P2XOrder(space, (Potential2HardSpherical)p2);
        }
        potentialMaster.addPotential(p2, new AtomType[]{species.getMoleculeType(),species.getMoleculeType()});

        Primitive primitive;
        if (space.D() == 1) {
            primitive = new PrimitiveCubic(space, 1.0/density);
            boundaryTarget = new BoundaryRectangularPeriodic(space, getRandom(), numAtoms/density);
            integratorTarget.setNullPotential(new P1HardPeriodic(space));
        } else {
            primitive = new PrimitiveFcc(space, 1);
            double v = primitive.unitCell().getVolume();
            primitive.scaleSize(Math.pow(v*density,-1.0/3.0));
            int n = (int)Math.round(Math.pow(numAtoms, 1.0/3.0));
            boundaryTarget = new BoundaryDeformableLattice(primitive, getRandom(), new int[]{n,n,n});
        }
        phaseTarget.setBoundary(boundaryTarget);

        lattice = new BravaisLattice(primitive);
        ConfigurationLattice config = new ConfigurationLattice(lattice);

        config.initializeCoordinates(phaseTarget);

        if (potentialMaster instanceof PotentialMasterList) {
            double neighborRange;
            if (space.D() == 1) {
                neighborRange = 1.01 / density;
            }
            else {
                //FCC
                double L = Math.pow(4.01/density, 1.0/3.0);
                neighborRange = L / Math.sqrt(2.0);
            }
            ((PotentialMasterList)potentialMaster).setRange(neighborRange);
            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
            ((PotentialMasterList)potentialMaster).getNeighborManager(phaseTarget).reset();
        }

        integratorTarget.setPhase(phaseTarget);

    
        // HARMONIC
        boundaryHarmonic = new BoundaryRectangularPeriodic(this);
        phaseHarmonic = new Phase(boundaryHarmonic);
        addPhase(phaseHarmonic);
        phaseHarmonic.getAgent(species).setNMolecules(numAtoms);

        IntegratorMC integratorHarmonic = new IntegratorMC(potentialMaster, random, 1.0);

        MCMoveHarmonic move = new MCMoveHarmonic(null, getRandom());
        integratorHarmonic.getMoveManager().addMCMove(move);
        integrators[0] = integratorHarmonic;
        
        if (space.D() == 1) {
            boundaryHarmonic = new BoundaryRectangularPeriodic(space, getRandom(), numAtoms/density);
        } else {
            int n = (int)Math.round(Math.pow(numAtoms, 1.0/3.0));
            boundaryHarmonic = new BoundaryDeformableLattice(primitive, getRandom(), new int[]{n,n,n});
        }
        phaseHarmonic.setBoundary(boundaryHarmonic);

        config.initializeCoordinates(phaseHarmonic);
        
        if(space.D() == 1) {
            normalModes = new NormalModes1DHR();
        } else {
            normalModes = new NormalModesFromFile(filename, space.D());
        }
        normalModes.setHarmonicFudge(harmonicFudge);
        
        WaveVectorFactory waveVectorFactory = normalModes.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(phaseHarmonic);
        move.setOmegaSquared(normalModes.getOmegaSquared(phaseHarmonic), waveVectorFactory.getCoefficients());
        move.setEigenVectors(normalModes.getEigenvectors(phaseHarmonic));
        move.setWaveVectors(waveVectorFactory.getWaveVectors());
        move.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
        move.setCoordinateDefinition(new CoordinateDefinitionLeaf(space));
        
        move.setPhase(phaseHarmonic);
        
        integratorHarmonic.setPhase(phaseHarmonic);
        
        if (potentialMaster instanceof PotentialMasterList) {
            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
            ((PotentialMasterList)potentialMaster).getNeighborManager(phaseHarmonic).reset();
        }

        
        // OVERLAP
        MeterHarmonicEnergy meterHarmonicEnergy = new MeterHarmonicEnergy(new CoordinateDefinitionLeaf(space), normalModes);
        meterHarmonicEnergy.setPhase(phaseTarget);
        MeterBoltzmannTarget meterTarget = new MeterBoltzmannTarget(integratorTarget, meterHarmonicEnergy);
        meters[1] = meterTarget;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, false), 1);

        MeterBoltzmannHarmonic meterHarmonic = new MeterBoltzmannHarmonic(move, potentialMaster);
        meterHarmonic.setTemperature(1.0);
        meters[0] = meterHarmonic;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, true), 0);
        
        setRefPref(1.0, 30);
        
        integratorOverlap = new IntegratorOverlap(null, random, new IntegratorPhase[]{integratorHarmonic, integratorTarget},
                accumulators);
        integratorOverlap.setDSVO(dsvo);
        // sadly, we have to ignore overlap in both phases since we expect "overlaps" in the harmonic phase.
        activityIntegrate = new ActivityIntegrate(integratorOverlap, false, true);
        
        getController().addAction(activityIntegrate);
    }

    public void setRefPref(double refPrefCenter, double span) {
        refPref = refPrefCenter;
        accumulators[0].setBennetParam(refPrefCenter,span);
        accumulators[1].setBennetParam(refPrefCenter,span);
    }

    public void setAccumulator(AccumulatorVirialOverlapSingleAverage newAccumulator, int iPhase) {
        accumulators[iPhase] = newAccumulator;
        if (accumulatorPumps[iPhase] == null) {
            accumulatorPumps[iPhase] = new DataPump(meters[iPhase],newAccumulator);
            accumulatorAAs[iPhase] = new IntervalActionAdapter(accumulatorPumps[iPhase]);
            integrators[iPhase].addListener(accumulatorAAs[iPhase]);
        }
        else {
            accumulatorPumps[iPhase].setDataSink(newAccumulator);
        }
        accumulatorAAs[iPhase].setActionInterval(1);
        if (integratorOverlap != null) {
            dsvo = new DataSourceVirialOverlap(accumulators[0],accumulators[1]);
            integratorOverlap.setDSVO(dsvo);
        }
    }
    
    public void setRefPref(double newRefPref) {
        System.out.println("setting ref pref (explicitly) to "+newRefPref);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(this,1,true),0);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(this,1,false),1);
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
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(this,1,true),0);
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(this,1,false),1);
                setRefPref(refPref,1);
            }
            catch (IOException e) {
                // file not there, which is ok.
            }
        }
        
        if (refPref == -1) {
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(this,41,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(this,41,false),1);
            setRefPref(1e20,100);
            activityIntegrate.setMaxSteps(initSteps);
            getController().run();
            getController().reset();

            int newMinDiffLoc = dsvo.minDiffLocation();
            refPref = accumulators[0].getBennetAverage(newMinDiffLoc)
                /accumulators[1].getBennetAverage(newMinDiffLoc);
            if (Double.isNaN(refPref) || refPref == 0 || Double.isInfinite(refPref)) {
                throw new RuntimeException("Simulation failed to find a valid ref pref");
            }
            System.out.println("setting ref pref to "+refPref);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(this,21,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(this,21,false),1);
            setRefPref(refPref,10);
            activityIntegrate.setMaxSteps(initSteps);
            getController().run();

            newMinDiffLoc = dsvo.minDiffLocation();
            refPref = accumulators[0].getBennetAverage(newMinDiffLoc)
                /accumulators[1].getBennetAverage(newMinDiffLoc);
            if (Double.isNaN(refPref) || refPref == 0 || Double.isInfinite(refPref)) {
                throw new RuntimeException("Simulation failed to find a valid ref pref");
            }
            System.out.println("setting ref pref to "+refPref);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(this,11,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(this,11,false),1);
            setRefPref(refPref,0.2);
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
        getController().run();
        getController().reset();

        if (refPref == -1) {
            int newMinDiffLoc = dsvo.minDiffLocation();
            refPref = accumulators[0].getBennetAverage(newMinDiffLoc)
                /accumulators[1].getBennetAverage(newMinDiffLoc);
            System.out.println("setting ref pref to "+refPref+" ("+newMinDiffLoc+")");
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(this,1,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(this,1,false),1);
            setRefPref(refPref,1);
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
     * @see SimOverlap.SimOverlapParam
     */
    public static void main(String[] args) {
        
        //set up simulation parameters
        SimOverlapParam params = new SimOverlapParam();
        String inputFilename = null;
        if (args.length > 0) {
            inputFilename = args[0];
            ReadParameters readParameters = new ReadParameters(inputFilename, params);
            readParameters.readParameters();
        }
        double density = params.density;
        long numSteps = params.numSteps;
        int numMolecules = params.numMolecules;
        double harmonicFudge = params.harmonicFudge;
        double temperature = params.temperature;
        int D = params.D;
        String filename = params.filename;
        if (filename.length() == 0) {
            filename = "normal_modes3D";
        }
        String refFileName = args.length > 0 ? filename+"_ref" : null;

        System.out.println("Running "+(D==1 ? "1D" : (D==3 ? "FCC" : "2D hexagonal")) +" hard sphere overlap simulation");
        System.out.println(numMolecules+" atoms at density "+density);
        System.out.println("harmonic fudge: "+harmonicFudge);
        System.out.println((numSteps/1000)+" total steps of 1000");
        System.out.println("output data to "+filename);

        //instantiate simulation
        SimOverlap sim = new SimOverlap(Space.getInstance(D), numMolecules, density, filename, harmonicFudge);
        
        //start simulation
        sim.integratorOverlap.setNumSubSteps(1000);
        numSteps /= 1000;
        
        
        sim.initRefPref(refFileName, numSteps/40);
        if (Double.isNaN(sim.refPref) || sim.refPref == 0 || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("Simulation failed to find a valid ref pref");
        }
        
        sim.equilibrate(refFileName, numSteps/10);
        if (Double.isNaN(sim.refPref) || sim.refPref == 0 || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("Simulation failed to find a valid ref pref");
        }
        
        System.out.println("equilibration finished");

        sim.integratorOverlap.getMoveManager().setEquilibrating(false);
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().run();

        System.out.println("final reference step frequency "+sim.integratorOverlap.getStepFreq0());
        
        double[][] omega2 = sim.normalModes.getOmegaSquared(sim.phaseTarget);
        double[] coeffs = sim.normalModes.getWaveVectorFactory().getCoefficients();
        double AHarmonic = 0;
        for(int i=0; i<omega2.length; i++) {
            for(int j=0; j<omega2[0].length; j++) {
                AHarmonic += coeffs[i]*Math.log(omega2[i][j]*coeffs[i]/(temperature*Math.PI));//coeffs in log?
            }
        }
        if (numMolecules % 2 == 0) {
            if (D == 1) {
                AHarmonic += Math.log(((D*numMolecules - 2)/2.0) / Math.pow(numMolecules,0.5*D));
            }
            else if (D == 3) {
                AHarmonic += Math.log(((D*numMolecules + 3)/2.0) / Math.pow(numMolecules,0.5*D));
            }
        }
        else {
            if (D == 1) {
                AHarmonic += Math.log(((D*numMolecules - 1)/2.0) / Math.pow(numMolecules,0.5*D));
            }
            else if (D == 3) {
                AHarmonic += Math.log(((D*numMolecules - 18)/2.0) / Math.pow(numMolecules,0.5*D));
            }
        }
        System.out.println("Harmonic-reference free energy: "+AHarmonic*temperature);

        double ratio = sim.dsvo.getDataAsScalar();
        double error = sim.dsvo.getError();
        System.out.println("ratio average: "+ratio+", error: "+error);
        System.out.println("free energy difference: "+(-Math.log(ratio))+", error: "+(error/ratio));
        System.out.println("target free energy: "+(AHarmonic-Math.log(ratio)));
        DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
        System.out.println("harmonic ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1]);
        System.out.println("harmonic   average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[0]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[0]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[0]);
        System.out.println("harmonic overlap average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1]);
        
        allYourBase = (DataGroup)sim.accumulators[1].getData(sim.accumulators[1].getNBennetPoints()-sim.dsvo.minDiffLocation()-1);
        System.out.println("target ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1]);
        System.out.println("target average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[0]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[0]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[0]);
        System.out.println("target overlap average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1]);

        if(D==1) {
            double AHR = -(numMolecules-1)*Math.log(numMolecules/density-numMolecules) + SpecialFunctions.lnFactorial(numMolecules) ;
            System.out.println("Hard-rod free energy: "+AHR);
        }
    }

    private static final long serialVersionUID = 1L;
    public IntegratorOverlap integratorOverlap;
    public DataSourceVirialOverlap dsvo;
    public IntegratorPhase[] integrators;
    public ActivityIntegrate activityIntegrate;
    public Phase phaseTarget, phaseHarmonic;
    public Boundary boundaryTarget, boundaryHarmonic;
    public BravaisLattice lattice;
    public NormalModes normalModes;
    public double refPref;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public DataPump[] accumulatorPumps;
    public IntervalActionAdapter[] accumulatorAAs;
    public DataSource[] meters;

    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules = 27;
        public double density = 1.04;
        public int D = 3;
        public long numSteps = 100000;
        public double harmonicFudge = .5;
        public String filename = "";
        public double temperature = 1.0;
    }
}
