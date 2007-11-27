package etomica.normalmode;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.DataPump;
import etomica.data.DataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.math.SpecialFunctions;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P1HardPeriodic;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
import etomica.potential.Potential2HardSpherical;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.Species;
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

    public SimOverlap(Space space, int numAtoms, double density, double temperature, String filename, double harmonicFudge) {
        super(space, true);

        integrators = new IntegratorBox[2];
        accumulatorPumps = new DataPump[2];
        meters = new DataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];

        SpeciesSpheresMono species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);

        // TARGET
        PotentialMasterList potentialMasterTarget = new PotentialMasterList(this);
        boxTarget = new Box(this);
        addBox(boxTarget);
        boxTarget.setNMolecules(species, numAtoms);

        IntegratorHard integratorTarget = new IntegratorHard(potentialMasterTarget, getRandom(), 4, 1.0);

        integratorTarget.setIsothermal(false);
        integrators[1] = integratorTarget;

        Potential2 p2 = new P2HardSphere(space, 1.0, true);
        if (space.D() == 1) {
            // don't need this for the target system, but we do need it for the reference
            p2 = new P2XOrder(space, (Potential2HardSpherical)p2);
        }
        potentialMasterTarget.addPotential(p2, new Species[]{species,species});

        if (space.D() == 1) {
            primitive = new PrimitiveCubic(space, 1.0/density);
            boundaryTarget = new BoundaryRectangularPeriodic(space, getRandom(), numAtoms/density);
            integratorTarget.setNullPotential(new P1HardPeriodic(space));
            nCells = new int[]{numAtoms};
            basis = new BasisMonatomic(space);
        } else {
            double L = Math.pow(4.0/density, 1.0/3.0);
            primitive = new PrimitiveCubic(space, L);
            int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
            nCells = new int[]{n,n,n};
            boundaryTarget = new BoundaryRectangularPeriodic(space, random, n * L);
            basis = new BasisCubicFcc();
        }
        boxTarget.setBoundary(boundaryTarget);

        CoordinateDefinitionLeaf coordinateDefinitionTarget = new CoordinateDefinitionLeaf(boxTarget, primitive, basis);
        coordinateDefinitionTarget.initializeCoordinates(nCells);

        if (potentialMasterTarget instanceof PotentialMasterList) {
            double neighborRange;
            if (space.D() == 1) {
                neighborRange = 1.01 / density;
            }
            else {
                //FCC
                double L = Math.pow(4.01/density, 1.0/3.0);
                neighborRange = L / Math.sqrt(2.0);
            }
            ((PotentialMasterList)potentialMasterTarget).setRange(neighborRange);
            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
            ((PotentialMasterList)potentialMasterTarget).getNeighborManager(boxTarget).reset();
        }

        integratorTarget.setBox(boxTarget);
    
        // HARMONIC
        if (space.D() == 1) {
            boundaryHarmonic = new BoundaryRectangularPeriodic(space, getRandom(), numAtoms/density);
        } else {
            double L = Math.pow(4.0/density, 1.0/3.0);
            int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
            boundaryHarmonic = new BoundaryRectangularPeriodic(space, random, n * L);
        }
        boxHarmonic = new Box(boundaryHarmonic);
        addBox(boxHarmonic);
        boxHarmonic.setNMolecules(species, numAtoms);

        IntegratorMC integratorHarmonic = new IntegratorMC(null, random, 1.0);
        integratorHarmonic.setBox(boxHarmonic);

        MCMoveHarmonic move = new MCMoveHarmonic(getRandom());
        integratorHarmonic.getMoveManager().addMCMove(move);
        integrators[0] = integratorHarmonic;
        
        CoordinateDefinitionLeaf coordinateDefinitionHarmonic = new CoordinateDefinitionLeaf(boxHarmonic, primitive, basis);
        coordinateDefinitionHarmonic.initializeCoordinates(nCells);

        if(space.D() == 1) {
            normalModes = new NormalModes1DHR();
        } else {
            normalModes = new NormalModesFromFile(filename, space.D());
        }
        normalModes.setHarmonicFudge(harmonicFudge);
        normalModes.setTemperature(temperature);
        
        WaveVectorFactory waveVectorFactory = normalModes.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(boxHarmonic);
        move.setOmegaSquared(normalModes.getOmegaSquared(boxHarmonic), waveVectorFactory.getCoefficients());
        move.setEigenVectors(normalModes.getEigenvectors(boxHarmonic));
        move.setWaveVectors(waveVectorFactory.getWaveVectors());
        move.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
        move.setCoordinateDefinition(coordinateDefinitionHarmonic);
        move.setTemperature(temperature);
        
        move.setBox(boxHarmonic);
        
        integratorHarmonic.setBox(boxHarmonic);
        
        if (potentialMasterTarget instanceof PotentialMasterList) {
            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
            ((PotentialMasterList)potentialMasterTarget).getNeighborManager(boxHarmonic).reset();
        }

        // OVERLAP
        integratorOverlap = new IntegratorOverlap(random, new IntegratorBox[]{integratorHarmonic, integratorTarget});
        MeterHarmonicEnergy meterHarmonicEnergy = new MeterHarmonicEnergy(coordinateDefinitionTarget, normalModes);
        meterHarmonicEnergy.setBox(boxTarget);
        MeterBoltzmannTarget meterTarget = new MeterBoltzmannTarget(integratorTarget, meterHarmonicEnergy);
        // lattice energy is 0
        meters[1] = meterTarget;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, false), 1);

        MeterBoltzmannHarmonic meterHarmonic = new MeterBoltzmannHarmonic(move, potentialMasterTarget);
        meterHarmonic.setTemperature(temperature);
        meters[0] = meterHarmonic;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, true), 0);
        
        setRefPref(1.0, 30);

        activityIntegrate = new ActivityIntegrate(integratorOverlap);
        
        getController().addAction(activityIntegrate);
    }

    public void setRefPref(double refPrefCenter, double span) {
        refPref = refPrefCenter;
        accumulators[0].setBennetParam(refPrefCenter,span);
        accumulators[1].setBennetParam(refPrefCenter,span);
    }

    public void setAccumulator(AccumulatorVirialOverlapSingleAverage newAccumulator, int iBox) {
        accumulators[iBox] = newAccumulator;
        if (accumulatorPumps[iBox] == null) {
            accumulatorPumps[iBox] = new DataPump(meters[iBox],newAccumulator);
            integrators[iBox].addIntervalAction(accumulatorPumps[iBox]);
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
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,true),0);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,false),1);
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
                setRefPref(refPref,1);
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

            setAccumulator(new AccumulatorVirialOverlapSingleAverage(41,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(41,false),1);
            setRefPref(1e40,40);
            activityIntegrate.setMaxSteps(initSteps);
            getController().actionPerformed();
            getController().reset();

            int newMinDiffLoc = dsvo.minDiffLocation();
            refPref = accumulators[0].getBennetAverage(newMinDiffLoc)
                /accumulators[1].getBennetAverage(newMinDiffLoc);
            if (Double.isNaN(refPref) || refPref == 0 || Double.isInfinite(refPref)) {
                throw new RuntimeException("Simulation failed to find a valid ref pref");
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
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,false),1);
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
        }
        if (inputFilename != null) {
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
            filename = "normal_modes_HS_"+D+"D_"+numMolecules;
        }
        String refFileName = args.length > 0 ? filename+"_ref" : null;

        System.out.println("Running "+(D==1 ? "1D" : (D==3 ? "FCC" : "2D hexagonal")) +" hard sphere overlap simulation");
        System.out.println(numMolecules+" atoms at density "+density);
        System.out.println("harmonic fudge: "+harmonicFudge);
        System.out.println((numSteps/1000)+" total steps of 1000");
        System.out.println("output data to "+filename);

        //instantiate simulation
        SimOverlap sim = new SimOverlap(Space.getInstance(D), numMolecules, density, temperature, filename, harmonicFudge);
        
        //start simulation
        sim.integratorOverlap.setNumSubSteps(1000);
        numSteps /= 1000;

//        StopWatcher stopWatcher = new StopWatcher("stop", sim, filename);
//        IntervalActionAdapter iaa = new IntervalActionAdapter(stopWatcher);
//        iaa.setActionInterval(100);
//        sim.integratorOverlap.addListener(iaa);

        sim.initRefPref(refFileName, numSteps/20);
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
        sim.getController().actionPerformed();

        System.out.println("final reference optimal step frequency "+sim.integratorOverlap.getStepFreq0()+" (actual: "+sim.integratorOverlap.getActualStepFreq0()+")");
        
        double[][] omega2 = sim.normalModes.getOmegaSquared(sim.boxTarget);
        double[] coeffs = sim.normalModes.getWaveVectorFactory().getCoefficients();
        double AHarmonic = 0;
        for(int i=0; i<omega2.length; i++) {
            for(int j=0; j<omega2[0].length; j++) {
                if (!Double.isInfinite(omega2[i][j])) {
                    AHarmonic += coeffs[i]*Math.log(omega2[i][j]*coeffs[i]/(temperature*Math.PI));//temperature should be here
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

        double ratio = sim.dsvo.getDataAsScalar();
        double error = sim.dsvo.getError();
        System.out.println("ratio average: "+ratio+", error: "+error);
        System.out.println("free energy difference: "+(-Math.log(ratio))+", error: "+(error/ratio));
        System.out.println("target free energy: "+(AHarmonic-Math.log(ratio)));
        DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
        System.out.println("harmonic ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1]);
        
        allYourBase = (DataGroup)sim.accumulators[1].getData(sim.accumulators[1].getNBennetPoints()-sim.dsvo.minDiffLocation()-1);
        System.out.println("target ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1]);

        if(D==1) {
            double AHR = -(numMolecules-1)*Math.log(numMolecules/density-numMolecules) + SpecialFunctions.lnFactorial(numMolecules) ;
            System.out.println("Hard-rod free energy: "+AHR);
        }
    }

    private static final long serialVersionUID = 1L;
    public CoordinateDefinition coordinateDefinition;
    public IntegratorOverlap integratorOverlap;
    public DataSourceVirialOverlap dsvo;
    public IntegratorBox[] integrators;
    public ActivityIntegrate activityIntegrate;
    public Box boxTarget, boxHarmonic;
    public Boundary boundaryTarget, boundaryHarmonic;
    public int[] nCells;
    public Basis basis;
    public NormalModes normalModes;
    public Primitive primitive;
    public double refPref;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public DataPump[] accumulatorPumps;
    public DataSource[] meters;

    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules = 32;
        public double density = 1.2;
        public int D = 3;
        public long numSteps = 100000;
        public double harmonicFudge = .1;
        public String filename = "S_d120_N32_s1000";
        public double temperature = 1.0;
    }
}
