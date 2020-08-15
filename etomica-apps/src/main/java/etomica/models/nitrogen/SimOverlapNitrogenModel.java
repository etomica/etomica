/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;


import etomica.action.activity.ActivityIntegrate2;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.data.IDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.normalmode.*;
import etomica.overlap.IntegratorOverlap;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.Kelvin;
import etomica.units.Pascal;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;

import java.io.*;

/**
 * Simulation to run sampling with the nitrogen model to find the free energy difference
 * 	between the harmonic reference system and the target system
 * 
 * The Bennett's Overlapping Sampling Simulation
 * 	- used to check for the computation time
 * 
 * @author Tai Boon Tan
 */
public class SimOverlapNitrogenModel extends Simulation {

    public SimOverlapNitrogenModel(Space _space, int numMolecules,double temperature, String filename, double scale) {
        super(_space);
        this.fname = filename;

        integrators = new IntegratorBox[2];
        accumulatorPumps = new DataPump[2];
        meters = new IDataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];

        double unitCellLength = scale * 5.661;
        nCell = (int) Math.round(Math.pow((numMolecules / 4), 1.0 / 3.0));

        Basis basisFCC = new BasisCubicFcc();
        Basis basis = new BasisBigCell(space, basisFCC, new int[]{nCell, nCell, nCell});

        ConformationNitrogen conformation = new ConformationNitrogen(space);
        species = new SpeciesN2(space);
        species.setConformation(conformation);
        addSpecies(species);

        potentialMasterTarget = new PotentialMaster();

        // TARGET
        boundaryTarget = new BoundaryDeformablePeriodic(space, nCell * unitCellLength);
        boxTarget = this.makeBox(boundaryTarget);
        boxTarget.setNMolecules(species, numMolecules);

        int[] nCells = new int[]{1, 1, 1};
        primitive = new PrimitiveCubic(space, nCell * unitCellLength);

        CoordinateDefinitionNitrogen coordDefTarget = new CoordinateDefinitionNitrogen(this, boxTarget, primitive, basis, space);
        coordDefTarget.setIsAlpha();
        coordDefTarget.setOrientationVectorAlpha(space);
        coordDefTarget.initializeCoordinates(nCells);

        double rCScale = 0.45;
        double rC = boxTarget.getBoundary().getBoxSize().getX(0) * rCScale;
        System.out.println("Truncation Radius (" + rCScale + " Box Length): " + rC);
        P2Nitrogen potential = new P2Nitrogen(space, rCScale);
        potential.setBox(boxTarget);
        System.out.println("Box Dimension(before): " + boxTarget.getBoundary().getBoxSize().toString());
        final Vector initBox = Vector.of(new double[]{boxTarget.getBoundary().getBoxSize().getX(0),
                boxTarget.getBoundary().getBoxSize().getX(1),
                boxTarget.getBoundary().getBoxSize().getX(2)});
        coordDefTarget.setInitVolume(initBox);

        P0LatticeEnergyCorrec p0correc = new P0LatticeEnergyCorrec(space);
        p0correc.setSpecies(species);
        p0correc.setBox(boxTarget);

        potentialMasterTarget.addPotential(potential, new ISpecies[]{species, species});
        potentialMasterTarget.addPotential(p0correc, new ISpecies[]{species, species});

        integratorTarget = new IntegratorMC(potentialMasterTarget, getRandom(), Kelvin.UNIT.toSim(temperature), boxTarget);
        MCMoveMoleculeCoupled move = new MCMoveMoleculeCoupled(potentialMasterTarget, getRandom(), space);
        move.setBox(boxTarget);
        move.setPotential(potential);

        MCMoveRotateMolecule3D rotate = new MCMoveRotateMolecule3D(potentialMasterTarget, getRandom(), space);
        rotate.setBox(boxTarget);

        MCMoveVolume mcMoveVolume = new MCMoveVolume(this, potentialMasterTarget, space);
        mcMoveVolume.setBox(boxTarget);
        mcMoveVolume.setPressure(Pascal.UNIT.toSim(0.0e9));

        integratorTarget = new IntegratorMC(this, potentialMasterTarget, boxTarget);
        integratorTarget.getMoveManager().addMCMove(move);
        integratorTarget.getMoveManager().addMCMove(rotate);
        integratorTarget.getMoveManager().addMCMove(mcMoveVolume);
        integrators[1] = integratorTarget;

        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMasterTarget, boxTarget);
        latticeEnergy = meterPE.getDataAsScalar();

        System.out.println("lattice energy per molecule in K: " + Kelvin.UNIT.fromSim(latticeEnergy) / numMolecules + "\n");

        // HARMONIC
        boundaryHarmonic = new BoundaryDeformablePeriodic(space, nCell * unitCellLength);
        boxHarmonic = this.makeBox(boundaryHarmonic);
        boxHarmonic.setNMolecules(species, numMolecules);

        integratorHarmonic = new IntegratorMC(null, random, 1.0, boxHarmonic); //null changed on 11/20/2009

        moveHarmonic = new MCMoveHarmonic(getRandom());
        integratorHarmonic.getMoveManager().addMCMove(moveHarmonic);
        integrators[0] = integratorHarmonic;

        CoordinateDefinitionNitrogen coordDefHarmonic = new CoordinateDefinitionNitrogen(this, boxHarmonic, primitive, basis, space);
        coordDefHarmonic.setIsAlpha();
        coordDefHarmonic.setOrientationVectorAlpha(space);
        coordDefHarmonic.initializeCoordinates(nCells);
        coordDefHarmonic.setInitVolume(initBox);

        normalModes = new NormalModesFromFile(filename, space.D());
        normalModes.setTemperature(Kelvin.UNIT.toSim(temperature));

        WaveVectorFactory waveVectorFactory = normalModes.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(boxHarmonic);
        moveHarmonic.setOmegaSquared(normalModes.getOmegaSquared());
        moveHarmonic.setEigenVectors(normalModes.getEigenvectors());
        moveHarmonic.setWaveVectors(waveVectorFactory.getWaveVectors());
        moveHarmonic.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
        moveHarmonic.setCoordinateDefinition(coordDefHarmonic);
        moveHarmonic.setTemperature(Kelvin.UNIT.toSim(temperature));
        //moveHarmonic.setModeNum(new int[] {4});

        moveHarmonic.setBox(boxHarmonic);

        // OVERLAP
        integratorOverlap = new IntegratorOverlap(new IntegratorBox[]{integratorHarmonic, integratorTarget});
        meterHarmonicEnergy = new MeterHarmonicEnergy(coordDefTarget, normalModes);
        MeterBoltzmannTarget meterTarget = new MeterBoltzmannTarget(new MeterPotentialEnergyFromIntegrator(integratorTarget), meterHarmonicEnergy);
        meterTarget.setTemperature(temperature);
        meterTarget.setLatticeEnergy(latticeEnergy);
        meters[1] = meterTarget;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, false), 1);

        MeterBoltzmannHarmonic meterHarmonic = new MeterBoltzmannHarmonic(moveHarmonic, potentialMasterTarget);
        meterHarmonic.setTemperature(Kelvin.UNIT.toSim(temperature));
        meterHarmonic.setLatticeEnergy(latticeEnergy);
        meters[0] = meterHarmonic;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, true), 0);

        setRefPref(1.0, 30);

        this.getController2().addActivity(new ActivityIntegrate2(integratorOverlap));
    }

    public void setRefPref(double refPrefCenter, double span) {
        refPref = refPrefCenter;
        accumulators[0].setBennetParam(refPrefCenter,span);
        accumulators[1].setBennetParam(refPrefCenter,span);

    }

    public void setAccumulator(AccumulatorVirialOverlapSingleAverage newAccumulator, int iBox) {

        accumulators[iBox] = newAccumulator;
    
        newAccumulator.setBlockSize(200); // setting the block size = 300
        
        if (accumulatorPumps[iBox] == null) {
            accumulatorPumps[iBox] = new DataPump(meters[iBox],newAccumulator);
            IntegratorListenerAction pumpListener = new IntegratorListenerAction(accumulatorPumps[iBox]);
            integrators[iBox].getEventManager().addListener(pumpListener);
            if (iBox == 1) {
            	if (boxTarget.getMoleculeList().size()==32){
            		
            	    pumpListener.setInterval(100);
            	
            	} else if (boxTarget.getMoleculeList().size()==108){
                
            	    pumpListener.setInterval(300);
            	} else 
            	    pumpListener.setInterval(500);
                    //pumpListener.setInterval(boxTarget.getMoleculeList().getMoleculeCount());
            }
        }
        else {
            accumulatorPumps[iBox].setDataSink(newAccumulator);
        }
        if (integratorOverlap != null && accumulators[0] != null && accumulators[1] != null) {
            dsvo = new DataSourceVirialOverlap(accumulators[0],accumulators[1]);
            integratorOverlap.setReferenceFracSource(dsvo);
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
            getController2().runActivityBlocking(new ActivityIntegrate2(integratorOverlap), initSteps/2);

            System.out.println("target equilibration finished");

            setAccumulator(new AccumulatorVirialOverlapSingleAverage(41,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(41,false),1);
            setRefPref(1,200);
getController2().runActivityBlocking(new ActivityIntegrate2(integratorOverlap), initSteps);


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

        }

    }
    
    public void equilibrate(String fileName, long initSteps) {
        // run a short simulation to get reasonable MC Move step sizes and
        // (if needed) narrow in on a reference preference
        for (int i=0; i<2; i++) {
            if (integrators[i] instanceof IntegratorMC) ((IntegratorMC)integrators[i]).getMoveManager().setEquilibrating(true);
        }
this.getController2().runActivityBlocking(new ActivityIntegrate2(this.integratorOverlap), initSteps);

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
     * @see SimOverlapNitrogenModel.SimOverlapParam
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
        int D = 3;
        double pressure = 0.0;
        long numSteps = params.numSteps;
        final int numMolecules = params.numMolecules;
        double temperature = params.temperature;
        double scale = params.scale;
        String filename = params.filename;
        if (filename.length() == 0) {
        	System.err.println("Need input files!!!");
        	if(temperature <1.0){
    			filename = "alphaN2_nA"+numMolecules+"_T0"+(int)(temperature*10);
    			
    		} else {
    			filename = "alphaN2_nA"+numMolecules+"_T"+Math.round(temperature);
    		}
        }
        //String refFileName = args.length > 0 ? filename+"_ref" : null;
        String refFileName = filename+"_ref";
        
    	System.out.println("Running alpha-N2 crystal structure overlap-sampling simulation with " + numSteps + " steps" );
		System.out.println("num Molecules: " + numMolecules+ " ; temperature: " + temperature
					+"K ; pressure: "+ pressure+"GPa");
		System.out.println("With volume scaling of " + scale);
		System.out.println((numSteps/1000)+" total steps of 1000");
        System.out.println("output data to "+filename +"\n");

        //instantiate simulation
        SimOverlapNitrogenModel sim = new SimOverlapNitrogenModel(Space.getInstance(D), numMolecules, temperature, filename, scale);
       
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
       
        sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integratorOverlap), numSteps);
        
        int totalCells = 1;
        for (int i=0; i<D; i++) {
            totalCells *= sim.nCell;
        }
        double  AHarmonic = Kelvin.UNIT.fromSim(CalcHarmonicA.doit(sim.normalModes, D, Kelvin.UNIT.toSim(temperature), numMolecules));
        System.out.println("Harmonic-reference free energy in K, A: "+AHarmonic + " " + AHarmonic/numMolecules);
        System.out.println(" ");
        
        System.out.println("final reference optimal step frequency "+sim.integratorOverlap.getIdealRefStepFraction()
        		+" (actual: "+sim.integratorOverlap.getRefStepFraction()+")");
              
        double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
        double ratio = ratioAndError[0];
        double error = ratioAndError[1];
        double uLatticeInfinite = -1023.896102;
		System.out.println("Ulattice energy(infinite) in K: "+ uLatticeInfinite);
		
        System.out.println("\nratio average: "+ratio+" ,error: "+error);
        System.out.println("free energy difference in K, deltaA: "+(-temperature*Kelvin.UNIT.fromSim(Math.log(ratio)))
        												  +" ,error: "+temperature*Kelvin.UNIT.fromSim((error/ratio)));
        System.out.println("target free energy in K, A: "+(AHarmonic-temperature*Kelvin.UNIT.fromSim(Math.log(ratio))));
        System.out.println("target free energy per particle in K, A/N: "+ (AHarmonic-temperature*Kelvin.UNIT.fromSim(Math.log(ratio)))/numMolecules 
        		+" ;error: "+temperature*Kelvin.UNIT.fromSim((error/ratio))/numMolecules);
        System.out.println("target Helmholtz free energy per particle in K, A/N: "+ ((AHarmonic-temperature*Kelvin.UNIT.fromSim(Math.log(ratio))+
        		Kelvin.UNIT.fromSim(sim.latticeEnergy))/numMolecules) 
        		+" ;error: "+temperature*Kelvin.UNIT.fromSim((error/ratio))/numMolecules);
        
        DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
        double betaFAW = -Math.log(((DataDoubleArray)allYourBase.getData(sim.accumulators[0].AVERAGE.index)).getData()[1]);
        System.out.println("harmonic ratio average: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].AVERAGE.index)).getData()[1]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].STANDARD_DEVIATION.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].ERROR.index)).getData()[1]);
        
        allYourBase = (DataGroup)sim.accumulators[1].getData(sim.dsvo.minDiffLocation());
        double betaFBW = -Math.log(((DataDoubleArray)allYourBase.getData(sim.accumulators[0].AVERAGE.index)).getData()[1]);
        System.out.println("target ratio average: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].AVERAGE.index)).getData()[1]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].STANDARD_DEVIATION.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulators[0].ERROR.index)).getData()[1]);
        
        
        long endTime = System.currentTimeMillis();
        System.out.println("End Time: " + endTime);
        System.out.println("Time taken: " + (endTime - startTime));
    
        
		if(false){
			SimulationGraphic simGraphic = new SimulationGraphic(sim);
		    simGraphic.getDisplayBox(sim.boxHarmonic).setPixelUnit(new Pixel(50));
		    simGraphic.makeAndDisplayFrame("Overlap Sampling Alpha-Phase Nitrogen Crystal Structure");
			MeterWorkHarmonicPhaseSpace meterHarmonicTarget = new MeterWorkHarmonicPhaseSpace(sim.moveHarmonic, sim.potentialMasterTarget);
			meterHarmonicTarget.setLatticeEnergy(sim.latticeEnergy);
			meterHarmonicTarget.setTemperature(Kelvin.UNIT.toSim(temperature));

			AccumulatorAverage energyAverage = new AccumulatorAverageCollapsing();
			DataPump energyPump = new DataPump(meterHarmonicTarget, energyAverage);

			IntegratorListenerAction energyListener = new IntegratorListenerAction(energyPump);
			energyListener.setInterval(1);
			sim.integratorHarmonic.getEventManager().addListener(energyListener);
sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integratorOverlap), numSteps);
		}
        
    }

    private static final long serialVersionUID = 1L;
    public IntegratorOverlap integratorOverlap;
    public DataSourceVirialOverlap dsvo;
    public IntegratorBox[] integrators;

    public Box boxTarget, boxHarmonic;
    public Boundary boundaryTarget, boundaryHarmonic;
    public NormalModes normalModes;
    public Primitive primitive, primitiveUnitCell;
    public double refPref;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public DataPump[] accumulatorPumps;
    public IDataSource[] meters;
    public String fname;
    protected MCMoveHarmonic moveHarmonic;
    protected PotentialMaster potentialMasterTarget;
    protected MeterHarmonicEnergy meterHarmonicEnergy;
    protected double latticeEnergy;
    protected SpeciesN2 species;
    protected int nCell;
    protected IntegratorMC integratorHarmonic, integratorTarget;
    
    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules =32;
        public long numSteps = 1000000;
        public String filename = "NPT32T35";
        public double temperature =35.0;
        public double scale = 1.0;
    }
}
