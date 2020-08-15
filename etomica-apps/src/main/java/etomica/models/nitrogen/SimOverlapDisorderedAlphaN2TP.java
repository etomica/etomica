/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.WriteConfiguration;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.nbr.list.molecule.PotentialMasterListMolecular;
import etomica.normalmode.BasisBigCell;
import etomica.normalmode.MCMoveMoleculeCoupled;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.species.Species;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;

import java.awt.*;
import java.io.File;

/**
 * Temperature-perturbation simulation for disordered alpha-phase Nitrogen
 * 
 * @author Tai Boon Tan
 */
public class SimOverlapDisorderedAlphaN2TP extends Simulation {

    public SimOverlapDisorderedAlphaN2TP(Space space, int numMolecules, double density, double temperature, double[] otherTemperatures,
    		double[] alpha, int numAlpha, double alphaSpan, long numSteps, boolean isBeta, boolean isBetaHCP, double rcScale) {
        super(space);

        int nC = (int) Math.pow(numMolecules / 3.999999999, 1.0 / 3.0);
        double a = Math.pow(4.0 / density, 1.0 / 3.0);
        System.out.println("Unit Cell Length, a: " + a);

        Basis basisFCC = new BasisCubicFcc();
        Basis basis = new BasisBigCell(space, basisFCC, new int[]{nC, nC, nC});

        Species species = new SpeciesN2(space);
        addSpecies(species);

        Boundary boundary = new BoundaryRectangularPeriodic(space, nC * a);
        box = this.makeBox(boundary);
        box.setNMolecules(species, numMolecules);

        int[] nCells = new int[]{1, 1, 1};
        Primitive primitive = new PrimitiveCubic(space, nC * a);

        coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
        coordinateDef.setIsAlpha();
        coordinateDef.setOrientationVectorAlpha(space);
        coordinateDef.initializeCoordinates(nCells);

        double rc = box.getBoundary().getBoxSize().getX(0) * rcScale;
        System.out.println("Truncation Radius (" + rcScale + " Box Length): " + rc);
        P2Nitrogen potential = new P2Nitrogen(space, rc);
        potential.setBox(box);

        //potentialMaster = new PotentialMaster();
        potentialMaster = new PotentialMasterListMolecular(this, space);
        potentialMaster.addPotential(potential, new ISpecies[]{species, species});

        int cellRange = 6;
        potentialMaster.setRange(rc);
        potentialMaster.setCellRange(cellRange);
        potentialMaster.getNeighborManager(box).reset();
        potential.setRange(Double.POSITIVE_INFINITY);

//        int potentialCells = potentialMaster.getNbrCellManager(box).getLattice().getSize()[0];
//        if (potentialCells < cellRange*2+1) {
//            throw new RuntimeException("oops ("+potentialCells+" < "+(cellRange*2+1)+")");
//        }
        int numNeigh = potentialMaster.getNeighborManager(box).getUpList(box.getMoleculeList().get(0))[0].size();
        System.out.println("numNeigh: " + numNeigh);

        move = new MCMoveMoleculeCoupled(potentialMaster, getRandom(), space);
        move.setBox(box);
        move.setPotential(potential);
        move.setDoExcludeNonNeighbors(true);
        //move.setStepSize(Kelvin.UNIT.toSim(temperature));
        //((MCMoveStepTracker)move.getTracker()).setNoisyAdjustment(true);

        integrator = new IntegratorMC(potentialMaster, getRandom(), Kelvin.UNIT.toSim(temperature), box);
        integrator.getMoveManager().addMCMove(move);
        if (isBetaHCP) {
            rotate = new MCMoveRotateMolecule3D(potentialMaster, getRandom(), space);
            rotate.setBox(box);
            integrator.getMoveManager().addMCMove(rotate);
        }

        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster, box);
        latticeEnergy = meterPE.getDataAsScalar();
        System.out.println("lattice energy per molecule (sim unit): " + latticeEnergy / numMolecules);
        System.out.println("lattice energy (sim unit): " + latticeEnergy);

        potential.setRange(rc);
        meter = new MeterTargetTPMolecule(potentialMaster, species, this, coordinateDef);
        meter.setLatticeEnergy(latticeEnergy);
        meter.setBetaPhase(true);
        meter.setTemperature(Kelvin.UNIT.toSim(temperature));
        meter.setOtherTemperatures(otherTemperatures);
        meter.setAlpha(alpha);
        meter.setAlphaSpan(alphaSpan);
        meter.setNumAlpha(numAlpha);
        potential.setRange(Double.POSITIVE_INFINITY);

        int numBlocks = 100;
        int interval = numMolecules;
        long blockSize = numSteps / (numBlocks * interval);
        if (blockSize == 0) blockSize = 1;
        System.out.println("block size " + blockSize + " interval " + interval);
        if (otherTemperatures.length > 1) {
            accumulator = new AccumulatorAverageCovariance(blockSize);
        } else {
            accumulator = new AccumulatorAverageFixed(blockSize);
        }
        accumulatorPump = new DataPumpListener(meter, accumulator, numMolecules);
        integrator.getEventManager().addListener(accumulatorPump);

        this.getController().addActivity(new ActivityIntegrate(integrator));

    }
    
    public void initialize(long initSteps) {
        // equilibrate off the lattice to avoid anomolous contributions
        System.out.println("\nEquilibration Steps: " + initSteps);
    	this.getController().runActivityBlocking(new ActivityIntegrate(this.integrator), initSteps);

        
        accumulator.reset();
    }
    
    public void initializeConfigFromFile(String fname){
        ConfigurationFile config = new ConfigurationFile(fname);
        config.initializeCoordinates(box);
    }
    
    public void writeConfiguration(String fname){
        WriteConfiguration writeConfig = new WriteConfiguration(space);
        writeConfig.setBox(box);
        writeConfig.setDoApplyPBC(false);
        writeConfig.setConfName(fname);
        writeConfig.actionPerformed();
        System.out.println("\n***output configFile: "+ fname);
    }
    
    /**
     * @param args filename containing simulation parameters
     * @see SimOverlapDisorderedAlphaN2TP.SimOverlapParam
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
        final int numMolecules = params.numMolecules;
        double temperature = params.temperature;
        double[] otherTemperatures = params.otherTemperatures;
        double[] alpha = params.alpha;
        int numAlpha = params.numAlpha;
        double alphaSpan = params.alphaSpan;
        boolean isBeta = params.isBeta;
        boolean isBetaHCP = params.isBetaHCP;
        double rcScale = params.rcScale;
        String configFileName = "configT"+temperature;
        
        System.out.println("Running disordered alpha-phase Nitrogen TP overlap simulation");
        System.out.println(numMolecules+" atoms at density "+density+" and temperature "+temperature + " K");
        System.out.print("perturbing into: ");
        for(int i=0; i<otherTemperatures.length; i++){
        	System.out.print(otherTemperatures[i]+" ");
        }
        System.out.print("\nwith alpha: ");
        for(int i=0; i<alpha.length; i++){
        	System.out.print(alpha[i]+" ");
        }
        System.out.println("\n"+numSteps+" steps");
        System.out.println("isBeta: " + isBeta);
        System.out.println("isBetaHCP: " + isBetaHCP);

        //instantiate simulation
        final SimOverlapDisorderedAlphaN2TP sim = new SimOverlapDisorderedAlphaN2TP(Space.getInstance(3), numMolecules, density, temperature, otherTemperatures, 
        		alpha, numAlpha, alphaSpan, numSteps, isBeta, isBetaHCP, rcScale);
        //start simulation

    	File configFile = new File(configFileName+".pos");
		if(configFile.exists()){
			System.out.println("\n***initialize coordinate from "+ configFile);
        	sim.initializeConfigFromFile(configFileName);
		} else {
			long initStep = (1+(numMolecules/500))*100*numMolecules;
			sim.initialize(initStep);
		}
        System.out.flush();
        
        if (false) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
//            simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(10));
//		  
//		    DiameterHashByType diameter = new DiameterHashByType(sim);
//			diameter.setDiameter(sim.species.getNitrogenType(), 3.1);
//			diameter.setDiameter(sim.species.getPType(), 0.0);
//			
//			simGraphic.getDisplayBox(sim.box).setDiameterHash(diameter);
            
            simGraphic.setPaintInterval(sim.box, 1000);           
            ColorScheme colorScheme = new ColorScheme() {
                public Color getAtomColor(IAtom a) {
                    if (allColors==null) {
                        allColors = new Color[768];
                        for (int i=0; i<256; i++) {
                            allColors[i] = new Color(255-i,i,0);
                        }
                        for (int i=0; i<256; i++) {
                            allColors[i+256] = new Color(0,255-i,i);
                        }
                        for (int i=0; i<256; i++) {
                            allColors[i+512] = new Color(i,0,255-i);
                        }
                    }
                    return allColors[(2*a.getLeafIndex()) % 768];
                }
                protected Color[] allColors;
            };
            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);
            
            DisplayTextBox timer = new DisplayTextBox();
            DataSourceCountSteps counter = new DataSourceCountSteps(sim.integrator);
            DataPumpListener counterPump = new DataPumpListener(counter, timer, 100);
            sim.integrator.getEventManager().addListener(counterPump);
            simGraphic.getPanel().controlPanel.add(timer.graphic());
            
            simGraphic.makeAndDisplayFrame("TP Alpha Nitrogen");
            return;
        }


        
        final long startTime = System.currentTimeMillis();
       
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator), numSteps);
        
        sim.writeConfiguration(configFileName);
        System.out.println("step size: " + sim.move.getStepSize());
        System.out.println("\nratio averages:\n");

        DataGroup data = (DataGroup)sim.accumulator.getData();
        IData dataErr = data.getData(sim.accumulator.ERROR.index);
        IData dataAvg = data.getData(sim.accumulator.AVERAGE.index);
        IData dataCorrelation = data.getData(sim.accumulator.BLOCK_CORRELATION.index);
        for (int i=0; i<otherTemperatures.length; i++) {
    
        	if(otherTemperatures[i] < 10.0){
        		System.out.printf("0"+ "%2.3f\n", otherTemperatures[i]);
                       
            } else {
            	System.out.printf("%2.3f\n", otherTemperatures[i]);
            }
            double[] iAlpha = sim.meter.getAlpha(i);
            for (int j=0; j<numAlpha; j++) {
                System.out.println("  "+iAlpha[j]+" "+dataAvg.getValue(i*numAlpha+j)
                        +" "+dataErr.getValue(i*numAlpha+j)
                        +" "+dataCorrelation.getValue(i*numAlpha+j));
            }
        }
        
        if (otherTemperatures.length == 2) {
            // we really kinda want each covariance for every possible pair of alphas,
            // but we're going to be interpolating anyway and the covariance is almost
            // completely insensitive to choice of alpha.  so just take the covariance for
            // the middle alphas.
            IData dataCov = data.getData(((AccumulatorAverageCovariance)sim.accumulator).BLOCK_COVARIANCE.index);
            System.out.print("covariance "+otherTemperatures[1]+" / "+otherTemperatures[0]+"   ");
            for (int i=0; i<numAlpha; i++) {
                i = (numAlpha-1)/2;
                double ivar = dataErr.getValue(i);
                ivar *= ivar;
                ivar *= (sim.accumulator.getBlockCount()-1);
                for (int j=0; j<numAlpha; j++) {
                    j = (numAlpha-1)/2;
                    double jvar = dataErr.getValue(numAlpha+j);
                    jvar *= jvar;
                    jvar *= (sim.accumulator.getBlockCount()-1);
                    System.out.print(dataCov.getValue(i*(2*numAlpha)+(numAlpha+j))/Math.sqrt(ivar*jvar)+" ");
                    break;
                }
                System.out.println();
                break;
            }
        }
        
      
        long endTime = System.currentTimeMillis();
        System.out.println("Time taken(s): " + (endTime - startTime)/1000.0);
    }

    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;

    public Box box;
    public Primitive primitive;
    public AccumulatorAverageFixed accumulator;
    public DataPumpListener accumulatorPump;
    public MeterTargetTPMolecule meter;
    protected PotentialMasterListMolecular potentialMaster;
    protected double latticeEnergy;
    protected SpeciesN2 species;
    protected CoordinateDefinitionNitrogen coordinateDef;
    protected P2Nitrogen potential;
    protected MCMoveRotateMolecule3D rotate;
    protected MCMoveMoleculeCoupled move;
    
    /**
     * Inner class for parameters understood by the SimOverlapBetaN2TP constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules = 864;
        public double density = 0.0220; //0.02204857502170207 (intial from literature with a = 5.661)
        public long numSteps = 100000;
        public double temperature = 0.02; // in unit Kelvin
        public double[] alpha = new double[]{1.0};
        public int numAlpha = 11;
        public double alphaSpan = 1;
        public double[] otherTemperatures = new double[]{0.03};
        public boolean isBeta = true;
        public boolean isBetaHCP = false;
        public double rcScale = 0.475;
    }
}
