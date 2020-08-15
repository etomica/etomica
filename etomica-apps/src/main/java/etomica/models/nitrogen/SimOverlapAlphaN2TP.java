/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.WriteConfiguration;

import etomica.action.activity.ActivityIntegrate2;
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
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveTetragonal;
import etomica.integrator.IntegratorListenerAction;
import etomica.nbr.list.molecule.PotentialMasterListMolecular;
import etomica.normalmode.BasisBigCell;
import etomica.normalmode.MCMoveMoleculeCoupled;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;

import java.awt.*;
import java.io.File;

/**
 * Temperature-perturbation simulation for alpha-phase Nitrogen
 * 
 * @author Tai Boon Tan
 */
public class SimOverlapAlphaN2TP extends Simulation {

    public SimOverlapAlphaN2TP(Space space, int[] nC, double density, double temperature, double[] otherTemperatures,
    		double[] alpha, int numAlpha, double alphaSpan, long numSteps, double rcScale, double constraintAngle, boolean noRotScale) {
        super(space);

        int numMolecules = nC[0] * nC[1] * nC[2] * 4;

        double a = Math.pow(4.0 / density, 1.0 / 3.0);
        System.out.println("Unit Cell Length, a: " + a);

        Basis basisFCC = new BasisCubicFcc();
        Basis basis = new BasisBigCell(space, basisFCC, new int[]{nC[0], nC[1], nC[2]});

        species = new SpeciesN2(space);
        addSpecies(species);

        double[] boxSize = new double[]{nC[0] * a, nC[1] * a, nC[2] * a};
        Boundary boundary = new BoundaryRectangularPeriodic(space, boxSize);
        box = this.makeBox(boundary);
        box.setNMolecules(species, numMolecules);

        int[] nCells = new int[]{1, 1, 1};
        primitive = new PrimitiveTetragonal(space, nC[0] * a, nC[2] * a);

        coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
        coordinateDef.setIsAlpha();
        coordinateDef.setOrientationVectorAlpha(space);
        coordinateDef.initializeCoordinates(nCells);
        double rC = box.getBoundary().getBoxSize().getX(0) * rcScale;
        System.out.println("Truncation Radius (" + rcScale + " Box Length): " + rC);
        potential = new P2Nitrogen(space, rC);
        potential.setBox(box);

        pRotConstraint = new PRotConstraint(space, coordinateDef, box);
        pRotConstraint.setConstraintAngle(constraintAngle);

        //potentialMaster = new PotentialMaster();
        potentialMaster = new PotentialMasterListMolecular(this, space);
        potentialMaster.addPotential(potential, new ISpecies[]{species, species});
        if (!noRotScale) {
            System.out.println("set constraint angle to = " + constraintAngle);
            potentialMaster.addPotential(pRotConstraint, new ISpecies[]{species});
        }

        int cellRange = 6;
        potentialMaster.setRange(rC);
        potentialMaster.setCellRange(cellRange);
        potentialMaster.getNeighborManager(box).reset();

        int potentialCells = potentialMaster.getNbrCellManager(box).getLattice().getSize()[0];
        if (potentialCells < cellRange * 2 + 1) {
            throw new RuntimeException("oops (" + potentialCells + " < " + (cellRange * 2 + 1) + ")");
        }

        int numNeigh = potentialMaster.getNeighborManager(box).getDownList(box.getMoleculeList().get(0))[0].size();
        System.out.println("numNeigh: " + numNeigh);

        MCMoveMoleculeCoupled move = new MCMoveMoleculeCoupled(potentialMaster, getRandom(), space);
        move.setBox(box);
        move.setPotential(potential);
        move.setDoExcludeNonNeighbors(true);

        rotate = new MCMoveRotateMolecule3DN2AveCosThetaConstraint(potentialMaster, getRandom(), space, coordinateDef, 0.8);
        rotate.setBox(box);

        integrator = new IntegratorMC(potentialMaster, getRandom(), Kelvin.UNIT.toSim(temperature), box);
        integrator.getMoveManager().addMCMove(move);
        integrator.getMoveManager().addMCMove(rotate);

        potential.setRange(Double.POSITIVE_INFINITY);

        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster, box);
        latticeEnergy = meterPE.getDataAsScalar();
        System.out.println("lattice energy per molecule (sim unit): " + latticeEnergy / numMolecules);

        potential.setRange(rC);
        meter = new MeterTargetTPMolecule(potentialMaster, species, this, coordinateDef);
        meter.setLatticeEnergy(latticeEnergy);
        meter.setTemperature(Kelvin.UNIT.toSim(temperature));
        meter.setOtherTemperatures(otherTemperatures);
        meter.setAlpha(alpha);
        meter.setAlphaSpan(alphaSpan);
        meter.setNumAlpha(numAlpha);
        potential.setRange(Double.POSITIVE_INFINITY);

        if (noRotScale) {
            System.out.println("**** NOT SCALING THE ROTATION!");
            meter.setBetaPhase(true);
        }

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
        accumulatorPump = new DataPumpListener(meter, accumulator, interval);
        integrator.getEventManager().addListener(accumulatorPump);

        this.getController().addActivity(new ActivityIntegrate2(integrator));

    }
    
    public void initialize(long initSteps) {
        // equilibrate off the lattice to avoid anomolous contributions
        System.out.println("\nEquilibration Steps: " + initSteps);
    	this.getController().runActivityBlocking(new ActivityIntegrate2(this.integrator), initSteps);

        
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
     * @see SimOverlapAlphaN2TP.SimOverlapParam
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
        int[] nC = params.nC;
        int numAlpha = params.numAlpha;
        double alphaSpan = params.alphaSpan;
        double rcScale = params.rcScale;
        double constraintAngle = params.constraintAngle;
        boolean noRotScale = params.noRotScale;
        String configFileName = "configT"+temperature;
        String filename = "alphaN2d"+density+"_T"+temperature+"Cons0.8";
		
        System.out.println("Running alpha-phase Nitrogen TP overlap simulation");
        System.out.println(numMolecules+" molecules at density "+density+" and temperature "+temperature + " K");
        System.out.print("perturbing into: ");
        for(int i=0; i<otherTemperatures.length; i++){
        	System.out.print(otherTemperatures[i]+" ");
        }
        System.out.print("\nwith alpha: ");
        for(int i=0; i<alpha.length; i++){
        	System.out.print(alpha[i]+" ");
        }
        System.out.println("\n"+numSteps+" steps");

        //instantiate simulation
        final SimOverlapAlphaN2TP sim = new SimOverlapAlphaN2TP(Space.getInstance(3), nC, density, temperature, otherTemperatures, 
        		alpha, numAlpha, alphaSpan, numSteps, rcScale, constraintAngle, noRotScale);
        
        //start simulation
        File configFile = new File(configFileName+".pos");
		if(configFile.exists()){
			System.out.println("\n***initialize coordinate from "+ configFile);
        	sim.initializeConfigFromFile(configFileName);
            sim.rotate.calcAveCosThetaInitial();    
		} else {
	        long initStep = (1+(numMolecules/500))*100*numMolecules;
	        sim.initialize(initStep);
		}
		System.out.flush();
	
		
        if (false) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
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
       
    	MeterRotationDistributionGroup meterRotation = new MeterRotationDistributionGroup(sim.box, sim.coordinateDef);
		IntegratorListenerAction meterRotationListener = new IntegratorListenerAction(meterRotation);
		meterRotationListener.setInterval(numMolecules);                                      
		sim.integrator.getEventManager().addListener(meterRotationListener);
        
		MeterOrientationOrderParameter meterOrientationOrderParameter = new MeterOrientationOrderParameter(sim.coordinateDef);
		
		AccumulatorAverage orderParameterAverage = new AccumulatorAverageFixed();
		DataPump orderParameterPump = new DataPump(meterOrientationOrderParameter, orderParameterAverage);
		IntegratorListenerAction orderParameterListener = new IntegratorListenerAction(orderParameterPump);
		orderParameterListener.setInterval(numMolecules);
		sim.integrator.getEventManager().addListener(orderParameterListener);
		
		
        sim.getController().runActivityBlocking(new ActivityIntegrate2(sim.integrator), numSteps);
        //MeterTargetTP.openFW("x"+numMolecules+".dat");
        //MeterTargetTP.closeFW();
        System.out.println("PRotConstraint counter: " + sim.pRotConstraint.counter);
        sim.writeConfiguration(configFileName);
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
        meterRotation.writeUdistribution(filename);
		double averageOrderParameter = orderParameterAverage.getData().getValue(sim.accumulator.AVERAGE.index);
//		System.out.println("orientational order parameter: "+ averageOrderParameter/(numMolecules*numMolecules));
		
		System.out.println("orientational order parameter: "+ 0.5*averageOrderParameter);
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
    protected PRotConstraint pRotConstraint;
    protected MCMoveRotateMolecule3DN2AveCosThetaConstraint rotate;
    
    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules = 864;
        public int[] nC = new int[]{6,6,6};
        public double density = 0.023; //0.02204857502170207 (intial from literature with a = 5.661)
        public long numSteps = 100000;
        public double temperature = 0.01; // in unit Kelvin
        public double[] alpha = new double[]{1.0};
        public int numAlpha = 11;
        public double alphaSpan = 1;
        public double[] otherTemperatures = new double[]{0.02};
        public double rcScale = 0.475;
        public double constraintAngle = 70;
        public boolean noRotScale = false;
    }
}
