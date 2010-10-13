package etomica.models.nitrogen;

import java.awt.Color;
import java.io.File;

import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.config.ConfigurationFile;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorRatioAverageCovariance;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.IData;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.nbr.list.molecule.BoxAgentSourceCellManagerListMolecular;
import etomica.normalmode.BasisBigCell;
import etomica.normalmode.MCMoveMoleculeCoupled;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.units.Degree;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;

/**
 * Temperature-perturbation simulation for beta-phase Nitrogen
 * 
 * @author Tai Boon Tan
 */
public class SimOverlapBetaN2TP extends Simulation {

    public SimOverlapBetaN2TP(Space space, int numMolecules, double density, double temperature, double[] otherTemperatures,
    		double[] alpha, int numAlpha, double alphaSpan, long numSteps, boolean isBeta, boolean isBetaHCP, boolean scaleRot) {
        super(space);
        
        BoxAgentSourceCellManagerListMolecular boxAgentSource = new BoxAgentSourceCellManagerListMolecular(this, null, space);
        BoxAgentManager boxAgentManager = new BoxAgentManager(boxAgentSource);
     
		double ratio = 1.631;
		double aDim = Math.pow(4.0/(Math.sqrt(3.0)*ratio*density), 1.0/3.0);
		double cDim = aDim*ratio;
		System.out.println("a: " + aDim + " ;cDim: " + cDim);
		int nC = (int)Math.pow(numMolecules/1.999999999, 1.0/3.0);
		
		Basis basisHCP = new BasisHcp();
		BasisBigCell basis = new BasisBigCell(space, basisHCP, new int[]{nC,nC,nC});
        
		species = new SpeciesN2(space);
		addSpecies(species);

		box = new Box(space);
		addBox(box);
		box.setNMolecules(species, numMolecules);
        
		IVector[] boxDim = new IVector[3];
		boxDim[0] = space.makeVector(new double[]{nC*aDim, 0, 0});
		boxDim[1] = space.makeVector(new double[]{-nC*aDim*Math.cos(Degree.UNIT.toSim(60)), nC*aDim*Math.sin(Degree.UNIT.toSim(60)), 0});
		boxDim[2] = space.makeVector(new double[]{0, 0, nC*cDim});
		
		int[] nCells = new int[]{1,1,1};
		Boundary boundary = new BoundaryDeformablePeriodic(space, boxDim);
		primitive = new PrimitiveHexagonal(space, nC*aDim, nC*cDim);
		
		coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		if(isBeta){coordinateDef.setIsBeta();}
		if(isBetaHCP){coordinateDef.setIsBetaHCP();}
        coordinateDef.setOrientationVectorBeta(space);
		coordinateDef.initializeCoordinates(nCells);
			
        box.setBoundary(boundary);
		double rCScale = 0.475;
		double rc = aDim*nC*rCScale;
		System.out.println("Truncation Radius (" + rCScale +" Box Length): " + rc);
		potential = new P2Nitrogen(space, rc);
		potential.setBox(box);

		pRotConstraint = new PRotConstraint(space,coordinateDef,box);
		pRotConstraint.setConstraintAngle(1);
		
		potentialMaster = new PotentialMaster();
		//potentialMaster = new PotentialMasterListMolecular(this, rc, boxAgentSource, boxAgentManager, new NeighborListManagerSlantyMolecular.NeighborListSlantyAgentSourceMolecular(rc, space), space);
	    potentialMaster.addPotential(potential, new ISpecies[]{species, species});
		if(scaleRot){potentialMaster.addPotential(pRotConstraint,new ISpecies[]{species} );}
		
//	    int cellRange = 6;
//        potentialMaster.setRange(rc);
//        potentialMaster.setCellRange(cellRange); 
//        potentialMaster.getNeighborManager(box).reset();
//        
//        int potentialCells = potentialMaster.getNbrCellManager(box).getLattice().getSize()[0];
//        if (potentialCells < cellRange*2+1) {
//            throw new RuntimeException("oops ("+potentialCells+" < "+(cellRange*2+1)+")");
//        }
//	
//        int numNeigh = potentialMaster.getNeighborManager(box).getUpList(box.getMoleculeList().getMolecule(0))[0].getMoleculeCount();
//        System.out.println("numNeigh: " + numNeigh);
		
		MCMoveMoleculeCoupled move = new MCMoveMoleculeCoupled(potentialMaster,getRandom(),space);
		move.setBox(box);
		move.setPotential(potential);
		
		MCMoveRotateMolecule3D rotate = new MCMoveRotateMolecule3D(potentialMaster, getRandom(), space);
		rotate.setBox(box);
			
		integrator = new IntegratorMC(potentialMaster, getRandom(), Kelvin.UNIT.toSim(temperature));
		integrator.getMoveManager().addMCMove(move);
		if(isBetaHCP|| scaleRot){integrator.getMoveManager().addMCMove(rotate);}
		integrator.setBox(box);
		
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(box);
        latticeEnergy = meterPE.getDataAsScalar();
        System.out.println("lattice energy per molecule (K): " + Kelvin.UNIT.fromSim(latticeEnergy)/numMolecules);
        System.out.println("lattice energy per molecule: " + latticeEnergy);
        meter = new MeterTargetTPMolecule(potentialMaster, species, space, this, coordinateDef);
        meter.setLatticeEnergy(latticeEnergy);
        if(!scaleRot){meter.setBetaPhase(true);}
        meter.setTemperature(Kelvin.UNIT.toSim(temperature));
        meter.setOtherTemperatures(otherTemperatures);
        meter.setAlpha(alpha);
        meter.setAlphaSpan(alphaSpan);
        meter.setNumAlpha(numAlpha);
       
        int numBlocks = 100;
        int interval = numMolecules;
        long blockSize = numSteps/(numBlocks*interval);
        if (blockSize == 0) blockSize = 1;
        System.out.println("block size "+blockSize+" interval "+interval);
        if (otherTemperatures.length > 1) {
            accumulator = new AccumulatorAverageCovariance(blockSize);
        }
        else {
            accumulator = new AccumulatorAverageFixed(blockSize);
        }
        accumulatorPump = new DataPumpListener(meter, accumulator, interval);
        integrator.getEventManager().addListener(accumulatorPump);
       
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
     
    }
    
    public void initialize(long initSteps) {
        // equilibrate off the lattice to avoid anomolous contributions
        System.out.println("\nEquilibration Steps: " + initSteps);
    	activityIntegrate.setMaxSteps(initSteps);
        getController().actionPerformed();
        getController().reset();
        
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
     * @see SimOverlapBetaN2TP.SimOverlapParam
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
        boolean scaleRot = params.scaleRot;
        String configFileName = "configT"+temperature;
        
        System.out.println("Running beta-phase Nitrogen TP overlap simulation");
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
        System.out.println("scaleRot: " + scaleRot);

        //instantiate simulation
        final SimOverlapBetaN2TP sim = new SimOverlapBetaN2TP(Space.getInstance(3), numMolecules, density, temperature, otherTemperatures, 
        		alpha, numAlpha, alphaSpan, numSteps, isBeta, isBetaHCP, scaleRot);
        //start simulation

    	File configFile = new File(configFileName+".pos");
		if(configFile.exists()){
			System.out.println("\n***initialize coordinate from "+ configFile);
        	sim.initializeConfigFromFile(configFileName);
		} else {
			long initStep = (1+(numMolecules/1000))*100*numMolecules;
			sim.initialize(initStep);
		}
        System.out.flush();
        
        if (false) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, sim.space, sim.getController());
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
       
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().actionPerformed();
        
        System.out.println("PRotConstraint counter: " + sim.pRotConstraint.counter);
        sim.writeConfiguration(configFileName);
        System.out.println("\nratio averages:\n");

        DataGroup data = (DataGroup)sim.accumulator.getData();
        IData dataErr = data.getData(AccumulatorAverage.StatType.ERROR.index);
        IData dataAvg = data.getData(AccumulatorAverage.StatType.AVERAGE.index);
        IData dataCorrelation = data.getData(AccumulatorRatioAverageCovariance.StatType.BLOCK_CORRELATION.index);
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
            IData dataCov = data.getData(AccumulatorAverageCovariance.StatType.BLOCK_COVARIANCE.index);
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
    public ActivityIntegrate activityIntegrate;
    public IBox box;
    public Primitive primitive;
    public AccumulatorAverageFixed accumulator;
    public DataPumpListener accumulatorPump;
    public MeterTargetTPMolecule meter;
    protected PotentialMaster potentialMaster;
    protected double latticeEnergy;
    protected SpeciesN2 species;
    protected CoordinateDefinitionNitrogen coordinateDef;
    protected P2Nitrogen potential;
    protected PRotConstraint pRotConstraint;
    
    /**
     * Inner class for parameters understood by the SimOverlapBetaN2TP constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules = 432;
        public double density = 0.025; //0.02204857502170207 (intial from literature with a = 5.661)
        public long numSteps = 100000;
        public double temperature = 1.0; // in unit Kelvin
        public double[] alpha = new double[]{1.0};
        public int numAlpha = 11;
        public double alphaSpan = 1;
        public double[] otherTemperatures = new double[]{1.5};
        public boolean isBeta = true;
        public boolean isBetaHCP = false;
        public boolean scaleRot = false;
    }
}
