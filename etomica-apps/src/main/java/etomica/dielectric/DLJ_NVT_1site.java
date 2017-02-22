package etomica.dielectric;

import java.awt.Color;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IPotentialMaster;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.DiameterHashByType;
import etomica.atom.DipoleSource;
import etomica.atom.IAtomOriented;
import etomica.atom.IAtomPositionDefinition;
import etomica.atom.IAtomTypeOriented;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPump;
import etomica.data.IData;
import etomica.data.meter.MeterDipoleSumSquared1site;
import etomica.data.meter.MeterDipoleSumSquaredMappedAverage;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayBoxCanvasG3DSys.OrientedSite;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.lattice.LatticeCubicBcc;
import etomica.listener.IntegratorListenerAction;
import etomica.potential.P2LJDipole;
import etomica.potential.P2ReactionFieldDipole;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Canonical ensemble Monte Carlo simulation (NVT)
 * calculate dielectric constant epsilon for dipolar LJ model
 * @author Weisong & shu
 * July 2015
 */

public class DLJ_NVT_1site extends Simulation {

	private static final long serialVersionUID = 1L;
	protected final IPotentialMaster potentialMaster;
	protected final IntegratorMC integrator;
	protected final MCMoveMolecule moveMolecule;//translation
	protected final MCMoveRotate rotateMolecule;//rotation, atomic
	protected final IBox box;
	protected SpeciesSpheresRotating species;
	private final static String APP_NAME = "dipolar LJ";
	private static final int PIXEL_SIZE = 15;
	public final ActivityIntegrate activityIntegrate;
    public Controller controller; 

    public static class DipoleSourceDLJ implements DipoleSource{//for potential reaction field
    	protected final IVectorMutable dipoleVector;
    	protected double dipoleStrength;
    	public DipoleSourceDLJ(ISpace space,double s){
    		dipoleStrength=s;
    		dipoleVector=space.makeVector();
    	}
    	
//    	public void setDipoleStrength(double dipoleStrength){
//    		this.dipoleStrength=dipoleStrength;
//    	}
		public IVector getDipole(IMolecule molecule) {
			IAtomList atomList = molecule.getChildList();
			if(atomList.getAtomCount() !=1){
				throw new RuntimeException("improper input of the molecule");
			}
	        IAtomOriented atom = (IAtomOriented)atomList.getAtom(0);
	    	dipoleVector.E(atom.getOrientation().getDirection());
			dipoleVector.TE(dipoleStrength);
			return dipoleVector; 
		}
    }
    
	//************************************* constructor ********************************************//
	public DLJ_NVT_1site(Space space, int numberMolecules,final double sigmaLJ,double epsilonLJ, double mu, 
			double dielectricOutside, double boxSize, double temperature,double truncation){
		super(space);		
//		setRandom(new RandomNumberGenerator(3));//Debug only
		species = new SpeciesSpheresRotating(space, new ElementSimple("A")); 
        addSpecies(species);
		box = new Box(space);
		addBox(box);
		box.setNMolecules(species, numberMolecules);
		box.getBoundary().setBoxSize(space.makeVector(new double[]{boxSize,boxSize,boxSize}));
	
		IAtomPositionDefinition positionDefinition = new IAtomPositionDefinition() {
			public IVector position(IMolecule molecule) {
				return molecule.getChildList().getAtom(0).getPosition();
			}
		};
		
		// dipolar LJ potential
		P2LJDipole pDLJ = new P2LJDipole(space,sigmaLJ,epsilonLJ,mu,truncation);
		// add reaction field potential
		DipoleSourceDLJ dipoleSourceDLJ = new DipoleSourceDLJ(space,mu);// add reaction field potential
        P2ReactionFieldDipole pRF = new P2ReactionFieldDipole(space,positionDefinition);
        pRF.setDipoleSource(dipoleSourceDLJ);
        pRF.setRange(truncation);
        pRF.setDielectric(dielectricOutside);
        
        
    	potentialMaster = new PotentialMaster();//   disable reaction field or dipole-dipole from here TODO???
        potentialMaster.addPotential(pDLJ, new ISpecies[] {species,species});
        potentialMaster.addPotential(pRF, new ISpecies[]{species, species});  
        potentialMaster.lrcMaster().addPotential(pRF.makeP0());
        
        // integrator from potential master
		integrator = new IntegratorMC(this, potentialMaster);
        // add mc move
        moveMolecule = new MCMoveMolecule(this, potentialMaster,space);//stepSize:1.0, stepSizeMax:15.0  ??????????????
        rotateMolecule = new MCMoveRotate(potentialMaster,random,space);

		activityIntegrate = new ActivityIntegrate(integrator);
		getController().addAction(activityIntegrate);

		//******************************** periodic boundary condition ******************************** //
		BoxImposePbc imposePbc = new BoxImposePbc(box, space);
		imposePbc.setApplyToMolecules(true);
		//**************************** integrator ****************************** //
		integrator.setTemperature(temperature);
	    integrator.setBox(box);
        integrator.getMoveManager().addMCMove(moveMolecule);
        integrator.getMoveManager().addMCMove(rotateMolecule);
		integrator.getEventManager().addListener(new IntegratorListenerAction(imposePbc));
		
		//******************************** initial configuration ******************************** //
		LatticeCubicBcc lattice = new LatticeCubicBcc(space);
		ConfigurationLattice configuration = new ConfigurationLattice(lattice, space);
		configuration.initializeCoordinates(box);
		
	}
		
	// **************************** simulation part **************************** //
	public static void main (String[] args){
		Param params = new Param();
		if (args.length > 0) {							
			ParseArgs.doParseArgs(params, args);
		} else {
			
		}
		final long startTime = System.currentTimeMillis();
		DateFormat date = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		Calendar cal = Calendar.getInstance();
		System.out.println("startTime : " + date.format(cal.getTime()));
		Space space = Space3D.getInstance();
		int steps = params.steps;
		boolean isGraphic = params.isGraphic;
		boolean mSquare = params.mSquare;
		boolean aEE = params.aEE;
		boolean sumtest = params.sumtest;
		
		double temperature = params.temperature;
	//	double temperatureK = params.temperatureK;
//		double temperature = Kelvin.UNIT.toSim(temperatureK);// convert Kelvin temperature to T(sim), essentially kT
		int numberMolecules = params.numberMolecules;
		double density = params.density;
		double sigmaLJ=params.sigmaLJ;
		double epsilonLJ=params.epsilonLJ;
		double dipoleStrength = Math.sqrt(params.dipoleStrength2);
		double dielectricOutside = params.dielectricOutside;
		//double densitySim = density * Constants.AVOGADRO * 1e-27;  // ?????? convert density to sim unit; in 1/(A)^3
		//System.out.println("Constants.AVOGADRO * 1e-27: "+Constants.AVOGADRO * 1e-27);
		double densitySim = density;
		double boxSize = Math.pow(numberMolecules/densitySim,(1.0/3.0)); 
		double truncation=boxSize* 0.49;
		

		System.out.println("******************* dipolar LJ, dielectric constant, NVT********************");
		System.out.println("number of molecules =\t"+numberMolecules);
		System.out.println("steps= "+ steps);
		System.out.println("density=\t"+density);
//		System.out.println("denisty(sim)="+densitySim);
		System.out.println("temperature=\t"+ temperature);
//		System.out.println("box size="+boxSize);
		System.out.println("truncation = "+truncation);
//		System.out.println("sigmaLJ="+sigmaLJ);
//		System.out.println("epsilonLJ="+epsilonLJ);
		System.out.println("dipoleStrength squared = "+params.dipoleStrength2);
//		System.out.println("dipoleStrength="+dipoleStrength);
//		System.out.println("dielectricOutside="+dielectricOutside);

		final DLJ_NVT_1site sim = new DLJ_NVT_1site(space,numberMolecules,sigmaLJ, epsilonLJ, dipoleStrength,
				dielectricOutside, boxSize,temperature,truncation);
		
    	if (isGraphic){
			SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
		    simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(PIXEL_SIZE));
	        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));	        
	        ((DiameterHashByType)((DisplayBox)simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species.getAtomType(0),1);
            ColorSchemeByType colorScheme = (ColorSchemeByType)simGraphic.getDisplayBox(sim.box).getColorScheme();
            colorScheme.setColor(sim.getSpecies(0).getAtomType(0), Color.red);
            OrientedSite[] sites = new OrientedSite[1];
            sites[0] = new OrientedSite(0.5, Color.BLUE, 0.2);
            ((DisplayBoxCanvasG3DSys)simGraphic.getDisplayBox(sim.box).canvas).setOrientationSites(
                    (IAtomTypeOriented)sim.getSpecies(0).getAtomType(0), sites);
            simGraphic.makeAndDisplayFrame(APP_NAME);
			simGraphic.getDisplayBox(sim.box).repaint();
	    	return ;
    	}
    	
    	sim.activityIntegrate.setMaxSteps(steps/5);// equilibration period
   		sim.getController().actionPerformed();
   		sim.getController().reset();
   		sim.integrator.getMoveManager().setEquilibrating(false);
//   		System.out.println("equilibration finished");
   		
   		int blockNumber = 1000;
   		int sampleAtInterval = numberMolecules;
   		int samplePerBlock = steps/sampleAtInterval/blockNumber;
//   		System.out.println("number of blocks is : "+blockNumber);
//   		System.out.println("sample per block is : "+samplePerBlock);
   		
   		MeterDipoleSumSquared1site dipoleSumSquaredMeter = null;
   		AccumulatorAverage dipoleSumSquaredAccumulator = null;
   		if(mSquare){
   			dipoleSumSquaredMeter = new MeterDipoleSumSquared1site(space,sim.box,dipoleStrength);  	
   			dipoleSumSquaredAccumulator = new AccumulatorAverageFixed(samplePerBlock);
   			DataPump dipolePump = new DataPump(dipoleSumSquaredMeter,dipoleSumSquaredAccumulator);
   			IntegratorListenerAction dipoleListener = new IntegratorListenerAction(dipolePump);
   			dipoleListener.setInterval(sampleAtInterval);
   			sim.integrator.getEventManager().addListener(dipoleListener);
   		}
   		
   		//AEE
   		DipoleSourceDLJ dipoleSourceDLJ = new DipoleSourceDLJ(space,dipoleStrength); 
   		MeterDipoleSumSquaredMappedAverage AEEMeter =  null;
   		AccumulatorAverageCovariance AEEAccumulator = null;
   		if(aEE){
         AEEMeter = new MeterDipoleSumSquaredMappedAverage(space, sim.box,sim, dipoleStrength,temperature,sim.potentialMaster);
		AEEMeter.setDipoleSource(dipoleSourceDLJ);
         AEEAccumulator = new AccumulatorAverageCovariance(samplePerBlock,true);
		DataPump AEEPump = new DataPump(AEEMeter,AEEAccumulator);
		IntegratorListenerAction AEEListener = new IntegratorListenerAction(AEEPump);
		
		
		AEEListener.setInterval(sampleAtInterval);
//		AEEListener.setInterval(1);// debug only
		
		sim.integrator.getEventManager().addListener(AEEListener);
        
   		}
		
        sim.activityIntegrate.setMaxSteps(steps);// equilibration period
        sim.getController().actionPerformed();
                
        //calculate dipoleSumSquared average
        double dipoleSumSquared = 0;
        double dipoleSumSquaredERR = 0;
        double dipoleSumCor = 0 ;
        if(mSquare){
        	dipoleSumSquared = ((DataDouble)((DataGroup)dipoleSumSquaredAccumulator.getData()).getData(dipoleSumSquaredAccumulator.AVERAGE.index)).x;
        	dipoleSumSquaredERR = ((DataDouble)((DataGroup)dipoleSumSquaredAccumulator.getData()).getData(dipoleSumSquaredAccumulator.ERROR.index)).x;
        	dipoleSumCor = ((DataDouble)((DataGroup)dipoleSumSquaredAccumulator.getData()).getData(dipoleSumSquaredAccumulator.BLOCK_CORRELATION.index)).x;
        }
        
        double AEE = 0;
        double AEEER =0;
        double AEECor = 0;
        if(aEE){
        double sum0 =  ((DataGroup)AEEAccumulator.getData()).getData(AEEAccumulator.AVERAGE.index).getValue(0); 
		double ERsum0 = ((DataGroup)AEEAccumulator.getData()).getData(AEEAccumulator.ERROR.index).getValue(0);
		AEECor = ((DataGroup)AEEAccumulator.getData()).getData(AEEAccumulator.BLOCK_CORRELATION.index).getValue(0);
		 AEE = sum0;
		 AEEER = ERsum0;
        }
		
        long endTime = System.currentTimeMillis();
        
        double totalTime = (endTime - startTime)/(1000.0*60.0);
        if(mSquare){
        System.out.println("-<M^2>*bt*bt:\t"+(-dipoleSumSquared/temperature/temperature)
        		+ " mSquareErr:\t" + (dipoleSumSquaredERR/temperature/temperature)
        		+ " mSquareDifficulty:\t"+(dipoleSumSquaredERR/temperature/temperature)*Math.sqrt(totalTime)
        		+ " dipolesumcor = " + dipoleSumCor );
        }
        if(aEE){
        	System.out.println("AEE_new:\t"+ (AEE) 
        		+ " AEEErr:\t" + AEEER 
        		+ " AEEDifficulty:\t"+ AEEER*Math.sqrt(totalTime)
        		+ " AEECor = " + AEECor );
        }
        
        if(sumtest){
        	double sum1 =  ((DataGroup)AEEAccumulator.getData()).getData(AEEAccumulator.AVERAGE.index).getValue(1);
        	double sum1Err =  ((DataGroup)AEEAccumulator.getData()).getData(AEEAccumulator.ERROR.index).getValue(1);
        	double sum1Cor = ((DataGroup)AEEAccumulator.getData()).getData(AEEAccumulator.BLOCK_CORRELATION.index).getValue(1);
//        	System.out.println("sum1:\t" + sum1 
//        		+ " sum1Err:\t" + sum1Err
//        		+ " sum1Difficulty:\t"+ sum1Err*Math.sqrt(totalTime)
//        		+ " sum1Cor:\t" + sum1Cor);
        	
        	IData covariance = ((DataGroup)AEEAccumulator.getData()).getData(AEEAccumulator.BLOCK_COVARIANCE.index);
    		covariance.getValue(1);
        	double ERsum0 = ((DataGroup)AEEAccumulator.getData()).getData(AEEAccumulator.ERROR.index).getValue(0);
        	double ERsum1 = ((DataGroup)AEEAccumulator.getData()).getData(AEEAccumulator.ERROR.index).getValue(1);
        	
        	double AEEtest = AEE + sum1*sum1;
    		double AEEerrtest = Math.sqrt(ERsum0*ERsum0 + 4*sum1*sum1*ERsum1*ERsum1 - 
    				2*ERsum0*sum1*2*ERsum1*covariance.getValue(1)/Math.sqrt(covariance.getValue(0)*covariance.getValue(3)));
        	
    		
        	System.out.println("AEEtest:\t" + AEEtest
        		+ " AEEerrtest:\t" + AEEerrtest );
        }
        System.out.println(  "Time taken (in mins): " + totalTime); 
	}
	
	// ******************* parameters **********************// 
	public static class Param extends ParameterBase {
		public boolean isGraphic = false;
		public boolean mSquare = false;
		public boolean aEE = true; 
		public boolean sumtest = true;
		public double temperature = 10;
		public int numberMolecules = 10;
		public double density = 0.325;
		public double dipoleStrength2 = 3;
		public double sigmaLJ = 1.0;
		public double epsilonLJ = 1.0;
		public double dielectricOutside = 1.0E11;
		public int steps = 1000000;
	}
}
