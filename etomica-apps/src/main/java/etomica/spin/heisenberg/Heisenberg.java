/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.spin.heisenberg;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.chem.elements.ElementSimple;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPump;
import etomica.data.IData;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxSpin2D;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.listener.IntegratorListenerAction;
import etomica.nbr.site.NeighborSiteManager;
import etomica.nbr.site.PotentialMasterSite;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.species.SpeciesSpheresMono;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.systems.LJ;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.RandomNumberGenerator;


/**
 * Compute  the dielectric constant of  2D Heisenberg model in conventional way and mapped averaging.
 * Conventional way: get the -bt*bt*<M^2> with M is the total dipole moment;
 * Mapped averaging: get the A_EE secondDerivative of free energy w.r.t electric field E when E is zero
 * Prototype for simulation of a more general magnetic system.
 *
 * @author David Kofke & Weisong Lin
 *
 */
public class Heisenberg extends Simulation {

	private static final String APP_NAME = "Heisenberg";


    /**
     * 2D heisenberg model in square lattice.
	 * @param space use to define vector/tensor
     * @param nCells total number of atoms = nCells*nCells
	 * @param interactionS  the J in heisenberg energy function: U = J*Cos(theta1-theta2)
	 * @param dipoleMagnitude is the strength of heisenberg dipole.
     */
    public Heisenberg(Space space, int nCells, double temperature, double interactionS, double dipoleMagnitude) {
		super(Space2D.getInstance());
		setRandom(new RandomNumberGenerator(1)); //debug only TODO



        potentialMaster = new PotentialMasterSite(this, nCells, space);
        box = new Box(space);
        addBox(box);
        int numAtoms = space.powerD(nCells);
        
        spins = new SpeciesSpheresRotating(space, new ElementSimple("A")); 
        
        addSpecies(spins);
        box.setNMolecules(spins, numAtoms);
        
        potential = new P2Spin(space,interactionS);
        // the electric field is zero for now, maybe need to add electric field in the future.
        field = new P1MagneticField(space,dipoleMagnitude);
        integrator = new IntegratorMC(this, potentialMaster);
        mcMove =new MCMoveRotate(potentialMaster, random, space);
        integrator.getMoveManager().addMCMove(mcMove);
        integrator.setTemperature(temperature);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        IAtomType type = spins.getLeafType();
//        potentialMaster.addPotential(field, new IAtomType[] {type});
        potentialMaster.addPotential(potential, new IAtomType[] {type, type});
        
        integrator.setBox(box);
        
    }

    private static final long serialVersionUID = 2L;
    public PotentialMasterSite potentialMaster; // difference betweet Pmaster pmastersite
    public IBox box;
    public SpeciesSpheresMono spins;
    public P2Spin potential;
    public P1MagneticField field;
    private IntegratorMC integrator;
    public MCMoveRotate mcMove;
    public final ActivityIntegrate activityIntegrate;
    
    
    public static void main(String[] args) {
    	Param params = new Param();
    	if (args.length > 0) {
    		ParseArgs.doParseArgs(params, args);
    	} else {

    	}
    	final long startTime = System.currentTimeMillis();
 		DateFormat date = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
 		Calendar cal = Calendar.getInstance();
 		System.out.println("startTime : " + date.format(cal.getTime()));
 		
    	boolean isGraphic = params.isGraphic;
    	double temperature = params.temperature;
        int nCells = params.nCells;
    	int numberMolecules = nCells*nCells;
    	boolean aEE = params.aEE;
    	boolean mSquare = params.mSquare;
    	int steps = params.steps;
    	double interactionS = params.interactionS;
    	double dipoleMagnitude = params.dipoleMagnitude;

    	System.out.println("steps= "+steps);
    	System.out.println("numberMolecules= "+numberMolecules);
    	System.out.println("temperature= " + temperature);

    	Space sp = Space2D.getInstance();
		Heisenberg sim = new Heisenberg(sp,nCells,temperature,interactionS,dipoleMagnitude);
    	
		MeterSpinMSquare meterMSquare = null;
		AccumulatorAverage dipoleSumSquaredAccumulator = null;

    	if(isGraphic){
    		meterMSquare   = new MeterSpinMSquare(sim.space,sim.box,dipoleMagnitude);
    		dipoleSumSquaredAccumulator = new AccumulatorAverageCollapsing(100);
    		DataPump dipolePump = new DataPump(meterMSquare,dipoleSumSquaredAccumulator);
    		IntegratorListenerAction dipoleListener = new IntegratorListenerAction(dipolePump);
    		dipoleListener.setInterval(10);
    		sim.integrator.getEventManager().addListener(dipoleListener);
    		
    		SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME, sim.space, sim.getController());
    		((SimulationRestart)simGraphic.getController().getReinitButton().getAction()).setConfiguration(null);
    		IAction repaintAction = simGraphic.getPaintAction(sim.box);
    		DisplayBox displayBox = simGraphic.getDisplayBox(sim.box);

    		simGraphic.remove(displayBox);
    		BoxAgentManager boxAgentManager = sim.potentialMaster.getCellAgentManager();
    		NeighborSiteManager neighborSiteManager = (NeighborSiteManager)boxAgentManager.getAgent(sim.box);
    		displayBox.setBoxCanvas(new DisplayBoxSpin2D(displayBox,neighborSiteManager, sim.space, sim.getController()));
    		simGraphic.add(displayBox);
    		DeviceSlider temperatureSlider = new DeviceSlider(sim.getController(), sim.integrator,"temperature");
    		temperatureSlider.setMinimum(0.5);
    		temperatureSlider.setMaximum(10.0);
    		temperatureSlider.setShowBorder(true);
    		LJ lj = new LJ();
    		temperatureSlider.setUnit(lj.temperature());
    		simGraphic.add(temperatureSlider);
    		temperatureSlider.setValue(sim.integrator.getTemperature());
    		DeviceSlider fieldSlider = new DeviceSlider(sim.getController(), sim.field, "h");
    		fieldSlider.setMinimum(-5.);
    		fieldSlider.setMaximum(+5.);
    		fieldSlider.setNMajor(5);
    		fieldSlider.setValue(0.0);
    		fieldSlider.setShowBorder(true);
    		fieldSlider.setLabel("Magnetic field");
    		simGraphic.add(fieldSlider);

    		DisplayTextBoxesCAE boxes = new DisplayTextBoxesCAE();
    		boxes.setAccumulator(dipoleSumSquaredAccumulator);
    		boxes.setLabel("Magnetization");
    		boxes.setLabelType(DisplayTextBox.LabelType.BORDER);
    		simGraphic.add(boxes);

    		simGraphic.getController().getReinitButton().setPostAction(repaintAction);
    		simGraphic.makeAndDisplayFrame(APP_NAME);
    	}//Graphics

		sim.activityIntegrate.setMaxSteps(steps/5);//TODO
		sim.getController().actionPerformed();
		sim.getController().reset();
    	int blockNumber = 100;
		int sampleAtInterval = numberMolecules;
		int samplePerBlock = steps/sampleAtInterval/blockNumber;

		
		System.out.println("equilibration finished");
        long equilibrationTime = System.currentTimeMillis();
		 System.out.println("equilibrationTime: " + (equilibrationTime - startTime)/(1000.0*60.0));
		
		//mSquare 
		if(mSquare){
			meterMSquare   = new MeterSpinMSquare(sim.space,sim.box,dipoleMagnitude);
			dipoleSumSquaredAccumulator = new AccumulatorAverageFixed(samplePerBlock);
			DataPump dipolePump = new DataPump(meterMSquare,dipoleSumSquaredAccumulator);
			IntegratorListenerAction dipoleListener = new IntegratorListenerAction(dipolePump);
			dipoleListener.setInterval(sampleAtInterval);
			sim.integrator.getEventManager().addListener(dipoleListener);
		}
		
		
		//AEE 
		MeterMappedAveraging AEEMeter = null;
		AccumulatorAverageCovariance AEEAccumulator = null;
		if(aEE){
			AEEMeter = new MeterMappedAveraging(sim.space, sim.box, sim,temperature,interactionS,dipoleMagnitude,sim.potentialMaster);
			AEEAccumulator = new AccumulatorAverageCovariance(samplePerBlock,true);
			DataPump AEEPump = new DataPump(AEEMeter,AEEAccumulator);
			IntegratorListenerAction AEEListener = new IntegratorListenerAction(AEEPump);
			AEEListener.setInterval(sampleAtInterval);
			sim.integrator.getEventManager().addListener(AEEListener);
		}
		
		sim.activityIntegrate.setMaxSteps(steps);	
		sim.getController().actionPerformed();

		
		//******************************** simulation start ******************************** //
		double dipoleSumSquared = 0 ;
		double dipoleSumSquaredERR = 0 ;
		double dipoleSumCor = 0;
		if(mSquare){
			 dipoleSumSquared = ((DataDouble)((DataGroup)dipoleSumSquaredAccumulator.getData()).getData(dipoleSumSquaredAccumulator.AVERAGE.index)).x;
			 dipoleSumSquaredERR = ((DataDouble)((DataGroup)dipoleSumSquaredAccumulator.getData()).getData(dipoleSumSquaredAccumulator.ERROR.index)).x;
			 dipoleSumCor = ((DataDouble)((DataGroup)dipoleSumSquaredAccumulator.getData()).getData(dipoleSumSquaredAccumulator.BLOCK_CORRELATION.index)).x;
		}
		
		
		double AEE = 0;
		 double AEEER = 0;
		 double AEECor = 0;
		if(aEE){
			AEE =  ((DataGroup)AEEAccumulator.getData()).getData(AEEAccumulator.AVERAGE.index).getValue(0); 
			AEEER = ((DataGroup)AEEAccumulator.getData()).getData(AEEAccumulator.ERROR.index).getValue(0);
			AEECor = ((DataGroup)AEEAccumulator.getData()).getData(AEEAccumulator.BLOCK_CORRELATION.index).getValue(0);
		IData covariance = ((DataGroup)AEEAccumulator.getData()).getData(AEEAccumulator.BLOCK_COVARIANCE.index);
		covariance.getValue(1);
		}
		
    	long endTime = System.currentTimeMillis();

    	double totalTime = (endTime - startTime)/(1000.0*60.0);
    	if(mSquare){
    		System.out.println("-<M^2>*bt*bt:\t"+(-dipoleSumSquared/temperature/temperature)
    				+ " mSquareErr:\t" + (dipoleSumSquaredERR/temperature/temperature)
    				+ " mSquareDifficulty:\t"+(dipoleSumSquaredERR/temperature/temperature)*Math.sqrt(totalTime)
    				+ " dipolesumCor= " + dipoleSumCor );
    		System.out.println(  "mSquare_Time: " + (endTime - startTime)/(1000.0*60.0)); 
    	}
    	
    	if(aEE){
    		System.out.println("AEE_new:\t"+ (AEE) 
    				+ " AEEErr:\t" + AEEER 
    				+ " AEEDifficulty:\t"+ AEEER*Math.sqrt(totalTime)
    				+ " AEECor= " + AEECor );
    		System.out.println(  "AEE_Time: " + (endTime - startTime)/(1000.0*60.0)); 
    	}
    	
    	

    }
    
    // ******************* parameters **********************//
    public static class Param extends ParameterBase {
    	public boolean isGraphic = false;
    	public boolean mSquare = true;
    	public boolean aEE = true; 
    	public double temperature = 10;// Kelvin
    	public int nCells = 1;//number of atoms is nCells*nCells
    	public double interactionS = 1;
    	public double dipoleMagnitude = 1;
    	public int steps = 1000000;
    }
}
