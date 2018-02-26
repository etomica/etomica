package etomica.rotation;

import java.awt.Color;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.*;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationFileBinary;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressure;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet.MyAgent;
import etomica.integrator.IntegratorVelocityVerletRattle;
import etomica.integrator.IntegratorVelocityVerletShake.BondConstraints;
import etomica.listener.IntegratorListenerAction;
import etomica.models.clathrates.MinimizationTIP4P.ChargeAgentSourceRPM;
import etomica.models.water.ConformationWaterTIP4P;
import etomica.models.water.SpeciesWater4P;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeAgentManager;
import etomica.molecule.MoleculePositionCOM;
import etomica.normalmode.*;
import etomica.potential.EwaldSummation;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.OrientationFull3D;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Calorie;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.units.Mole;
import etomica.units.Pixel;
import etomica.units.dimensions.Null;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;
/**
 @author  Weisong Lin
 */
public class WaterDADB extends Simulation {

	public IntegratorVelocityVerletRattle integrator;
	public PotentialMaster potentialMaster;
	public SpeciesWater4P species;
	public Box box;
	public ActivityIntegrate ai;
	protected Potential2SoftSphericalLS potentialLJLS;
	protected final MoleculeAgentManager latticeCoordinates;

	public WaterDADB(final Space space, double temperature, int numCells, double rCutRealES, double rCutLJ, boolean isIce, double kCut, double shakeTol, boolean unitCells) {
		super(space);
		setRandom(new RandomMersenneTwister(2));
//		if (precision ==1.0e-5 ){
//			precision_s = 3.047059472445871 ;
//		}	else if (precision == 5.0e-5){
//			precision_s = 2.800672811371045;}
//		else {throw new RuntimeException("improper precision value!");}
		box = new Box(space);
		addBox(box);
		species = new SpeciesWater4P(getSpace(), true);
		addSpecies(species);
		box.setNMolecules(species, 46*numCells*numCells*numCells);
		box.setDensity(46/12.03/12.03/12.03);
		ChargeAgentSourceRPM agentSource = new ChargeAgentSourceRPM(species, isIce);
		AtomLeafAgentManager<EwaldSummation.MyCharge> atomAgentManager = new AtomLeafAgentManager<EwaldSummation.MyCharge>(agentSource, box, EwaldSummation.MyCharge.class);
		int [] nC = new int[] {numCells, numCells, numCells};
		if(unitCells){
			numCells = 1;
		}
		ConfigurationFile config = new ConfigurationFile( numCells + "ncFinalPos");

		if(unitCells){
			ConfigurationFileBinary.replicate(config, box,nC ,space);
		}
		else{
			config.initializeCoordinates(box);
		}
		double a0 = box.getBoundary().getBoxSize().getX(0);
		double[] rC = new double[]{a0, a0, a0};
		double sigma , epsilon ; //TIP4P
		if(isIce){
			sigma = 3.1668; epsilon = Kelvin.UNIT.toSim(106.1);//TIP4P/Ice
		}else{//TIP4P
			double A = 600E3; // kcal A^12 / mol
			double C = 610.0; // kcal A^6 / mol
			double s6 = A/C;
			sigma = Math.pow(s6, 1.0/6.0);
			epsilon = Mole.UNIT.fromSim(Calorie.UNIT.toSim(C/s6*1000))/4.0;
		}
		latticeCoordinates = new MoleculeAgentManager(this,box,new MoleculeSiteSource(space, new MoleculePositionCOM(space), new WaterOrientationDefinition(space)));
        EwaldSummationLattice potentialES = new EwaldSummationLattice(box, atomAgentManager, space, kCut ,rCutRealES,latticeCoordinates);
		P2LennardJones potentialLJ = new P2LennardJones(space, sigma, epsilon);
		potentialLJLS = new Potential2SoftSphericalLS(space,rCutLJ,rC,potentialLJ,latticeCoordinates);
//		potentialLJ =  new P2SoftSphericalTruncated(space, potentialLJ, rC);
		potentialMaster = new PotentialMaster();
		potentialMaster.addPotential(potentialES, new AtomType[0]);
		potentialMaster.addPotential(potentialLJLS, new AtomType[]{species.getOxygenType(), species.getOxygenType()});


		double timeInterval = 0.002;
		int maxIterations = 100;
		integrator = new IntegratorVelocityVerletRattle(this, potentialMaster, space);
		integrator.setShakeTolerance(shakeTol);
		double lOH = ConformationWaterTIP4P.bondLengthOH;
		double lHH = Math.sqrt(2*lOH*lOH*(1-Math.cos(ConformationWaterTIP4P.angleHOH)));
		double lOM = ConformationWaterTIP4P.rOM;
		double lMH = Math.sqrt(lOH*lOH+lOM*lOM-2*lOH*lOM*Math.cos(0.5*ConformationWaterTIP4P.angleHOH));
		BondConstraints bondConstraints = new BondConstraints(new int[][]{{0,2},{1,2},{0,1}}, new double[]{lOH, lOH, lHH}){

			Vector vectorsum = space.makeVector();
			Vector ovector = space.makeVector();
			Vector centermass = space.makeVector();
			Vector h1vector = space.makeVector();
			Vector h2vector = space.makeVector();
			Vector mvector = space.makeVector();
			Vector newforce = space.makeVector();
			Vector newtorque = space.makeVector();

			public void redistributeForces(IMolecule molecule, AtomLeafAgentManager agentManager) {
				IAtomList leafList = molecule.getChildList();
				Vector h1 = leafList.getAtom(0).getPosition();
				Vector h2 = leafList.getAtom(1).getPosition();
				Vector o = leafList.getAtom(2).getPosition();
				Vector m = leafList.getAtom(3).getPosition();
				Vector h1torque = space.makeVector();
				Vector h2torque = space.makeVector();
				Vector h1h2 = space.makeVector();
				Vector om = space.makeVector();
				Vector otorque = space.makeVector();
				Vector mtorque = space.makeVector();
				Vector totaltorque = space.makeVector();
				Vector totalforce = space.makeVector();
				double hmass = leafList.getAtom(0).getType().getMass();
				double omass = leafList.getAtom(2).getType().getMass();
				double hmassPersent = hmass/(2*hmass + omass);
				double omassPersent = omass/(2*hmass + omass);
				centermass.Ea1Tv1(hmass, h1);
				centermass.PEa1Tv1(hmass, h2) ;
				centermass.PEa1Tv1(omass, o) ;
				centermass.TE(1/(2*hmass + omass));
				mvector.Ev1Mv2(m, centermass);
				h1vector.Ev1Mv2(h1, centermass);
				h2vector.Ev1Mv2(h2, centermass);
				ovector.Ev1Mv2(o, centermass);
				vectorsum.E(h1vector);
				vectorsum.PE(h2vector);
				vectorsum.PEa1Tv1(-2.0, ovector);
				Vector h1force = ((MyAgent)agentManager.getAgent(leafList.getAtom(0))).force();
				Vector h2force = ((MyAgent)agentManager.getAgent(leafList.getAtom(1))).force();
				Vector oforce = ((MyAgent)agentManager.getAgent(leafList.getAtom(2))).force();
				Vector mforce = ((MyAgent)agentManager.getAgent(leafList.getAtom(3))).force();

				//test for total force
//				totalforce.E(h1force);
//				totalforce.PE(h2force);
//				totalforce.PE(mforce);
//				totalforce.PE(oforce);
//				System.out.println(" totalforcebeforce = "+totalforce);
//				//test for total torque
//				h1torque.E(h1force);
//				h1torque.XE(h1vector);
//				h2torque.E(h2force);
//				h2torque.XE(h2vector);
//				otorque.E(oforce);
//				otorque.XE(ovector);
//				mtorque.E(mforce);
//				mtorque.XE(mvector);
//				totaltorque.E(h1torque);
//				totaltorque.PE(h2torque);
//				totaltorque.PE(otorque);
//				totaltorque.PE(mtorque);
//				System.out.println("totaltorque = "+totaltorque);




				//distribute  translation and rotation caused by mforce to h1 h2 & o
				newtorque.E(mforce);
				newtorque.XE(mvector);
				newforce.E(0);
				if(Math.abs(vectorsum.getX(0)) > 0.1){
					newforce.setX(1, -newtorque.getX(2)/vectorsum.getX(0));
					newforce.setX(2, newtorque.getX(1)/vectorsum.getX(0));
				}
				else if(Math.abs(vectorsum.getX(1)) > 0.1){

					newforce.setX(0, newtorque.getX(2)/vectorsum.getX(1));
					newforce.setX(2, -newtorque.getX(0)/vectorsum.getX(1));
				}
				else {
					newforce.setX(0, -newtorque.getX(1)/vectorsum.getX(2));
					newforce.setX(1, newtorque.getX(0)/vectorsum.getX(2));
				}

				h1force.PE(newforce);
				h2force.PE(newforce);
				oforce.PEa1Tv1(-2.0, newforce);
				h1force.PEa1Tv1(hmassPersent, mforce);
				h2force.PEa1Tv1(hmassPersent, mforce);
				oforce.PEa1Tv1(omassPersent, mforce);
				mforce.ME(mforce);

				//redistribute the forces s.t. only translation
//				totalforce.E(h1force);
//				totalforce.PE(h2force);
//				totalforce.PE(mforce);
//				totalforce.PE(oforce);
//				h1force.E(0);
//				h2force.E(0);
//				oforce.E(0);
//				mforce.E(0);
//				h1force.PEa1Tv1(hmassPersent, totalforce);
//				h2force.PEa1Tv1(hmassPersent, totalforce);
//				oforce.PEa1Tv1(omassPersent, totalforce);

				//redistribute s.t. only rotation
//				totalforce.E(h1force);
//				totalforce.PE(h2force);
//				totalforce.PE(mforce);
//				totalforce.PE(oforce);
//				h1force.PEa1Tv1(-1.0*hmassPersent, totalforce);
//				h2force.PEa1Tv1(-1.0*hmassPersent, totalforce);
//				oforce.PEa1Tv1(-1.0*omassPersent, totalforce);

				//redistribute such that only K3 freedom
//				oforce.E(0);
//				h1h2.Ev1Mv2(h2, h1);
//				om.Ev1Mv2(m, o);
//				h1h2.normalize();
//				om.normalize();
//				h1force.PEa1Tv1(-1.0*h1force.dot(om),om);
//				h1force.PEa1Tv1(-1.0*h1force.dot(h1h2), h1h2);
//				h2force.PEa1Tv1(-1.0*h2force.dot(om),om);
//				h2force.PEa1Tv1(-1.0*h2force.dot(h1h2), h1h2);
//				totalforce.E(h1force);
//				totalforce.PE(h2force);
//				h1force.PEa1Tv1(-0.5, totalforce);
//				h2force.PEa1Tv1(-0.5, totalforce);


				//redistribute s.t. only k1 freedom with only in h1h2 direction
//				h1h2.Ev1Mv2(h2, h1);
//				om.Ev1Mv2(m, o);
//				h1h2.XE(om);
//				h1h2.normalize();
//				h1force.PEa1Tv1(-1.0*h1force.dot(h1h2),h1h2);
//				h2force.PEa1Tv1(-1.0*h2force.dot(h1h2),h1h2);
//				oforce.PEa1Tv1(-1.0*oforce.dot(h1h2), h1h2);

				//redistribute s.t. only k1&k2 freedom
//				h1h2.Ev1Mv2(h2, h1);
//				om.Ev1Mv2(m, o);
//				h1h2.XE(om);
//				h1h2.normalize();
//				h1force.PEa1Tv1(-1.0*h1force.dot(h1h2),h1h2);
//				h2force.PEa1Tv1(-1.0*h2force.dot(h1h2),h1h2);
//				oforce.PEa1Tv1(-1.0*h2force.dot(h1h2), o);
//


				//test for total force
//				totalforce.E(h1force);
//				totalforce.PE(h2force);
//				totalforce.PE(mforce);
//				totalforce.PE(oforce);
//				System.out.println(" retributiedtotalforce = " + totalforce);

				//test for total torque
//				h1torque.E(h1force);
//				h1torque.XE(h1vector);
//				h2torque.E(h2force);
//				h2torque.XE(h2vector);
//				otorque.E(oforce);
//				otorque.XE(ovector);
//				mtorque.E(mforce);
//				mtorque.XE(mvector);
//				totaltorque.E(h1torque);
//				totaltorque.PE(h2torque);
//				totaltorque.PE(mtorque);
//				totaltorque.PE(otorque);
//				System.out.println("totaltorque = "+totaltorque);
//				System.(2);

			}

			public void relaxMolecule(IMolecule molecule) {
				IAtomList leafList = molecule.getChildList();

				Vector h1 = leafList.getAtom(0).getPosition();
				Vector h2 = leafList.getAtom(1).getPosition();
				Vector o = leafList.getAtom(2).getPosition();
				Vector m = leafList.getAtom(3).getPosition();
				m.E(h1);
				m.PE(h2);
				m.PEa1Tv1(-2, o);
				m.TE(0.15/Math.sqrt(m.squared())); // TIP4P
//            	m.TE(0.1577/Math.sqrt(m.squared())); TIP4P/Ice
				m.PE(o);

			}

		};
		integrator.getShakeAgentManager().setAgent(species, bondConstraints);
		integrator.setTimeStep(timeInterval);
		integrator.setMaxIterations(maxIterations);
		integrator.setBox(box);
//        integrator.setOrientAtom((IAtom)((IMolecule)box.getMoleculeList(speciesOrient).getAtom(0)).getChildList().getAtom(0));
		integrator.setIsothermal(true);
		integrator.setTemperature(Kelvin.UNIT.toSim(temperature));
		integrator.setThermostatInterval(100);
		integrator.setThermostatNoDrift(true);
		try {
			integrator.reset();
		}
		catch (ConfigurationOverlapException e){}
		ai = new ActivityIntegrate(integrator);
//        System.out.println("using rigid with dt="+dt);
		getController().addAction(ai);
//        System.out.println("h1 at "+((IAtomPositioned)box.getLeafList().getAtom(0)).getPosition());
//        System.out.println("o at "+((IAtomPositioned)box.getLeafList().getAtom(2)).getPosition());

	}

	public static void main(String[] args) {
//    	System.out.println(Kelvin.UNIT.toSim(1));
//    	System.out.println(Joule.UNIT.toSim(1)/Constants.AVOGADRO);
//    	System.exit(2);
		final long startTime = System.currentTimeMillis();
		WaterDADBParam waterDADBParam = new WaterDADBParam();
		ParseArgs.doParseArgs(waterDADBParam, args);
		final double temperature = waterDADBParam.temperature;
		int numCells = waterDADBParam.numCells;
		int numSteps = waterDADBParam.numSteps;
		double rCutRealES = waterDADBParam.rCutRealES;
		double rCutLJ = waterDADBParam.rCutLJ;
		boolean isIce = waterDADBParam.isIce;
		double kCut = waterDADBParam.kCut;
		double shakeTol = waterDADBParam.shakeTol;
		boolean uniteCells = waterDADBParam.unitCells;
		final WaterDADB sim = new WaterDADB(Space3D.getInstance(),temperature,numCells,rCutRealES,rCutLJ,isIce,kCut,shakeTol,uniteCells);
		MeterPotentialEnergy meterPE2 = new MeterPotentialEnergy(sim.potentialMaster);
		meterPE2.setBox(sim.box);
		final double latticeEnergy;
		latticeEnergy = meterPE2.getDataAsScalar();
		System.out.println("latticeEnergy = " + latticeEnergy);

//      try{
//      	EwaldSummation.fileWriter.close();
//  	}
//  	catch (IOException e){
//  		throw new RuntimeException(e);
//  	}
//      System.exit(2);



		if (false) {
			SpeciesSpheresMono guestSpecies = new SpeciesSpheresMono(sim,sim.space);
			sim.addSpecies(guestSpecies);
			sim.box.setNMolecules(guestSpecies, 8);
//        	sim.box.getMoleculeList(guestSpecies).getMolecule(0).getChildList().getAtom(0).getPosition().E(new double [] {6, 6, 6});
//        	sim.box.getMoleculeList(guestSpecies).getMolecule(1).getChildList().getAtom(0).getPosition().E(new double [] {6, -6, -6});
//        	sim.box.getMoleculeList(guestSpecies).getMolecule(2).getChildList().getAtom(0).getPosition().E(new double [] {6, -6, 6});
//        	sim.box.getMoleculeList(guestSpecies).getMolecule(3).getChildList().getAtom(0).getPosition().E(new double [] {6, 6, -6});
//        	sim.box.getMoleculeList(guestSpecies).getMolecule(4).getChildList().getAtom(0).getPosition().E(new double [] {-6, 6, 6});
//        	sim.box.getMoleculeList(guestSpecies).getMolecule(5).getChildList().getAtom(0).getPosition().E(new double [] {-6, -6, -6});
//        	sim.box.getMoleculeList(guestSpecies).getMolecule(6).getChildList().getAtom(0).getPosition().E(new double [] {-6, 6, -6});
//        	sim.box.getMoleculeList(guestSpecies).getMolecule(7).getChildList().getAtom(0).getPosition().E(new double [] {-6, -6, 6});
			sim.box.getMoleculeList(guestSpecies).getMolecule(0).getChildList().getAtom(0).getPosition().E(new double [] {0, 0, 0});

			sim.ai.setSleepPeriod(2);
			SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "Rattle", 1, sim.space, sim.getController());
			((ColorSchemeByType)graphic.getDisplayBox(sim.box).getColorScheme()).setColor(sim.species.getHydrogenType(), Color.WHITE);
			((ColorSchemeByType)graphic.getDisplayBox(sim.box).getColorScheme()).setColor(sim.species.getOxygenType(), Color.RED);
			((DiameterHashByType)graphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species.getMType(), 0.1);
			((DiameterHashByType)graphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(guestSpecies.getAtomType(0), 2.5);

			((ColorSchemeByType)graphic.getDisplayBox(sim.box).getColorScheme()).setColor(guestSpecies.getAtomType(0), Color.ORANGE);

			((DisplayBoxCanvasG3DSys)graphic.getDisplayBox(sim.box).canvas).setBackgroundColor(Color.WHITE);
			((DisplayBoxCanvasG3DSys)graphic.getDisplayBox(sim.box).canvas).setBoundaryFrameColor(Color.black);
			//((DiameterHashByType)graphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species.getOxygenType(), .1);
//            ((DiameterHashByType)graphic.getDisplayBox(box).getDiameterHash()).setDiameter(hType, 1);
//            MeterEnergy meterE = new MeterEnergy(sim.potentialMaster, sim.box);
			MeterKineticEnergy meterE = new MeterKineticEnergy();
			meterE.setBox(sim.box);
			AccumulatorHistory history = new AccumulatorHistory(new HistoryCollapsingAverage());
			history.setTimeDataSource(new DataSourceCountTime(sim.integrator));
			DataProcessor processor = new DataProcessor() {
				protected DataInfoDouble dataInfo;
				protected DataDouble data = new DataDouble();
				protected DataTag tag = new DataTag();
				public DataPipe getDataCaster(IDataInfo inputDataInfo) {
					return null;
				}
				protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
					dataInfo = new DataInfoDouble("foo", Null.DIMENSION);
					dataInfo.addTags(inputDataInfo.getTags());
					dataInfo.addTag(tag);
					return dataInfo;
				}
				protected IData processData(IData inputData) {
					data.x = inputData.getValue(0) /temperature;
					return data;
				}
			};
			DataPump pump = new DataPump(meterE, processor);
			processor.setDataSink(history);
			DisplayPlot ePlot = new DisplayPlot();
//            ePlot.setUnit(new SimpleUnit(Energy.DIMENSION,46*Joule.UNIT.toSim(1)/Constants.AVOGADRO*1000,"energy","symbol",true));
			history.setDataSink(ePlot.getDataSet().makeDataSink());
			IntegratorListenerAction pumpListener = new IntegratorListenerAction(pump);
			pumpListener.setInterval(10);
			sim.integrator.getEventManager().addListener(pumpListener);
			ePlot.setLabel("Energy");
			graphic.add(ePlot);
			graphic.getDisplayBox(graphic.getSimulation().getBox(0)).setPixelUnit(new Pixel(25));
			graphic.makeAndDisplayFrame();
			System.out.println(meterE.getDataAsScalar()/(46*Joule.UNIT.toSim(1)/Constants.AVOGADRO*1000));
			System.out.println(meterE.getDataAsScalar());
			return ;
		}
		final MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster);
		meterPE.setBox(sim.box);

		int blockSize = numSteps >= 1000? (numSteps/1000):1;

		MeterDADBWaterTIP4P meterDADB = new MeterDADBWaterTIP4P(sim.space,meterPE,sim.potentialMaster, Kelvin.UNIT.toSim(temperature),sim.latticeCoordinates);
//        MeterDADB.justU = true;
		MeterPotentialEnergy meterPotentialEnergy = new MeterPotentialEnergy(sim.potentialMaster);

		MeterPressure meterPressure = new MeterPressure(sim.space);
		meterPressure.setBox(sim.box);
		meterPressure.setPotentialMaster(sim.potentialMaster);
		meterPressure.setTemperature(temperature);


		meterPotentialEnergy.setBox(sim.box);


		AccumulatorAverageFixed accumulatorAverageFixed4 = new AccumulatorAverageFixed(blockSize);
		DataPumpListener dataPumpListener4 = new DataPumpListener(meterPotentialEnergy, accumulatorAverageFixed4, 10);

		AccumulatorAverageFixed accumulatorAverageFixed2 = new AccumulatorAverageFixed(blockSize);
		DataPumpListener dataPumpListener2 = new DataPumpListener(meterPotentialEnergy, accumulatorAverageFixed2, 10);
		MeterKineticEnergy meterKineticEnergy = new MeterKineticEnergy();
		meterKineticEnergy.setBox(sim.box);
		AccumulatorAverageFixed accumulatorAverageFixed3 = new AccumulatorAverageFixed(blockSize);
		DataPumpListener dataPumpListener3 = new DataPumpListener(meterKineticEnergy, accumulatorAverageFixed3, 10);
		AccumulatorAverageFixed accumulatorAverageFixed = new AccumulatorAverageFixed(blockSize);
		DataPumpListener dataPumpListener = new DataPumpListener(meterDADB, accumulatorAverageFixed, 10);

		sim.ai.setMaxSteps(numSteps/10);
		sim.getController().actionPerformed();
//        meterDADB.debug();
//        System.exit(2);

		sim.ai.setMaxSteps(numSteps);


        sim.integrator.getEventManager().addListener(dataPumpListener);
		sim.integrator.getEventManager().addListener(dataPumpListener2);
//        sim.integrator.getEventManager().addListener(dataPumpListener3);


//        sim.integrator.getEventManager().addListener(dataPumpListener4);


		sim.getController().reset();
		sim.getController().actionPerformed();
//        double KEaverage = accumulatorAverageFixed3.getData(accumulatorAverageFixed3.AVERAGE).getValue(0);
//        double KEerror = accumulatorAverageFixed3.getData(accumulatorAverageFixed3.ERROR).getValue(0);
//        double KEcorrelation = accumulatorAverageFixed3.getData(accumulatorAverageFixed3.BLOCK_CORRELATION).getValue(0);
//        System.out.println("KEaverage =" +  KEaverage);
//        System.out.println("KEerror =" + KEerror);
//        System.out.println("KEcorrelation = " + KEcorrelation);


//        double Paverage = accumulatorAverageFixed4.getData(accumulatorAverageFixed4.AVERAGE).getValue(0);
//        double Perror = accumulatorAverageFixed4.getData(accumulatorAverageFixed4.ERROR).getValue(0);
//        double Pcorrelation = accumulatorAverageFixed4.getData(accumulatorAverageFixed4.BLOCK_CORRELATION).getValue(0);


		long endTime = System.currentTimeMillis();
		DateFormat date = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		Calendar cal = Calendar.getInstance();
		System.out.println(date.format(cal.getTime()));
		System.out.println("Time taken (in mins): " + (endTime - startTime)/(1000.0*60.0));
		System.out.println("numSteps = " + numSteps);

//
//        System.out.println("Paverage =" +  Paverage);
//        System.out.println("Perror =" + Perror);
//        System.out.println("Pcorrelation = " + Pcorrelation);


		double average = accumulatorAverageFixed.getData(accumulatorAverageFixed.AVERAGE).getValue(0);
		double error = accumulatorAverageFixed.getData(accumulatorAverageFixed.ERROR).getValue(0);
		double correlation = accumulatorAverageFixed.getData(accumulatorAverageFixed.BLOCK_CORRELATION).getValue(0);

//        System.out.println("rCutLJ = "  + rCutLJ);
//        System.out.println("rCutRealES = " + rCutRealES);
		System.out.println("temperature = " + temperature );
//        System.out.println("hh and om are made perpendicular" );
//        System.out.println("shakeTol = " + shakeTol);
		IMoleculeList molecules = sim.box.getMoleculeList();
//        System.out.println("for only translation case");
		System.out.println("beginLE = " + latticeEnergy);
		System.out.println("1.5*nkT = " + Kelvin.UNIT.toSim((molecules.getMoleculeCount()*1.5*temperature)));
		System.out.println("average = " +  average);
		System.out.println("error = " + error);
		System.out.println("correlation = " + correlation);
		double PEAverage = accumulatorAverageFixed2.getData(accumulatorAverageFixed2.AVERAGE).getValue(0);
		double PEArror = accumulatorAverageFixed2.getData(accumulatorAverageFixed2.ERROR).getValue(0);
		double PECorrelation = accumulatorAverageFixed2.getData(accumulatorAverageFixed2.BLOCK_CORRELATION).getValue(0);
        System.out.println("PEaverage = " + (PEAverage - latticeEnergy - Kelvin.UNIT.toSim(((molecules.getMoleculeCount() - 1) * 3 * temperature))));
		System.out.println("PEerror = " + PEArror);
		System.out.println("PEcorrelation = " + PECorrelation);

		ConfigurationFile config = new ConfigurationFile( numCells + "ncFinalPos");
		config.initializeCoordinates(sim.box);
		final double endLatticeEnergy = meterPE2.getDataAsScalar();
		System.out.println("endLE = " + endLatticeEnergy);


	}

	public static class WaterOrientationDefinition implements etomica.normalmode.MoleculeSiteSource.MoleculeOrientationDefinition{
		protected final OrientationFull3D or;
		protected final Vector v1,v2;

		public WaterOrientationDefinition(Space space){
			or = new OrientationFull3D(space);
			v1 = space.makeVector();
			v2 = space.makeVector();

		}
		public IOrientation getOrientation(IMolecule molecule) {
			IAtomList leafList = molecule.getChildList();
			Vector h1 = leafList.getAtom(0).getPosition();
			Vector h2 = leafList.getAtom(1).getPosition();
			Vector o = leafList.getAtom(2).getPosition();
			Vector m = leafList.getAtom(3).getPosition();
			v1.Ev1Mv2(m, o);
			v1.normalize();
			v2.Ev1Mv2(h2, h1);
			v2.normalize();
			or.setDirections( v1,  v2);
			//v1 is a0 and v2 is a
			return or;
		}
	}

	public static class WaterDADBParam extends ParameterBase {
		public int numCells = 1;
		public int numSteps = 100000;
		public double temperature = 100;
		public double rCutLJ = 11;
		public double rCutRealES = 11;
		public double kCut = 1.5;
		public boolean isIce =  false;
		public double shakeTol = 1e-12;
		public boolean unitCells = false;

	}
}
