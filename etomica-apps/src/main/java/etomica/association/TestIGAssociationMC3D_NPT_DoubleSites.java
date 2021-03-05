/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.action.BoxInflate;
import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.*;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2HardAssociationConeSW;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.Degree;
import etomica.util.ParameterBase;

/**
 * Simple ideal gas + S-W Association Monte Carlo NPT simulation in 3D.
 * average density = N*<1/V>
 * Initial configurations at http://rheneas.eng.buffalo.edu/etomica/tests/
 */
public class TestIGAssociationMC3D_NPT_DoubleSites extends Simulation {
    
    public IntegratorMC integrator;
    public MCMoveAtomMonomer mcMoveAtomMonomer;
    public MCMoveAtomSmer mcMoveAtomSmer;
    public MCMoveRotateAssociated mcMoveRotate;
    public SpeciesGeneral species;
    public Box box;
    public P2HardAssociationConeSW potential;
    public MCMoveSmer mcMoveSmer;
    public MCMoveSmerRotate mcMoveSmerRotate;
    public MCMoveVolumeAssociated mcMoveVolume;

    public MCMoveBiasUB mcMoveBiasUB;
    public AssociationManager associationManagerOriented;
    public BiasVolumeSphereOrientedDoubleSites bvso;
    public AssociationHelperDouble associationHelper;
    double epsilon = 1.0;
        
    
    public TestIGAssociationMC3D_NPT_DoubleSites(int numAtoms, double pressure, double density, double wellConstant, double temperature,double truncationRadius,int maxChainLength, boolean useUB, long numSteps) {
        super(Space3D.getInstance());

        species = SpeciesSpheresRotating.create(space, new ElementSimple(this), false,true);//Species in which molecules are made of a single atom of type OrientedSphere
        addSpecies(species);

        PotentialMasterCell potentialMaster = new PotentialMasterCell(this);

        double sigma = 1.0;
        //setRandom(new RandomNumberGenerator(3));

        System.out.println("Double site Association");
        System.out.println("pressure = " + pressure);
        System.out.println("initial density = " + density);
        System.out.println("association strength = " + wellConstant + "*epsilon");
        System.out.println("temperature = " + temperature);
        System.out.println("numSteps = " + numSteps);
        System.out.println("maximum chain length= " + maxChainLength);
        box = this.makeBox();
        integrator = new IntegratorMC(this.getRandom(), potentialMaster, box);
        integrator.setTemperature(temperature);
        mcMoveAtomMonomer = new MCMoveAtomMonomer(random, potentialMaster, space);//Standard Monte Carlo atom-displacement trial move
        mcMoveAtomMonomer.setMaxLength(maxChainLength);
        mcMoveAtomSmer = new MCMoveAtomSmer(random, potentialMaster, space);
        mcMoveAtomSmer.setMaxLength(maxChainLength);
        mcMoveRotate = new MCMoveRotateAssociated(potentialMaster, random, space);//Performs a rotation of an atom (not a molecule) that has an orientation coordinate
        mcMoveRotate.setMaxLength(maxChainLength);
        bvso = new BiasVolumeSphereOrientedDoubleSites(space, random);
        System.out.println("biasVolume= 2*(partialVolume/totalVolume)^2");
        bvso.setTheta(Degree.UNIT.toSim(27.0));
        bvso.setBiasSphereInnerRadius(0.8);
        bvso.setBox(box);
        associationManagerOriented = new AssociationManager(box, potentialMaster, bvso);//define and track atom associations
        associationHelper = new AssociationHelperDouble(space, box, associationManagerOriented);
        box.setNMolecules(species, numAtoms);
        mcMoveBiasUB = new MCMoveBiasUB(potentialMaster, bvso, random, space);
        mcMoveBiasUB.setMaxLength(maxChainLength);//only allow the formation up to maxChainLengh-mer
        mcMoveAtomMonomer.setAssociationManager(associationManagerOriented);
        mcMoveAtomSmer.setAssociationManager(associationManagerOriented);
        mcMoveRotate.setAssociationManager(associationManagerOriented);
        mcMoveBiasUB.setAssociationManager(associationManagerOriented);

        //mcMoveAtom.setStepSize(0.2*sigma);
        ((MCMoveStepTracker) mcMoveAtomMonomer.getTracker()).setNoisyAdjustment(true);
        ((MCMoveStepTracker) mcMoveAtomSmer.getTracker()).setNoisyAdjustment(true);
        ((MCMoveStepTracker) mcMoveRotate.getTracker()).setNoisyAdjustment(true);
        integrator.getMoveManager().addMCMove(mcMoveAtomMonomer);
        integrator.getMoveManager().addMCMove(mcMoveAtomSmer);
        integrator.getMoveManager().addMCMove(mcMoveRotate);
        if (useUB) {
            integrator.getMoveManager().addMCMove(mcMoveBiasUB);
            System.out.println("with BiasUB");
        } else {
            System.out.println("without BiasUB");
        }
        integrator.getMoveEventManager().addListener(associationManagerOriented);
        integrator.getMoveManager().setEquilibrating(true);
        this.getController().addActivity(new ActivityIntegrate(integrator), numSteps);
        //actionIntegrate.setSleepPeriod(1);
        BoxInflate inflater = new BoxInflate(box, space);//Performs actions that cause volume of system to expand or contract
        inflater.setTargetDensity(density);
        inflater.actionPerformed();

        System.out.println("truncation distance of potential = " + truncationRadius);
        if (truncationRadius > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large.  Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        potential = new P2HardAssociationConeSW(space, sigma, epsilon, truncationRadius, wellConstant);
        potential.setInnerWellCutoffFactor(0.8);
        potentialMaster.setCellRange(3);
        potentialMaster.setRange(potential.getRange());
        mcMoveSmer = new MCMoveSmer(random, potentialMaster, space, potential);
        mcMoveSmerRotate = new MCMoveSmerRotate(random, potentialMaster, space, potential);
        mcMoveVolume = new MCMoveVolumeAssociated(random, potentialMaster, space);
        //MCMoveVolumeAssociated.dodebug =false;
        mcMoveVolume.setAssociationManager(associationManagerOriented);
        mcMoveSmer.setAssociationManager(associationManagerOriented);
        mcMoveSmerRotate.setAssociationManager(associationManagerOriented);
        mcMoveVolume.setPressure(pressure);

        AtomType leafType = species.getLeafType();
        potentialMaster.addPotential(potential, new AtomType[]{leafType, leafType});
        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());
        integrator.getMoveManager().addMCMove(mcMoveSmer);
        integrator.getMoveManager().addMCMove(mcMoveSmerRotate);
        integrator.getMoveManager().addMCMove(mcMoveVolume);

        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box);
        associationManagerOriented.initialize();
        potentialMaster.getNbrCellManager(box).assignCellAll();
        potentialMaster.getNbrCellManager(box).setDoApplyPBC(true);
    }
 
    public static void main(String[] args) {
    	VirialAssociatingFluidParam params = new VirialAssociatingFluidParam();
    	
    	int numAtoms = params.numAtoms;
    	double pressure = params.pressure;
    	double density = params.density;
    	double wellConstant = params.wellConstant;
        double temperature = params.temperature;
        double truncationRadius = params.truncationRadius;//truncation distance of potential, default = 3.0*sigma
        int maxChainLength = params.maxChainLength;
        boolean useUB = params.useUB;
        long numSteps = params.numSteps;
        if (args.length > 0) {
            numAtoms = Integer.parseInt(args[0]);
            pressure = Double.parseDouble(args[1]);
            density = Double.parseDouble(args[2]);
            wellConstant = Double.parseDouble(args[3]);
            temperature = Double.parseDouble(args[4]);
            truncationRadius = Double.parseDouble(args[5]);
            useUB = Integer.parseInt(args[6])!=0;
            numSteps = Long.parseLong(args[7]);
            
        }
        final TestIGAssociationMC3D_NPT_DoubleSites sim = new TestIGAssociationMC3D_NPT_DoubleSites(numAtoms, pressure, density, wellConstant, temperature, truncationRadius,maxChainLength, useUB, numSteps);
        IAction energyDiffActionEq = new IAction() {
            protected final Vector dr = sim.space.makeVector();

            public void actionPerformed() {
//				if (sim.integrator.getStepCount() > 25300){
//					IntegratorMC.dodebug = true;
//				}
				if (sim.integrator.getStepCount()%100 == 0){
					IAtomList leafList = sim.box.getLeafList();
					for (int i = 0; i <leafList.size(); i+=1){
						IAtomOriented atomi = (IAtomOriented)sim.box.getLeafList().get(i);
						AtomArrayList bondList = (AtomArrayList)sim.associationManagerOriented.getAssociatedAtoms(atomi);
						if (bondList.size() == 2){
				    		IAtom atom0 = bondList.get(0);
				    		IAtom atom1 = bondList.get(1);
				    		double innerRadius = 0.8;
				        	double minDistance = 2*(innerRadius*innerRadius)*(1+Math.cos(Degree.UNIT.toSim(27.0)));
				    		dr.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());//dr = distance from the atom0 to atom1
				        	sim.box.getBoundary().nearestImage(dr);
				        	if (dr.squared() < minDistance){
				        		System.out.println(sim.integrator.getStepCount()+ " steps");
				        		System.out.println("atomi "+atomi+ "bondList "+bondList);
				        		System.out.println("dr*dr "+dr.squared());
								throw new RuntimeException();
				        	}
				    	}
						if (bondList.size() > 2) {
							System.out.println(sim.integrator.getStepCount()+ " steps");
							System.out.println("atomi "+atomi+ "bondList "+bondList);
							throw new RuntimeException();
						}
						for (int j =0; j<i;j+=1){
							IAtomOriented atomj = (IAtomOriented)sim.box.getLeafList().get(j);
							boolean isBonded1 = bondList.indexOf(atomj)>-1;
							boolean isBonded2 = sim.bvso.isAssociated(atomi, atomj);
							if (isBonded1 != isBonded2){
								System.out.println(sim.integrator.getStepCount()+ " steps");
								System.out.println("Wrong. isBonded1: "+isBonded1+" isBonded2: "+isBonded2);
								System.out.println("atomi= "+atomi+" atomj= "+atomj);
								System.out.println("position of atomi: "+atomi.getPosition()+" orientation of atomi: "+atomi.getOrientation().getDirection());
								System.out.println("position of atomj: "+atomj.getPosition()+" orientation of atomj: "+atomj.getOrientation().getDirection());
								throw new RuntimeException();
							}
						}
					}
				}


			}
		};
        IntegratorListenerAction energyDiffListenerEq = new IntegratorListenerAction(energyDiffActionEq,1);
        sim.integrator.getEventManager().addListener(energyDiffListenerEq);
        System.out.println("equilibrium period = " +numSteps/5);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps / 5));
System.out.println("equilibrium finished");

MeterDensity rhoMeter = new MeterDensity(sim.box);
        AccumulatorAverage rhoAccumulator = new AccumulatorAverageFixed(10);//Accumulator that keeps statistics for averaging and error analysis
        DataPump rhoPump = new DataPump(rhoMeter,rhoAccumulator);
        IntegratorListenerAction listener = new IntegratorListenerAction(rhoPump);
        listener.setInterval(50);
        sim.integrator.getEventManager().addListener(listener);
        final MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed(10);
        DataPump energyManager = new DataPump(energyMeter, energyAccumulator);
        energyAccumulator.setBlockSize(50);
        IntegratorListenerAction energyListener = new IntegratorListenerAction(energyManager);
        energyListener.setInterval(1000);
        sim.integrator.getEventManager().addListener(energyListener);

        AccumulatorAverage smerAccumulator = new AccumulatorAverageFixed(10);
        MeterDimerMoleFraction smerMeter = new MeterDimerMoleFraction(sim.getSpace(), sim.box);
        DataPump dimerPump = new DataPump(smerMeter,smerAccumulator);
        IntegratorListenerAction smerListener = new IntegratorListenerAction(dimerPump);
        smerListener.setInterval(50);
        sim.integrator.getEventManager().addListener(smerListener);
        smerMeter.setAssociationManager(sim.associationManagerOriented);

        AccumulatorAverage energy2Accumulator = new AccumulatorAverageFixed(10);
        final MeterPotentialEnergy energy2Meter = new MeterPotentialEnergy(sim.integrator.getPotentialMaster());//true energy
        energy2Meter.setBox(sim.box);
        DataPump energy2Pump = new DataPump(energy2Meter,energy2Accumulator);
        IntegratorListenerAction energy2Listener = new IntegratorListenerAction(energy2Pump);
        energy2Listener.setInterval(1000);
        sim.integrator.getEventManager().addListener(energy2Listener);

        if (false) {
        	SimulationGraphic graphic = new SimulationGraphic(sim,SimulationGraphic.TABBED_PANE);
        	AccumulatorHistory densityHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
            rhoAccumulator.addDataSink(densityHistory, new StatType[]{rhoAccumulator.MOST_RECENT});
            DisplayPlot rhoPlot = new DisplayPlot();
        	densityHistory.setDataSink(rhoPlot.getDataSet().makeDataSink());
        	rhoPlot.setLabel("density");
        	graphic.add(rhoPlot);
        	AccumulatorHistory energyHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
            energyAccumulator.addDataSink(energyHistory, new StatType[]{energyAccumulator.MOST_RECENT});
            DisplayPlot energyPlot = new DisplayPlot();
        	energyHistory.setDataSink(energyPlot.getDataSet().makeDataSink());
        	energyPlot.setLabel("energy");
        	graphic.add(energyPlot);
        	AccumulatorHistory smerHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
            smerAccumulator.addDataSink(smerHistory, new StatType[]{smerAccumulator.MOST_RECENT});
            DisplayPlot smerPlot = new DisplayPlot();
        	smerHistory.setDataSink(smerPlot.getDataSet().makeDataSink());
        	smerPlot.setLabel("smer fraction");
        	graphic.add(smerPlot);
        	DataSourceCountSteps stepCounter = new DataSourceCountSteps(sim.integrator);
        	AccumulatorHistory energy2History = new AccumulatorHistory(new HistoryCollapsingAverage());
        	energyHistory.setTimeDataSource(stepCounter);
        	energy2History.setTimeDataSource(stepCounter);
            energy2Accumulator.addDataSink(energy2History, new StatType[]{energy2Accumulator.MOST_RECENT});
            //DisplayPlot energy2Plot = new DisplayPlot();
        	energy2History.setDataSink(energyPlot.getDataSet().makeDataSink());
//        	energy2Plot.setLabel("energy2");
//        	graphic.add(energy2Plot);
        	ColorSchemeSmer colorScheme = new ColorSchemeSmer(sim.associationHelper,sim.box,sim.getRandom());
        	graphic.getDisplayBox(sim.box).setColorScheme(colorScheme);
        	graphic.makeAndDisplayFrame();
        	return;
        }
        IAction energyDiffAction = new IAction() {

			public void actionPerformed() {
				IAtomOriented atom253 = (IAtomOriented)sim.box.getLeafList().get(253);
				IAtomOriented atom63 = (IAtomOriented)sim.box.getLeafList().get(63);
				AtomArrayList bondList = (AtomArrayList)sim.associationManagerOriented.getAssociatedAtoms(atom253);
				boolean isBonded1 = bondList.indexOf(atom63)>-1;
				boolean isBonded2 = sim.bvso.isAssociated(atom253, atom63);
				if (isBonded1 != isBonded2){
					System.out.println(sim.integrator.getStepCount()+ " steps");
					System.out.println("Wrong. isBonded1: "+isBonded1+" isBonded2: "+isBonded2);
					System.out.println("position of atom253: "+atom253.getPosition()+" orientation of atom253: "+atom253.getOrientation().getDirection());
					System.out.println("position of atom63: "+atom63.getPosition()+" orientation of atom63: "+atom63.getOrientation().getDirection());
					throw new RuntimeException();
				}
				if (sim.integrator.getStepCount() > 1){
					System.out.println(sim.integrator.getStepCount()+ " steps");
					System.out.println("energy= "+energyMeter.getDataAsScalar()+" true energy= "+energy2Meter.getDataAsScalar());
					System.out.println("position of atom253: "+atom253.getPosition()+" orientation of atom263: "+atom253.getOrientation().getDirection());
					System.out.println("position of atom63: "+atom63.getPosition()+" orientation of atom63: "+atom63.getOrientation().getDirection());
					//MCMoveVolumeAssociated.dodebug =true;
					IntegratorMC.dodebug = true;
					double energyDifference = energyMeter.getDataAsScalar()-energy2Meter.getDataAsScalar();
					if (sim.integrator.getStepCount() > 3807000 || Math.abs(energyDifference)> 1E-7){
						System.exit(1);
					}
				}
				else if (sim.integrator.getStepCount()%1000 == 0){
					double energyDifference = energyMeter.getDataAsScalar()-energy2Meter.getDataAsScalar();
					if (Math.abs(energyDifference)> 1E-7){
						System.out.println(sim.integrator.getStepCount()+ " steps");
						System.out.println("energy= "+energyMeter.getDataAsScalar()+" true energy= "+energy2Meter.getDataAsScalar());
						throw new RuntimeException();
					}
				}


			}
		};
        IntegratorListenerAction energyDiffListener = new IntegratorListenerAction(energyDiffAction,1000);
        //sim.integrator.getEventManager().addListener(energyDiffListener)

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));
        //Meter for measurement of the total molecule number density((number of molecules)/(volume of box)) in a box
        
        System.out.println("numAtom=" +numAtoms);
        double avgDensity = ((DataDouble) ((DataGroup) rhoAccumulator.getData()).getData(rhoAccumulator.AVERAGE.index)).x;//average density
        System.out.println("average density= " +avgDensity);
        double Z = pressure/(avgDensity*sim.integrator.getTemperature());
        double avgPE = ((DataDouble) ((DataGroup) energyAccumulator.getData()).getData(energyAccumulator.AVERAGE.index)).x;
        System.out.println("average energy= "+avgPE);
        avgPE /= numAtoms;
        System.out.println("Z="+Z);
        System.out.println("PE/epsilon="+avgPE);
        double avgDimerFraction = ((DataDouble) ((DataGroup) smerAccumulator.getData()).getData(smerAccumulator.AVERAGE.index)).x;
        System.out.println("average fraction of smer= "+avgDimerFraction);
        double temp = sim.integrator.getTemperature();
        double Cv = ((DataDouble) ((DataGroup) energyAccumulator.getData()).getData(energyAccumulator.STANDARD_DEVIATION.index)).x;
        Cv /= temp;
        Cv *= Cv/numAtoms;
        System.out.println("Cv/k="+Cv);
        
        if (Double.isNaN(Z) || Math.abs(Z+0.25) > 0.15) {
            System.exit(1);
        }
        if (Double.isNaN(avgPE) || Math.abs(avgPE+4.56) > 0.03) {
            System.exit(1);
        }
        if (Double.isNaN(Cv) || Math.abs(Cv-0.61) > 0.45) {  // actual average seems to be 0.51
            System.exit(1);
        }
          
    }
    public static class VirialAssociatingFluidParam extends ParameterBase {
		public int numAtoms = 256;
		public double pressure = 0.001;
		public double density = 0.001;
		public double wellConstant = 16.0;
		public double temperature = 2.0;
		public double truncationRadius =6.0;
		public int maxChainLength = Integer.MAX_VALUE;//default: Integer.MAX_VALUE
		public boolean useUB = true;
		public long numSteps = 200000;
	}

}
