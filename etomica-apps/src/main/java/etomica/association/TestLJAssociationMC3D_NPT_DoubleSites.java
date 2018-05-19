/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.action.BoxInflate;
import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomType;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
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
import etomica.potential.P2HardAssociationConeDoubleSites;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.Degree;
import etomica.util.ParameterBase;

/**
 * Simple Lennard-Jones + S-W Association Monte Carlo NPT simulation in 3D.
 * average density = N*<1/V>
 * Initial configurations at http://rheneas.eng.buffalo.edu/etomica/tests/
 * @author Hye Min Kim
 */
public class TestLJAssociationMC3D_NPT_DoubleSites extends Simulation {
    
    public IntegratorMC integrator;
    public MCMoveAtomMonomer mcMoveAtomMonomer;
    public MCMoveAtomSmer mcMoveAtomSmer;
    public MCMoveRotateAssociated mcMoveRotate;
    public SpeciesSpheresRotating species;
    public Box box;
    public P2HardAssociationConeDoubleSites potential;
    public MCMoveSmer mcMoveSmer;
    public MCMoveSmerRotate mcMoveSmerRotate;
    public MCMoveVolumeAssociated mcMoveVolume;
    public ActivityIntegrate actionIntegrator;
    public MCMoveBiasUB mcMoveBiasUB;
    public AssociationManager associationManagerOriented;
    public BiasVolumeSphereOrientedDoubleSites bvso;
    public AssociationHelperDouble associationHelper;
    double epsilon = 1.0;
        
    
    public TestLJAssociationMC3D_NPT_DoubleSites(int numAtoms, double pressure, double density, double wellConstant, double temperature,double truncationRadius,int maxChainLength, boolean useUB, long numSteps) {
        super(Space3D.getInstance());

        species = new SpeciesSpheresRotating(this, space);//Species in which molecules are made of a single atom of type OrientedSphere
        addSpecies(species);

        PotentialMasterCell potentialMaster = new PotentialMasterCell(this, space);

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
        integrator = new IntegratorMC(this, potentialMaster, box);
        integrator.setTemperature(temperature);
        mcMoveAtomMonomer = new MCMoveAtomMonomer(this, potentialMaster, space);//Standard Monte Carlo atom-displacement trial move
        mcMoveAtomMonomer.setMaxLength(maxChainLength);
        mcMoveAtomSmer = new MCMoveAtomSmer(this, potentialMaster, space);
        mcMoveAtomSmer.setMaxLength(maxChainLength);
        mcMoveRotate = new MCMoveRotateAssociated(potentialMaster, random, space);//Performs a rotation of an atom (not a molecule) that has an orientation coordinate
        mcMoveRotate.setMaxLength(maxChainLength);
        bvso = new BiasVolumeSphereOrientedDoubleSites(space, random);
        System.out.println("biasVolume= 2*(partialVolume/totalVolume)^2");
        bvso.setTheta(Degree.UNIT.toSim(27.0));
        bvso.setBiasSphereInnerRadius(0.0);
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
        actionIntegrator = new ActivityIntegrate(integrator);
        //actionIntegrate.setSleepPeriod(1);
        actionIntegrator.setMaxSteps(numSteps);
        getController().addAction(actionIntegrator);
        BoxInflate inflater = new BoxInflate(box, space);//Performs actions that cause volume of system to expand or contract
        inflater.setTargetDensity(density);
        inflater.actionPerformed();

        System.out.println("truncation distance of potential = " + truncationRadius);
        if (truncationRadius > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large.  Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        potential = new P2HardAssociationConeDoubleSites(space, sigma, epsilon, truncationRadius, wellConstant);
        potentialMaster.setCellRange(3);
        potentialMaster.setRange(potential.getRange());
        mcMoveSmer = new MCMoveSmer(this, potentialMaster, space, potential);
        mcMoveSmerRotate = new MCMoveSmerRotate(this, potentialMaster, space, potential);
        mcMoveVolume = new MCMoveVolumeAssociated(this, potentialMaster, space);
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
        final TestLJAssociationMC3D_NPT_DoubleSites sim = new TestLJAssociationMC3D_NPT_DoubleSites(numAtoms, pressure, density, wellConstant, temperature, truncationRadius,maxChainLength, useUB, numSteps);
        sim.actionIntegrator.setMaxSteps(numSteps/5);//equilibrium period
        IAction energyDiffActionEq = new IAction() {
    		
			public void actionPerformed() {
				//if (sim.integrator.getStepCount()%1000 == 0){
					IAtomList leafList = sim.box.getLeafList();
					for (int i = 0; i <leafList.size(); i+=1){
						IAtomOriented atomi = (IAtomOriented)sim.box.getLeafList().get(i);
						for (int j =0; j<i;j+=1){
							IAtomOriented atomj = (IAtomOriented)sim.box.getLeafList().get(j);
							AtomArrayList bondList = (AtomArrayList)sim.associationManagerOriented.getAssociatedAtoms(atomi);
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
				
		
//			}
		};
        IntegratorListenerAction energyDiffListenerEq = new IntegratorListenerAction(energyDiffActionEq,100000);
        sim.integrator.getEventManager().addListener(energyDiffListenerEq);
        System.out.println("equilibrium period = " +numSteps/5);
        sim.getController().actionPerformed();
        System.out.println("equilibrium finished");
        sim.getController().reset();
        
        sim.actionIntegrator.setMaxSteps(numSteps);
        MeterDensity rhoMeter = new MeterDensity(sim.space);//Meter for measurement of the total molecule number density((number of molecules)/(volume of box)) in a box 
        rhoMeter.setBox(sim.box);
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
            rhoAccumulator.addDataSink(densityHistory, new StatType[]{AccumulatorAverage.MOST_RECENT});
            DisplayPlot rhoPlot = new DisplayPlot();
        	densityHistory.setDataSink(rhoPlot.getDataSet().makeDataSink());
        	rhoPlot.setLabel("density");
        	graphic.add(rhoPlot);
        	AccumulatorHistory energyHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
            energyAccumulator.addDataSink(energyHistory, new StatType[]{AccumulatorAverage.MOST_RECENT});
            DisplayPlot energyPlot = new DisplayPlot();
        	energyHistory.setDataSink(energyPlot.getDataSet().makeDataSink());
        	energyPlot.setLabel("energy");
        	graphic.add(energyPlot);
        	AccumulatorHistory smerHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
            smerAccumulator.addDataSink(smerHistory, new StatType[]{AccumulatorAverage.MOST_RECENT});
            DisplayPlot smerPlot = new DisplayPlot();
        	smerHistory.setDataSink(smerPlot.getDataSet().makeDataSink());
        	smerPlot.setLabel("smer fraction");
        	graphic.add(smerPlot);
        	DataSourceCountSteps stepCounter = new DataSourceCountSteps(sim.integrator);
        	AccumulatorHistory energy2History = new AccumulatorHistory(new HistoryCollapsingAverage());
        	energyHistory.setTimeDataSource(stepCounter);
        	energy2History.setTimeDataSource(stepCounter);
            energy2Accumulator.addDataSink(energy2History, new StatType[]{AccumulatorAverage.MOST_RECENT});
            //DisplayPlot energy2Plot = new DisplayPlot();
        	energy2History.setDataSink(energyPlot.getDataSet().makeDataSink());
//        	energy2Plot.setLabel("energy2");
//        	graphic.add(energy2Plot);
        	ColorSchemeSmer colorScheme = new ColorSchemeSmer(sim.associationHelper,sim.box,sim.getRandom());
        	graphic.getDisplayBox(sim.box).setColorScheme(colorScheme);
        	graphic.makeAndDisplayFrame();
        	sim.actionIntegrator.setMaxSteps(Long.MAX_VALUE);
        	return;
        }
        IAction energyDiffAction = new IAction() {
		
			public void actionPerformed() {
				IAtomOriented atom207 = (IAtomOriented)sim.box.getLeafList().get(207);
				IAtomOriented atom58 = (IAtomOriented)sim.box.getLeafList().get(58);
				AtomArrayList bondList = (AtomArrayList)sim.associationManagerOriented.getAssociatedAtoms(atom207);
				boolean isBonded1 = bondList.indexOf(atom58)>-1;
				boolean isBonded2 = sim.bvso.isAssociated(atom207, atom58);
				if (isBonded1 != isBonded2){
					System.out.println(sim.integrator.getStepCount()+ " steps");
					System.out.println("Wrong. isBonded1: "+isBonded1+" isBonded2: "+isBonded2);
					System.out.println("position of atom207: "+atom207.getPosition()+" orientation of atom207: "+atom207.getOrientation().getDirection());
					System.out.println("position of atom58: "+atom58.getPosition()+" orientation of atom58: "+atom58.getOrientation().getDirection());
					throw new RuntimeException();
				}
				if (sim.integrator.getStepCount() > 3805999){
					System.out.println(sim.integrator.getStepCount()+ " steps");
					System.out.println("energy= "+energyMeter.getDataAsScalar()+" true energy= "+energy2Meter.getDataAsScalar());
					System.out.println("position of atom207: "+atom207.getPosition()+" orientation of atom207: "+atom207.getOrientation().getDirection());
					System.out.println("position of atom58: "+atom58.getPosition()+" orientation of atom58: "+atom58.getOrientation().getDirection());
					//MCMoveVolumeAssociated.dodebug =true;
					IntegratorMC.dodebug = true;
					double energyDifference = energyMeter.getDataAsScalar()-energy2Meter.getDataAsScalar();
					if (sim.integrator.getStepCount() > 3807000 || Math.abs(energyDifference)> 1E-7){
						System.exit(1);
					}
				} else if (sim.integrator.getStepCount()%1000 == 0){
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
        //sim.integrator.getEventManager().addListener(energyDiffListener);
        
        sim.getController().actionPerformed();
        
        System.out.println("numAtom=" +numAtoms);
        double avgDensity = ((DataDouble) ((DataGroup) rhoAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index)).x;//average density
        System.out.println("average density= " +avgDensity);
        double Z = pressure/(avgDensity*sim.integrator.getTemperature());
        double avgPE = ((DataDouble) ((DataGroup) energyAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index)).x;
        System.out.println("average energy= "+avgPE);
        avgPE /= numAtoms;
        System.out.println("Z="+Z);
        System.out.println("PE/epsilon="+avgPE);
        double avgDimerFraction = ((DataDouble) ((DataGroup) smerAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index)).x;
        System.out.println("average fraction of smer= "+avgDimerFraction);
        double temp = sim.integrator.getTemperature();
        double Cv = ((DataDouble) ((DataGroup) energyAccumulator.getData()).getData(AccumulatorAverage.STANDARD_DEVIATION.index)).x;
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
		public int numAtoms = 1024;
		public double pressure = 0.17;
		public double density = 0.2;
		public double wellConstant = 16.0;
		public double temperature = 2.0;
		public double truncationRadius =6.0;
		public int maxChainLength = Integer.MAX_VALUE;//default: Integer.MAX_VALUE
		public boolean useUB = true;
		public long numSteps = 200000000;
	}

}
