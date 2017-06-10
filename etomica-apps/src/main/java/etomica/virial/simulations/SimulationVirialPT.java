/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.data.*;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorPT;
import etomica.integrator.mcmove.MCMove;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresMono;
import etomica.species.SpeciesSpheresRotating;
import etomica.virial.*;
import etomica.virial.paralleltempering.MCMoveSwapCluster;

/**
 * Generic simulation using Mayer sampling to evaluate cluster integrals
 */
public class SimulationVirialPT extends Simulation {

	/**
	 * Constructor for simulation to determine the ratio bewteen reference and target Clusters
	 */
	public SimulationVirialPT(Space space, SpeciesFactory speciesFactory, 
			double[] temperature, ClusterWeight.Factory sampleClusterFactory, 
			ClusterAbstract refCluster, ClusterAbstract[] targetClusters) {
		super(space);
        PotentialMaster potentialMaster = new PotentialMaster();
		int nMolecules = refCluster.pointCount();
		species = speciesFactory.makeSpecies(space);//SpheresMono(this,AtomLinker.FACTORY);
        addSpecies(species);

        box = new BoxCluster[temperature.length];
        integrator = new IntegratorMC[temperature.length];
        meter = new MeterVirial[temperature.length];
        accumulator = new DataAccumulator[temperature.length];
        accumulatorPump = new DataPumpListener[temperature.length];
        mcMoveMulti = new MCMoveBoxStep[temperature.length];
        mcMoveRotate = new MCMoveBoxStep[temperature.length];
        meterAccept = new IEtomicaDataSource[temperature.length-1];
        meterAcceptP = new IEtomicaDataSource[temperature.length-1];
        
        // Parallel tempering would sorta work without separate instances of the clusters
        // but the value caching based on coordinate pair set ID would be confused because
        // the coordinate pair sets from different boxs could have the same ID and different
        // values
        sampleCluster = new ClusterWeight[temperature.length];
        allValueClusters = new ClusterAbstract[temperature.length][targetClusters.length+1];
        integratorPT = new IntegratorPT(getRandom(),MCMoveSwapCluster.FACTORY, space);
//        integratorPT.setSwapInterval(2);
        integratorPT.setEventInterval(1);
        ai = new ActivityIntegrate(integratorPT);
        getController().addAction(ai);
        
        for (int iTemp=0; iTemp<temperature.length; iTemp++) {
            allValueClusters[iTemp][0] = refCluster.makeCopy();
            allValueClusters[iTemp][0].setTemperature(temperature[iTemp]);
            for (int i=0; i<targetClusters.length; i++) {
                allValueClusters[iTemp][i+1] = targetClusters[i].makeCopy();
                allValueClusters[iTemp][i+1].setTemperature(temperature[iTemp]);
            }

            sampleCluster[iTemp] = sampleClusterFactory.makeWeightCluster(allValueClusters[iTemp]);
            sampleCluster[iTemp].setTemperature(temperature[iTemp]);
            box[iTemp] = new BoxCluster(sampleCluster[iTemp], space);
            box[iTemp].setNMolecules(species, nMolecules);
            
            integrator[iTemp] = new IntegratorMC(this, potentialMaster);
            integrator[iTemp].setTemperature(temperature[iTemp]);
            integrator[iTemp].setBox(box[iTemp]);
            integrator[iTemp].getMoveManager().setEquilibrating(false);
            integratorPT.addIntegrator(integrator[iTemp]);
            
            MCMoveManager moveManager = integrator[iTemp].getMoveManager();
            
            if (species instanceof SpeciesSpheresMono || species instanceof SpeciesSpheresRotating) {
                mcMoveMulti[iTemp] = new MCMoveClusterAtomMulti(random, space);
                moveManager.addMCMove(mcMoveMulti[iTemp]);
            }
            else {
                mcMoveMulti[iTemp] = new MCMoveClusterMoleculeMulti(this, space);
                moveManager.addMCMove(mcMoveMulti[iTemp]);
                mcMoveRotate[iTemp] = new MCMoveClusterRotateMoleculeMulti(getRandom(), space);
                moveManager.addMCMove(mcMoveRotate[iTemp]);
            }
            
            ConfigurationCluster configuration = new ConfigurationCluster(space);
            configuration.initializeCoordinates(box[iTemp]);
            
            setMeter(iTemp,new MeterVirial(allValueClusters[iTemp]));
//            meter[iTemp].getDataInfo().setLabel("Target/Refernce Ratio "+iTemp);
            setAccumulator(iTemp,new AccumulatorRatioAverage());

            if(iTemp>0) {
                MCMove swapMove = integratorPT.getMoveManager().getMCMoves().get(iTemp-1);
                meterAccept[iTemp-1] = new DataSourceAcceptanceRatio(swapMove);
                meterAcceptP[iTemp-1] = new DataSourceAcceptanceProbability(swapMove);
            }
        }
	}
	
	public MeterVirial[] meter;
    public IEtomicaDataSource[] meterAccept;
    public IEtomicaDataSource[] meterAcceptP;
	public DataAccumulator[] accumulator;
	public DataPumpListener[] accumulatorPump;
	public ISpecies species;
	public ActivityIntegrate ai;
	public IntegratorMC[] integrator;
	public BoxCluster[] box;
    public ClusterAbstract[][] allValueClusters;
    public ClusterWeight[] sampleCluster;
    public MCMoveBoxStep[] mcMoveRotate;
    public MCMoveBoxStep[] mcMoveMulti;
    public IntegratorPT integratorPT;

	public void setMeter(int i, MeterVirial newMeter) {
		if (accumulator[i] != null) { 
			if (accumulatorPump[i] != null) {
                integrator[i].getEventManager().removeListener(accumulatorPump[i]);
				accumulatorPump[i] = null;
			}
			accumulator[i] = null;
		}
		meter[i] = newMeter;
		if (meter[i] != null) {
			meter[i].setBox(box[i]);
		}
	}

	public void setAccumulator(int i, DataAccumulator newAccumulator) {
		accumulator[i] = newAccumulator;
		if (accumulatorPump[i] == null) {
			accumulatorPump[i] = new DataPumpListener(meter[i],accumulator[i]);
		}
		else {
			accumulatorPump[i].setDataSink(accumulator[i]);
		}

		integrator[i].getEventManager().addListener(accumulatorPump[i]);
	}
}

