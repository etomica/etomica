/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.controller.Activity;
import etomica.action.controller.Controller;
import etomica.atom.AtomTypeOriented;
import etomica.data.*;
import etomica.data.types.DataGroup;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.potential.IPotential2;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.NeighborManagerIntra;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.BoxCluster;
import etomica.virial.ConfigurationCluster;
import etomica.virial.MeterVirial;
import etomica.virial.cluster.ClusterAbstract;
import etomica.virial.cluster.ClusterWeight;
import etomica.virial.mcmove.*;

/**
 * Generic simulation using Mayer sampling to evaluate cluster integrals
 */
public class SimulationVirial extends Simulation {


    public IDataSource meter;
    public AccumulatorRatioAverageCovariance accumulator;
    public IDataSink dataSink;
    public DataPumpListener accumulatorPump;
    public ISpecies[] species;

    public IntegratorMC integrator;
    public BoxCluster box;
    public ClusterAbstract[] allValueClusters;
    public ClusterWeight sampleCluster;
    public MCMoveBoxStep mcMoveTranslate;
    public MCMoveBoxStep mcMoveRotate;
    public MCMoveBoxStep mcMoveWiggle;
    public double temperature;
    public boolean doWiggle;
    public int[] newSeeds;
    public int[] numMolecules;
    protected double boxLength;
    protected PotentialMasterBonding.FullBondingInfo bondingInfo;
    protected IPotential2[][] pairPotentials;
    protected boolean initialized;

    public SimulationVirial(Space space, ISpecies[] species, int[] nMolecules, double temperature, ClusterWeight aSampleCluster, ClusterAbstract refCluster, ClusterAbstract[] targetClusters) {
        super(space);
        this.species = species;
        this.numMolecules = nMolecules;
        this.temperature = temperature;
        this.sampleCluster = aSampleCluster;
        allValueClusters = new ClusterAbstract[targetClusters.length+1];
        allValueClusters[0] = refCluster;
        System.arraycopy(targetClusters,0,allValueClusters,1,targetClusters.length);
    }

    public void setDoWiggle(boolean newDoWiggle) {
        this.doWiggle = newDoWiggle;
    }

    public void setSeeds(int[] newSeeds) {
        this.seeds = newSeeds;
    }

    public void setBoxLength(double length) {
        boxLength = length;
    }

    public void setIntraPairPotentials(IPotential2[][] pairPotentials) {
        if (initialized) throw new RuntimeException("too late");
        this.pairPotentials = pairPotentials;
    }

    public void setBondingInfo(PotentialMasterBonding.FullBondingInfo bondingInfo) {
        this.bondingInfo = bondingInfo;
    }

    public void init() {
        if (initialized) throw new RuntimeException("you can only call me once");
        // we aren't actually initialized yet, but we will be unless we crash.
        // if we crash, we shouldn't get called again!
        initialized = true;

        if (seeds != null) {
            setRandom(new RandomMersenneTwister(seeds));
        }
        for (int i = 0; i < species.length; i++) {
            addSpecies(species[i]);
        }

        box = new BoxCluster(sampleCluster, space, boxLength);
        addBox(box);

        for (int i = 0; i < species.length; i++) {
            box.setNMolecules(species[i], numMolecules[i]);
        }

        PotentialCompute pc;
        if (pairPotentials != null) {
            PotentialMasterBonding pmBonding = new PotentialMasterBonding(getSpeciesManager(), box, bondingInfo);
            PotentialComputePair pcPair = new PotentialComputePair(getSpeciesManager(), box, new NeighborManagerIntra(box, bondingInfo), pairPotentials);
            pc = new PotentialComputeAggregate(pmBonding, pcPair);
        }
        else if (bondingInfo != null){
            pc = new PotentialMasterBonding(getSpeciesManager(), box, bondingInfo);
        }
        else {
            pc = new PotentialComputeAggregate();
        }

        // temperature isn't going to mean anything here, but pass it anyway
        integrator = new IntegratorMC(pc, random, temperature, box);
        integrator.getMoveManager().setEquilibrating(false);
        integrator.setEventInterval(1);
        getController().addActivity(new ActivityIntegrate(integrator));

        boolean doRotate = false;
        boolean multiAtomic = false;
        for (int i=0; i<species.length; i++) {
            ISpecies sp = species[i];
            if (sp.getLeafAtomCount() == 1 && sp.getLeafType() instanceof AtomTypeOriented) {
                doRotate = true;
            }
            if (sp.getLeafAtomCount() > 1) {
                multiAtomic = true;
            }
        }

        if (!multiAtomic) {
            mcMoveTranslate = new MCMoveClusterAtomMulti(random, box);
            if (doRotate) {
                mcMoveRotate = new MCMoveClusterAtomRotateMulti(random, box);
                integrator.getMoveManager().addMCMove(mcMoveRotate);
            }
        } else {
            mcMoveTranslate = new MCMoveClusterMoleculeMulti(random, box);
            mcMoveRotate = new MCMoveClusterRotateMoleculeMulti(getRandom(), box);
            mcMoveRotate.setStepSize(Math.PI);
            if (doWiggle) {
                mcMoveWiggle = new MCMoveClusterWiggleMulti(random, pc, box);
                integrator.getMoveManager().addMCMove(mcMoveWiggle);
            }
            integrator.getMoveManager().addMCMove(mcMoveRotate);
        }
        integrator.getMoveManager().addMCMove(mcMoveTranslate);

        ConfigurationCluster configuration = new ConfigurationCluster(space);
        configuration.initializeCoordinates(box);

        setMeter(new MeterVirial(allValueClusters));
        ((MeterVirial) meter).setBox(box);
        setAccumulator(new AccumulatorRatioAverageCovariance());
    }

    public void setMeter(IDataSource newMeter) {
        meter = newMeter;
        if (accumulator != null) {
            if (accumulatorPump != null) {
                integrator.getEventManager().removeListener(accumulatorPump);
                accumulatorPump = null;
            }
            if (meter != null) {
                setAccumulator(accumulator);
            }
        }
	}

	public void setAccumulator(IDataSink newAccumulator) {
	    dataSink = newAccumulator;
	    if (newAccumulator instanceof AccumulatorRatioAverageCovariance) {
	        accumulator = (AccumulatorRatioAverageCovariance)newAccumulator;
	    }
	    else {
	        accumulator = null;
	    }
		if (accumulatorPump == null) {
			accumulatorPump = new DataPumpListener(meter,newAccumulator);
            integrator.getEventManager().addListener(accumulatorPump);
		}
		else {
			accumulatorPump.setDataSink(newAccumulator);
		}
	}

	public void setAccumulatorBlockSize(long newBlockSize) {
	    accumulator.setBlockSize(newBlockSize);
	}

	public void equilibrate(long initSteps) {
        // run a short simulation to get reasonable MC Move step sizes and
        // (if needed) narrow in on a reference preference
        getController().runActivityBlocking(makeEquilibrationActivity(initSteps));
    }

    private Activity makeEquilibrationActivity(long initSteps) {
        ActivityIntegrate ai = new ActivityIntegrate(integrator, initSteps);
        return new Activity() {
            @Override
            public void runActivity(Controller.ControllerHandle handle) {
                integrator.getMoveManager().setEquilibrating(true);
                ai.runActivity(handle);
                integrator.getMoveManager().setEquilibrating(false);
                if (accumulator != null) {
                    accumulator.reset();
                }
            }
        };
    }

    public Controller.ActivityHandle<?> addEquilibration(long initSteps) {
        return getController().addActivity(makeEquilibrationActivity(initSteps));
    }


    public void printResults(double refIntegral) {
        printResults(refIntegral, null);
    }

    public void printResults(double refIntegral, String[] extraNames) {

        DataGroup allYourBase = (DataGroup)accumulator.getData();
        IData averageData = allYourBase.getData(accumulator.AVERAGE.index);
        IData stdevData = allYourBase.getData(accumulator.STANDARD_DEVIATION.index);
        IData errorData = allYourBase.getData(accumulator.ERROR.index);
        IData correlationData = allYourBase.getData(accumulator.BLOCK_CORRELATION.index);
        IData ratioData = allYourBase.getData(accumulator.RATIO.index);
        IData ratioErrorData = allYourBase.getData(accumulator.RATIO_ERROR.index);
        IData covarianceData = allYourBase.getData(accumulator.BLOCK_COVARIANCE.index);

        System.out.println();
        System.out.print(String.format("reference average: %20.15e stdev: %9.4e error: %9.4e cor: %6.4f\n",
                averageData.getValue(0), stdevData.getValue(0), errorData.getValue(0), correlationData.getValue(0)));

        System.out.print(String.format("target average: %20.15e stdev: %9.4e error: %9.4e cor: %6.4f\n",
                averageData.getValue(1), stdevData.getValue(1), errorData.getValue(1), correlationData.getValue(1)));

        int nData = averageData.getLength();
        int nCovData = covarianceData.getLength();
        if (nData*nData != nCovData) {
            // we need to know these to grab the right elements of covarianceData
            throw new RuntimeException("unexpected number of data values");
        }
        double correlationCoef = covarianceData.getValue(1)/Math.sqrt(covarianceData.getValue(0)*covarianceData.getValue(nData+1));
        correlationCoef = (Double.isNaN(correlationCoef) || Double.isInfinite(correlationCoef)) ? 0 : correlationCoef;

        System.out.println();
        if (averageData.getValue(0) != 0) {
            System.out.print(String.format("ratio average: %20.15e  error: %9.4e  cor: %6.4f\n", ratioData.getValue(1), ratioErrorData.getValue(1), correlationCoef));
            System.out.print(String.format("abs average: %20.15e  error: %9.4e\n", ratioData.getValue(1) * refIntegral, ratioErrorData.getValue(1) * Math.abs(refIntegral)));
        }

        int n = averageData.getLength();
        double[] var = new double[n];
        for (int i=0; i<n; i++) {
            var[i] = covarianceData.getValue((i) * n + (i));
        }

        for (int i=2; i<averageData.getLength(); i++) {
            String name = String.format("%d", i);
            if (extraNames != null) {
                name = extraNames[i-2];
            }
            System.out.print(String.format("target %s average: %20.15e stdev: %9.4e error: %9.4e cor: %6.4f tcor: ",
                    name, averageData.getValue(i), stdevData.getValue(i), errorData.getValue(i), correlationData.getValue(i)));
            for (int j=1; j<i; j++) {
                double c = var[i]*var[j] == 0 ? 0 : covarianceData.getValue((i)*n+(j))/Math.sqrt(var[i]*var[j]);
                System.out.print(String.format(" %20.18f", c));
            }
            System.out.println();
            if (averageData.getValue(0) != 0) {
                correlationCoef = covarianceData.getValue(i) / Math.sqrt(covarianceData.getValue(0) * covarianceData.getValue((i - 1) * nData + i));
                correlationCoef = (Double.isNaN(correlationCoef) || Double.isInfinite(correlationCoef)) ? 0 : correlationCoef;
                System.out.print(String.format("ratio %s average: %20.15e  error: %9.4e  cor: %6.4f\n", name, ratioData.getValue(i), ratioErrorData.getValue(i), correlationCoef));
                System.out.print(String.format("full %s average: %20.15e  error: %9.4e\n", name, ratioData.getValue(i) * refIntegral, ratioErrorData.getValue(i) * Math.abs(refIntegral)));
                System.out.println();
            }
        }
    }

    @Override
    public Integrator getIntegrator() {
        return integrator;
    }
}

