/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedvirial;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.data.meter.MeterPressure;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

public class MappedVirialLJ extends Simulation {
    
    public SpeciesSpheresMono species;
    public IBox box;
    public IntegratorMC integrator;
    public MCMoveAtom move;
    public ActivityIntegrate activityIntegrate;
    public P2SoftSphericalTruncated p2Truncated;
    
    public MappedVirialLJ(ISpace _space, int numAtoms, double temperature, double density, double rc) {
        super(_space);
        PotentialMasterCell potentialMaster = new PotentialMasterCell(this, rc, space);
        potentialMaster.lrcMaster().setEnabled(false);
        
        //controller and integrator
	    integrator = new IntegratorMC(potentialMaster, random, temperature);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        move = new MCMoveAtom(random, potentialMaster, space);
        integrator.getMoveManager().addMCMove(move);

	    //species and potentials
	    species = new SpeciesSpheresMono(this, space);
        addSpecies(species);
        
	    P2LennardJones potential = new P2LennardJones(space);
        p2Truncated = new P2SoftSphericalTruncated(space, potential, rc);
	    potentialMaster.addPotential(p2Truncated, new IAtomType[]{species.getLeafType(), species.getLeafType()});
	    
        //construct box
	    box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
        
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();
        
        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        integrator.setBox(box);
        potentialMaster.setCellRange(2);

        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());
    }
    
    public static void main(String[] args) {
        
        LJMDParams params = new LJMDParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.density = 0.025;
        }

        int numAtoms = params.numAtoms;
        double temperature = params.temperture;
        double density = params.density;
        long numSteps = params.numSteps;
        double rc = params.rc;
        
        ISpace space = Space.getInstance(3);

        MappedVirialLJ sim = new MappedVirialLJ(space, numAtoms, temperature, density, rc);
        
        sim.activityIntegrate.setMaxSteps(numSteps/10);
        sim.getController().actionPerformed();
        sim.getController().reset();
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.integrator.getMoveManager().setEquilibrating(false);
        
        int nBins = 10000;
        int numBlocks = 100;
        long numSamples = numSteps/numAtoms;
        long samplesPerBlock = numSamples/numBlocks;
        if (samplesPerBlock == 0) samplesPerBlock = 1;

        double[] cutoff = new double[]{0.95, 0.97, 1.0, 1.1, 1.2, 1.4, 1.492, 1.6, 2.0, 2.5, 3.0, 4.0};
        final MeterMappedVirial meterMappedVirial = new MeterMappedVirial(space, sim.integrator.getPotentialMaster(), 
                sim.p2Truncated, sim.box, nBins, cutoff);
        meterMappedVirial.setTemperature(sim.integrator.getTemperature());
        final AccumulatorAverageFixed accMappedVirial = new AccumulatorAverageFixed(samplesPerBlock);
        DataPumpListener pumpMappedVirial = new DataPumpListener(meterMappedVirial, accMappedVirial, numAtoms);
        sim.integrator.getEventManager().addListener(pumpMappedVirial);

        final MeterPressure meterP = new MeterPressure(space);
        meterP.setIntegrator(sim.integrator);
        final AccumulatorAverageFixed accP = new AccumulatorAverageFixed(samplesPerBlock);
        DataPumpListener pumpP = new DataPumpListener(meterP, accP, numAtoms);
        sim.integrator.getEventManager().addListener(pumpP);

        
        sim.getController().actionPerformed();
        
        IData mappedAvg = accMappedVirial.getData(accMappedVirial.AVERAGE);
        IData mappedErr = accMappedVirial.getData(accMappedVirial.ERROR);
        IData mappedCor = accMappedVirial.getData(accMappedVirial.BLOCK_CORRELATION);
        for (int i=0; i<cutoff.length; i++) {
            double avg = mappedAvg.getValue(i);
            double err = mappedErr.getValue(i);
            double cor = mappedCor.getValue(i);
            System.out.print(String.format("xc: %6.3f   avg: %12.5e   err: %11.4e   cor: %4.2f\n", cutoff[i], avg, err, cor));
        }

        double pAvg = accP.getData(accP.AVERAGE).getValue(0);
        double pErr = accP.getData(accP.ERROR).getValue(0);
        double pCor = accP.getData(accP.BLOCK_CORRELATION).getValue(0);
        System.out.print(String.format("Pressure     avg: %12.5e   err: %11.4e   cor: %4.2f\n", pAvg, pErr, pCor));

    }
    
    public static class LJMDParams extends ParameterBase {
        public int numAtoms = 100;
        public double temperture = 1.0;
        public double density = 0.01;
        public long numSteps = 1000000;
        public double rc = 4;
    }
}
