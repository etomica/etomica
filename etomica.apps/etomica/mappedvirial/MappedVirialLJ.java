/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedvirial;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.data.meter.MeterPressure;
import etomica.data.meter.MeterRDF;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.listener.IntegratorListenerAction;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.Histogram;
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
        
        potentialMaster.getNbrCellManager(box).assignCellAll();
        
        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());
    }
    
    public static void main(String[] args) throws IOException {
        
        LJMDParams params = new LJMDParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.temperature = 2;
            params.density = 0.1;
            params.numSteps = 100000;
            params.rc = 2.5;
            params.numAtoms = 5000;
            params.functionsFile = "sim";
        }

        int numAtoms = params.numAtoms;
        double temperature = params.temperature;
        double density = params.density;
        long numSteps = params.numSteps;
        double rc = params.rc;
        String functionsFile = params.functionsFile;

        System.out.println("Virial mapped average");
        System.out.println(numAtoms+" atoms, "+numSteps+" steps");
        System.out.println("density: "+density);
        System.out.println("temperature: "+temperature);
        System.out.println("cutoff: "+rc);
        System.out.println("x0: "+0.95*rc);

        ISpace space = Space.getInstance(3);

        MappedVirialLJ sim = new MappedVirialLJ(space, numAtoms, temperature, density, rc);
        
        long t1 = System.currentTimeMillis();
        
        sim.activityIntegrate.setMaxSteps(numSteps/10);
        sim.getController().actionPerformed();
        sim.getController().reset();
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.integrator.getMoveManager().setEquilibrating(false);
        
        int nBins = 1000000;
        int numBlocks = 100;
        long numSamples = numSteps/numAtoms;
        long samplesPerBlock = numSamples/numBlocks;
        if (samplesPerBlock == 0) samplesPerBlock = 1;

        double[] cutoff0 = new double[]{0.9, 0.95, 0.97, 1.0, 1.1, 1.2, 1.4, 1.492, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};
        double[] cutoff = null;
        for (int i=0; i<cutoff0.length; i++) {
            if (cutoff0[i]*1.15 > params.rc) {
                cutoff = Arrays.copyOfRange(cutoff0, 0, i);
                break;
            }
        }
        if (cutoff==null) cutoff = cutoff0;
        final MeterMappedVirial meterMappedVirial = new MeterMappedVirial(space, sim.integrator.getPotentialMaster(), 
                rc, sim.box, nBins, cutoff);
        meterMappedVirial.getPotentialCalculation().setTemperature(sim.integrator.getTemperature(), sim.p2Truncated);
        final AccumulatorAverageFixed accMappedVirial = new AccumulatorAverageFixed(samplesPerBlock);
        DataPumpListener pumpMappedVirial = new DataPumpListener(meterMappedVirial, accMappedVirial, numAtoms);
        sim.integrator.getEventManager().addListener(pumpMappedVirial);

        final MeterPressure meterP = new MeterPressure(space);
        meterP.setIntegrator(sim.integrator);
        final AccumulatorAverageFixed accP = new AccumulatorAverageFixed(samplesPerBlock);
        DataPumpListener pumpP = new DataPumpListener(meterP, accP, numAtoms);
        sim.integrator.getEventManager().addListener(pumpP);

        MeterMeanForce meterF = null;
        MeterRDF meterRDF = null;
        if (functionsFile != null) {
            int nbins = (int)Math.round(rc/0.01);
            meterF = new MeterMeanForce(space, sim.integrator.getPotentialMaster(), sim.p2Truncated, sim.box, nbins);
//            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterF, numAtoms));
    
            meterRDF = new MeterRDF(space); //, sim.integrator.getPotentialMaster(), sim.box);
            meterRDF.setBox(sim.box);
            meterRDF.getXDataSource().setNValues(nbins);
            meterRDF.getXDataSource().setXMax(rc);
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterRDF, numAtoms));
        }

        sim.getController().actionPerformed();
        
        IData mappedAvg = accMappedVirial.getData(accMappedVirial.AVERAGE);
        IData mappedErr = accMappedVirial.getData(accMappedVirial.ERROR);
        IData mappedCor = accMappedVirial.getData(accMappedVirial.BLOCK_CORRELATION);
        double[] epsilon = meterMappedVirial.getPotentialCalculation().getEpsilon();
        for (int i=0; i<cutoff.length; i++) {
            double avg = mappedAvg.getValue(i);
            double err = mappedErr.getValue(i);
            double cor = mappedCor.getValue(i);
            double q = meterMappedVirial.getPotentialCalculation().getQ(i);
            System.out.print(String.format("xc: %6.3f   eps: %6.3f   q: %18.10e   avg: %13.6e   err: %11.4e   cor: % 4.2f\n", cutoff[i], epsilon[i], q, avg, err, cor));
        }

        
        double pAvg = accP.getData(accP.AVERAGE).getValue(0);
        double pErr = accP.getData(accP.ERROR).getValue(0);
        double pCor = accP.getData(accP.BLOCK_CORRELATION).getValue(0);
        System.out.print(String.format("Pressure     avg: %13.6e   err: %11.4e   cor: % 4.2f\n", pAvg, pErr, pCor));

        if (functionsFile != null) {
            FileWriter fw = new FileWriter(functionsFile+"_gr.dat");
            IData rdata = meterRDF.getIndependentData(0);
            IData gdata = meterRDF.getData();
            for (int i=0; i<rdata.getLength(); i++) {
                double r = rdata.getValue(i);
                double g = gdata.getValue(i);
                double e = Math.exp(-sim.p2Truncated.u(r*r)/temperature);
                fw.write(String.format("%5.3f %22.15e %22.15e\n", r, g, e));
            }
            fw.close();
            
            fw = new FileWriter(functionsFile+"_mf.dat");
            rdata = meterF.getIndependentData(0);
            IData fdata = meterF.getData();
            Histogram f2hist = meterF.getHistogram2();
            double[] f2 = f2hist.getHistogram();
            for (int i=0; i<rdata.getLength(); i++) {
                double r = rdata.getValue(i);
                double mf = fdata.getValue(i);
                if (Double.isNaN(mf)) continue;
                double sdf = Math.sqrt(f2[i]-mf*mf);
                double pf = -sim.p2Truncated.du(r*r)/r;
                fw.write(String.format("%5.3f %22.15e %22.15e %22.15e\n", r, mf, sdf, pf));
            }
            fw.close();
        }

        long t2 = System.currentTimeMillis();
        System.out.println("time: "+(t2-t1)*0.001);
    }
    
    public static class LJMDParams extends ParameterBase {
        public int numAtoms = 100;
        public double temperature = 1.0;
        public double density = 0.01;
        public long numSteps = 1000000;
        public double rc = 4;
        public String functionsFile = null;
    }
}
