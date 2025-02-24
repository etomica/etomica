/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterStructureFactor;
import etomica.graphics.*;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.LatticeCubicBcc;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeCubicSimple;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeField;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

import java.awt.*;

/**
 * Simple Lennard-Jones Monte Carlo simulation in 3D.
 */
public class LJMC3DFM extends Simulation {
    public IntegratorMC integrator;
//    public IntegratorVelocityVerlet integrator;
    public MCMoveAtom mcMoveAtom;
    public SpeciesGeneral species;
    public Box box;
    public PotentialComputePair potentialMaster, potentialMasterMeter;
    public PotentialComputePair potentialMasterSS;
    public PotentialComputeField pc1Meter, pc1;
    public P2LennardJones potential;
    public double temperature;

    public LJMC3DFM(int numAtoms, double rho, double temperature, double rc, double kSine, double lambda, Structure strc) {
        super(Space3D.getInstance());
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
//        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);
        box = this.makeBox();
        int cellRange = 2;
        NeighborCellManager neighborManager = new NeighborCellManager(getSpeciesManager(), box, cellRange, BondingInfo.noBonding());
//        NeighborListManager neighborManager = new NeighborListManager(getSpeciesManager(), box, 2, 4, BondingInfo.noBonding());
//        potentialMasterSS = new PotentialComputePair(getSpeciesManager(), box, neighborManager);
        potentialMaster = new PotentialComputePair(getSpeciesManager(), box, neighborManager);
        potentialMasterMeter = new PotentialComputePair(getSpeciesManager(), box, neighborManager);

        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(rho);
        inflater.actionPerformed();

        ConfigurationLattice config;
        String strcString;
        if (strc == Structure.SC) {
            config = new ConfigurationLattice(new LatticeCubicSimple(space), space);
            strcString = "SC";
        } else if (strc == Structure.BCC) {
            config = new ConfigurationLattice(new LatticeCubicBcc(space), space);
            strcString = "BCC";
        } else { //FCC
            config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
            strcString = "FCC";
        }

        config.initializeCoordinates(box);
        System.out.println(" L/2: " + box.getBoundary().getBoxSize().getX(0)/2.0);

        double sigma = 1.0;
        double eps = 1.0*lambda;
        potential = new P2LennardJones(sigma, eps);
        P2SoftSphericalTruncatedForceShifted potentialTruncated = new P2SoftSphericalTruncatedForceShifted(potential, rc);
        AtomType leafType = species.getLeafType();
        potentialMaster.setPairPotential(leafType, leafType, potentialTruncated);
        potentialMaster.doAllTruncationCorrection = false;

        P2LennardJones potentialMeter = new P2LennardJones(1,1);
        P2SoftSphericalTruncatedForceShifted potentialTruncatedMeter = new P2SoftSphericalTruncatedForceShifted(potentialMeter, rc);
        potentialMasterMeter.setPairPotential(leafType, leafType, potentialTruncatedMeter);
        potentialMasterMeter.doAllTruncationCorrection = false;

        double a;
        if (strc == Structure.SC) {
            a = Math.pow(1.0/rho, 1.0/3.0);
        } else if (strc == Structure.BCC) {
            a = Math.pow(2.0/rho, 1.0/3.0);
        } else {
            a = Math.pow(4.0/rho, 1.0/3.0);
        }

        Vector shift = space.makeVector();
        shift.E(box.getLeafList().get(0).getPosition());

//        P1Sinusoidal p1Sinusoidal = new P1Sinusoidal(getSpace(), a, lambda*kSine, strcString, shift);
        P1Sinusoidal p1Sinusoidal = new P1Sinusoidal(getSpace(), a, (1-lambda)*kSine, strcString, shift);
        pc1 = new PotentialComputeField(getSpeciesManager(), box);
        pc1.setFieldPotential(leafType, p1Sinusoidal);

        PotentialComputeAggregate.localStorageDefault = true;
        PotentialComputeAggregate pcAgg = new PotentialComputeAggregate(pc1, potentialMaster);

        integrator = new IntegratorMC(pcAgg, random, temperature, box);
        mcMoveAtom = new MCMoveAtom(random, pcAgg, box);
        integrator.getMoveManager().addMCMove(mcMoveAtom);

        P1Sinusoidal p1SinusoidalMeter = new P1Sinusoidal(getSpace(), a, kSine, strcString, shift);
        pc1Meter = new PotentialComputeField(getSpeciesManager(), box);
        pc1Meter.setFieldPotential(species.getLeafType(), p1SinusoidalMeter);

//        double r_nn = Double.NaN;
//        if (strc == Structure.SC ) r_nn = a;
//        if (strc == Structure.BCC) r_nn = a/Math.sqrt(2);
//        if (strc == Structure.FCC) r_nn = a/Math.sqrt(3);
//        double k1 = 2*Math.PI/r_nn;
//
//        MeterStructureFactor mSF = new MeterStructureFactor(box, k1*1.001);
//        mSF.getData();
//
//        for (int i=0;i<mSF.getPhaseAngles().length; i++) {
//            if(mSF.getWaveVectors()[i].squared() == k1*k1) {
//                double kx = mSF.getWaveVectors()[i].getX(0);
//                double ky = mSF.getWaveVectors()[i].getX(1);
//                double kz = mSF.getWaveVectors()[i].getX(2);
//                double k = Math.sqrt(mSF.getWaveVectors()[i].squared());
//                System.out.println(mSF.getPhaseAngles()[i] + " [ " + kx + " " + ky + " " + kz + " ]" + " " + k);
//            }
//        }
//        System.exit(0);
    }

    public static void main(String[] args) {
        long t1 = System.nanoTime();
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        int numAtoms = params.numAtoms;
        int numSteps = params.numSteps;
        double temperature = params.temperature;
        double rho = params.rho;
        double rc = params.rc;
        boolean isGraphic = params.isGraphic;
        Structure strc = params.strc;
        double alpha = params.alpha;
        double lambda = params.lambda;
        double kSine = params.kSine;
        double dbeta = params.dbeta;

        LJMC3DFM sim = new LJMC3DFM(numAtoms, rho, temperature, rc, kSine, lambda, strc);
        System.out.println(" LJ");
        System.out.println(" N: " + numAtoms);
        System.out.println(" density: " + rho);
        System.out.println(" T: " + temperature);
        System.out.println(" rc: " + rc);
        System.out.println(" steps: " +  numSteps);
        System.out.println(" kSine: " + kSine);
        System.out.println(" lambda: " + lambda);
        System.out.println(" Structure: " + strc);

        if (isGraphic) {
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.setPaintInterval(sim.box, 1000);
            ColorScheme colorScheme = new ColorScheme() {
                protected Color[] allColors;

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
            };
            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);

            ((DiameterHashByType)((DisplayBox)simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species().getAtomType(0),1.0);

            DisplayTextBox timer = new DisplayTextBox();
            DataSourceCountSteps counter = new DataSourceCountSteps(sim.integrator);
            DataPumpListener counterPump = new DataPumpListener(counter, timer, 100);
            sim.integrator.getEventManager().addListener(counterPump);
            simGraphic.getPanel().controlPanel.add(timer.graphic());

            DisplayPlot pePlot = new DisplayPlot();
            simGraphic.add(pePlot);

            simGraphic.makeAndDisplayFrame();
            return;
        }
        System.out.flush();

//        MeterFM meterFM = new MeterFM(sim.potentialMaster, temperature, sim.box, dbeta);
//        MeterFM meterFM = new MeterFM(sim.potentialMaster, sim.potentialMasterSS, temperature, sim.box, alpha, dbeta);

//        MeterPotentialEnergy meterPESinusoidal = new MeterPotentialEnergy(sim.pc1Meter);
//        MeterPotentialEnergy meterPELJ = new MeterPotentialEnergy(sim.potentialMaster);
//        MeterFEP meterFEP = new MeterFEP(sim.potentialMaster, 1.0/temperature);
//        MeterdUdlambda meterDUdlambda = new MeterdUdlambda(sim.potentialMasterMeter, sim.pc1Meter);

//        sim.integrator.reset();
//        System.out.println(" uLat: "+sim.integrator.getPotentialEnergy()/numAtoms);
//        meterFM.setUlat(sim.integrator.getPotentialEnergy());



        // Equilibaration
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps / 10));
        System.out.println(" Done with equilibration ...");



        int interval = numAtoms;
        int blocks = 100;
        long blockSize = params.numSteps / (interval * blocks);
        System.out.println(" numBlocks: " + blocks + " blocksize: " + blockSize + " interval: " + interval);
        System.out.println();

//        MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
//        AccumulatorAverageFixed energyAccumulator = new AccumulatorAverageFixed(blockSize);
//        DataPumpListener energyPump = new DataPumpListener(energyMeter, energyAccumulator, interval);
//        sim.integrator.getEventManager().addListener(energyPump);

//        AccumulatorAverageCovariance energyFMAccumulator = new AccumulatorAverageCovariance(blockSize);
//        DataPumpListener energyFMPump = new DataPumpListener(meterFM, energyFMAccumulator, interval);
//        sim.integrator.getEventManager().addListener(energyFMPump);

//        AccumulatorAverageFixed accumulator = new AccumulatorAverageFixed(blockSize);
//        DataPumpListener accumulatorPump = new DataPumpListener(meterPESinusoidal, accumulator, interval);
//        sim.integrator.getEventManager().addListener(accumulatorPump);


//        AccumulatorAverageFixed accumulatorLJ = new AccumulatorAverageFixed(blockSize);
//        DataPumpListener accumulatorPumpLJ = new DataPumpListener(meterPELJ, accumulatorLJ, interval);


//        if (lambda == 1) {
//            sim.integrator.getEventManager().addListener(accumulatorPumpLJ);
//        }

// Short sim
//        AccumulatorAverageFixed accumulatorFEP = new AccumulatorAverageFixed(blockSize);

//        if (lambda == 1) {
//            double EnShift = 0, errEnShift = 0;
//            long numStepsShort = numSteps / 10;
//            System.out.println(" Short sim for Covariance: " + numStepsShort + " steps");
//            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numStepsShort));
//            System.out.println(" Done with " + numStepsShort + " steps of short run-1");
//
//            double U0 = accumulatorLJ.getData(accumulatorLJ.AVERAGE).getValue(0);
//            double errU0 = accumulatorLJ.getData(accumulatorLJ.ERROR).getValue(0);
//            double corU0 = accumulatorLJ.getData(accumulatorLJ.BLOCK_CORRELATION).getValue(0);
//            System.out.println(" U0:  " + U0 + " " + errU0 + " " + corU0);
//            meterFEP.setU0(U0);
//
//            DataPumpListener accumulatorPumpFEP = new DataPumpListener(meterFEP, accumulatorFEP, interval);
//            sim.integrator.getEventManager().addListener(accumulatorPumpFEP);
//        }

//        AccumulatorAverageFixed accumulatorDUdlambda = new AccumulatorAverageFixed(blockSize);
//        DataPumpListener accumulatorPumpDUdlambda = new DataPumpListener(meterDUdlambda, accumulatorDUdlambda, interval);
//        sim.integrator.getEventManager().addListener(accumulatorPumpDUdlambda);
//
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps));

//        if (lambda == 1) {
//            double avgFEP = accumulatorFEP.getData(accumulatorFEP.AVERAGE).getValue(0);
//            double errFEP = accumulatorFEP.getData(accumulatorFEP.ERROR).getValue(0);
//            double corFEP = accumulatorFEP.getData(accumulatorFEP.BLOCK_CORRELATION).getValue(0);
//            System.out.println();
//            System.out.println(" Step1: LS+S -> S");
//            System.out.println(" dA1: " + -temperature*Math.log(avgFEP)/numAtoms + "  " + errFEP/avgFEP/numAtoms + "    " + corFEP);
//        }

//        double avg = accumulator.getData(accumulator.AVERAGE).getValue(0);
//        double err = accumulator.getData(accumulator.ERROR).getValue(0);
//        double cor = accumulator.getData(accumulator.BLOCK_CORRELATION).getValue(0);
//        System.out.println();
//        System.out.println(" Step2: LJ -> LS+S");
//        System.out.println(" dA2: " + avg/numAtoms + "  " + err/numAtoms + "    " + cor);



//        double avg = accumulatorDUdlambda.getData(accumulatorDUdlambda.AVERAGE).getValue(0);
//        double err = accumulatorDUdlambda.getData(accumulatorDUdlambda.ERROR).getValue(0);
//        double cor = accumulatorDUdlambda.getData(accumulatorDUdlambda.BLOCK_CORRELATION).getValue(0);
//        System.out.println();
//        System.out.println(" dUdl: " + avg/numAtoms + "  " + err/numAtoms + "    " + cor);



//        double avg2 = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(1);
//        double err2 = energyFMAccumulator.getData(energyFMAccumulator.ERROR).getValue(1);
//        double cor2 = energyFMAccumulator.getData(energyFMAccumulator.BLOCK_CORRELATION).getValue(1);
//
//        double avg3 = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(2);
//        double err3 = energyFMAccumulator.getData(energyFMAccumulator.ERROR).getValue(2);
//        double cor3 = energyFMAccumulator.getData(energyFMAccumulator.BLOCK_CORRELATION).getValue(2);

//        System.out.println(" m: " + avg3 + "  " + err3 + "    " + cor3);
//        System.out.println(" d: " + avg2 + "  " + err2 + "    " + cor2);
//        System.out.println(" errRatio: " + err1/err3);





//        double avg = energyAccumulator.getData(energyAccumulator.AVERAGE).getValue(0)/numAtoms;
//        double err = energyAccumulator.getData(energyAccumulator.ERROR).getValue(0)/numAtoms;
//        double cor = energyAccumulator.getData(energyAccumulator.BLOCK_CORRELATION).getValue(0);
//        System.out.println("e_int: " + avg + " " + err + " " + cor);


//        double avg_conv = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(0);
//        double err_conv = energyFMAccumulator.getData(energyFMAccumulator.ERROR).getValue(0);
//        double cor_conv = energyFMAccumulator.getData(energyFMAccumulator.BLOCK_CORRELATION).getValue(0);
//        System.out.println("e_conv: " + avg_conv + " " + err_conv + " " + cor_conv);
//
//        double avg_fm = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(1);
//        double err_fm = energyFMAccumulator.getData(energyFMAccumulator.ERROR).getValue(1);
//        double cor_fm = energyFMAccumulator.getData(energyFMAccumulator.BLOCK_CORRELATION).getValue(1);
//        System.out.println("e_fm: " + avg_fm + " " + err_fm + " " + cor_fm);
//
//        double avg_fm2 = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(2);
//        double err_fm2 = energyFMAccumulator.getData(energyFMAccumulator.ERROR).getValue(2);
//        double cor_fm2 = energyFMAccumulator.getData(energyFMAccumulator.BLOCK_CORRELATION).getValue(2);
//        System.out.println("diff: " + avg_fm2 + " " + err_fm2 + " " + cor_fm2);


//        double sumF2 = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(0);
//        double sumLambda = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(1);
//        double sumLambda2 = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(2);
//        double sumFHF = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(3);
//        System.out.println("sumF2  +  sumLambda  +  sumLambda2  +  sumFHF");
//        System.out.println(sumF2 +" "+ sumLambda +" "+ sumLambda2 +" "+ sumFHF);
//
//        double beta1 = 1.0/temperature + dbeta;
//        double alphaOpt = (beta1*sumF2-sumLambda)/(beta1*sumFHF+sumLambda2);
//        double beta0 = 1.0/temperature;
//        double dAlphaOpt = sumF2/(beta0*sumFHF+sumLambda2);
//        System.out.println("alphaOpt: " + alphaOpt);
//        System.out.println("dAlphaOpt: " + dAlphaOpt);
//
//
//
//        double avg_dphi = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(6);
//        double err_dphi = energyFMAccumulator.getData(energyFMAccumulator.ERROR).getValue(6);
//        double cor_dphi = energyFMAccumulator.getData(energyFMAccumulator.BLOCK_CORRELATION).getValue(6);
//        System.out.println("<dphi>/dalpha: " + avg_dphi + " " + err_dphi + " " + cor_dphi);
//        System.out.println();
//
//        x[0] = U0/numAtoms; // <U/N>
//        x[1] = sumF2/numAtoms; // <F2/N>
//        x[2] = sumHxx/numAtoms; // <trH/N>
//        x[3] = sumHxx/sumF2; //<trH/F2>



//        var_sumFdrHMAhalf = energyFMAccumulator.getData(energyFMAccumulator.COVARIANCE).getValue(18);
//        var_sumFdrHMA = energyFMAccumulator.getData(energyFMAccumulator.COVARIANCE).getValue(24);
//        cov_sumFdrHMA_sumFdrHMAhalf = energyFMAccumulator.getData(energyFMAccumulator.COVARIANCE).getValue(19);
//        corU_F2H = cov_sumFdrHMA_sumFdrHMAhalf/Math.sqrt(var_sumFdrHMAhalf*var_sumFdrHMA);
//        System.out.println(" cor: " + corU_F2H);
//        System.out.println(" opt: " + Math.sqrt(var_sumFdrHMAhalf/var_sumFdrHMA)*corU_F2H);


//        System.out.println("  u_hma: " + avg_fb + " " + err_fb + " " + cor_fb);
//        double varU = energyFMAccumulator.getData(energyFMAccumulator.COVARIANCE).getValue(0);
//        double varF2H = energyFMAccumulator.getData(energyFMAccumulator.COVARIANCE).getValue(4);
//
//        double covU_F2H = energyFMAccumulator.getData(energyFMAccumulator.COVARIANCE).getValue(1);
//        double corU_F2H = covU_F2H/Math.sqrt(varU*varF2H);
//        System.out.println();
//        System.out.println(" corU_F2H: " + corU_F2H);
//
//        System.out.println(err_u +" "+ err_F2 +" "+  corU_F2H);
//        System.out.println(" opt: " + err_u/err_F2*corU_F2H);



//        //www
//        double avg_h = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(2);
//        double err_h = energyFMAccumulator.getData(energyFMAccumulator.ERROR).getValue(2);
//        double cor_h = energyFMAccumulator.getData(energyFMAccumulator.BLOCK_CORRELATION).getValue(2);
//        System.out.println(" avg_h: " + avg_h + " " + err_h + " " + cor_h);
//
//        double avg_h_F2 = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(3);
//        double err_h_F2 = energyFMAccumulator.getData(energyFMAccumulator.ERROR).getValue(3);
//        double cor_h_F2 = energyFMAccumulator.getData(energyFMAccumulator.BLOCK_CORRELATION).getValue(3);
//        System.out.println(" avg_h_F2: " + avg_h + " " + err_h_F2 + " " + cor_h_F2);
//
//        double f1 = avg_F2/avg_h;
//        double f2 = 1/avg_h_F2;
//        System.out.println();
//        System.out.println("f1: " + f1);
//        System.out.println("f2: " + f2);
//        System.out.println();
//
//        double err_f2 = err_h_F2/avg_h_F2/avg_h_F2;
//        System.out.println(" alpha=sd_u/sd_f1: " + err_u/err_f2);
//        System.out.println();
//
//        double u_fb1 = avg_u - alpha*(f1 - temperature);
//        double u_fb2 = avg_u - alpha*(f2 - temperature);
//        System.out.println(" u_conv: " + avg_u + "  " + err_u + "  " + cor_u);
//        System.out.println("  u_fb1: " + u_fb1);
//        System.out.println("  u_fb2: " + u_fb2);
//
//        System.out.println();
//        double varU = energyFMAccumulator.getData(energyFMAccumulator.COVARIANCE).getValue(0);
//        double varF = energyFMAccumulator.getData(energyFMAccumulator.COVARIANCE).getValue(5);
//        double varH  = energyFMAccumulator.getData(energyFMAccumulator.COVARIANCE).getValue(10);
//        double varHF  = energyFMAccumulator.getData(energyFMAccumulator.COVARIANCE).getValue(15);
//
//        double covUF = energyFMAccumulator.getData(energyFMAccumulator.COVARIANCE).getValue(1);
//        double covUH = energyFMAccumulator.getData(energyFMAccumulator.COVARIANCE).getValue(2);
//        double covU_HF = energyFMAccumulator.getData(energyFMAccumulator.COVARIANCE).getValue(3);
//        double covFH = energyFMAccumulator.getData(energyFMAccumulator.COVARIANCE).getValue(6);
//        double covF_HF = energyFMAccumulator.getData(energyFMAccumulator.COVARIANCE).getValue(7);
//        double covH_HF = energyFMAccumulator.getData(energyFMAccumulator.COVARIANCE).getValue(11);
//
//        double corUF = covUF/Math.sqrt(varU*varF);
//        double corUH = covUH/Math.sqrt(varU*varH);
//        double corU_HF = covU_HF/Math.sqrt(varU*varHF);
//        double corFH = covFH/Math.sqrt(varF*varH);
//        double corF_HF = covF_HF/Math.sqrt(varF*varHF);
//        double corH_HF = covH_HF/Math.sqrt(varH*varHF);
//
//        System.out.println(" corUF: " + corUF);
//        System.out.println(" corUH: " + corUH);
//        System.out.println(" corU_HF: " + corU_HF);
//        System.out.println(" corFH: " + corFH);
//        System.out.println(" corF_HF: " + corF_HF);
//        System.out.println(" corH_HF: " + corH_HF);
//        www



//        double e_f = avg_conv - 1.0/2.0/alpha*temperature*avg_h + 1.0/2.0/alpha*avg_h_f2*avg_f2;
//        double avg_hma = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(8);
//        double err_hma = energyFMAccumulator.getData(energyFMAccumulator.ERROR).getValue(8);
//        double cor_hma = energyFMAccumulator.getData(energyFMAccumulator.BLOCK_CORRELATION).getValue(8);
//        System.out.println("  hma: " + avg_hma + " " + err_hma + " " + cor_hma);
//        System.out.println();
//
//        double avg = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(7);
//        double err = energyFMAccumulator.getData(energyFMAccumulator.ERROR).getValue(7);
//        double cor = energyFMAccumulator.getData(energyFMAccumulator.BLOCK_CORRELATION).getValue(7);
//        System.out.println("U0: " + avg + " " + err + " " + cor);

//        System.out.println();
//        double avg_u = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(0);
//        double err_u = energyFMAccumulator.getData(energyFMAccumulator.ERROR).getValue(0);
//        double cor_u = energyFMAccumulator.getData(energyFMAccumulator.BLOCK_CORRELATION).getValue(0);
//        System.out.println(" u_conv: " + avg_u + " " + err_u + " " + cor_u);
//
//        double avg_bf2 = energyFMAccumulator.getData(energyFMAccumulator.AVERAGE).getValue(1);
//        double err_bf2 = energyFMAccumulator.getData(energyFMAccumulator.ERROR).getValue(1);
//        double cor_bf2 = energyFMAccumulator.getData(energyFMAccumulator.BLOCK_CORRELATION).getValue(1);
//        System.out.println("   bF2: " + avg_bf2 + " " + err_bf2 + " " + cor_bf2);



//        double var_U = energyFMAccumulator.getData(energyFMAccumulator.COVARIANCE).getValue(0);
//        double var_bF2 = energyFMAccumulator.getData(energyFMAccumulator.COVARIANCE).getValue(4);
//        double var_H  = energyFMAccumulator.getData(energyFMAccumulator.COVARIANCE).getValue(8);
//
//        double covUbF2 = energyFMAccumulator.getData(energyFMAccumulator.COVARIANCE).getValue(1);
//        double covUH = energyFMAccumulator.getData(energyFMAccumulator.COVARIANCE).getValue(2);
//        System.out.println();
//        System.out.println(" cor_ubf2: " + covUbF2/Math.sqrt(var_U*var_bF2));
//        System.out.println(" cor_uh: " + covUH/Math.sqrt(var_U*var_H));
//
//        double zeroTest = covUH - covUbF2 + (avg_bf2+avg_bf2_short)*temperature/numAtoms;
//        System.out.println("covUH - covUbF2 + (avg_bf2+avg_bf2_short)*temperature/numAtoms");
//        System.out.println(covUH + "  " + covUbF2 + "  " + (avg_bf2+avg_bf2_short)*temperature/numAtoms + " ||| " + (covUbF2-(avg_bf2+avg_bf2_short)*temperature/numAtoms));
//        System.out.println("zeroTest: " + zeroTest);

        System.out.println("Move acceptance: " + sim.mcMoveAtom.getTracker().acceptanceProbability());
        long t2 = System.nanoTime();
        System.out.println("\ntime: " + (t2 - t1)/1.0e9/60.0 + " min");
    }

    public enum Structure {SC, BCC, FCC};

    public static class SimParams extends ParameterBase {
        //FCC: n=4-8: 256 500 864 1372 2048
        //FCC: n=5-10: 250 432 686 1024 1458 2000
        public Structure strc = Structure.SC;
        int n = 6;
        public int numAtoms = 1*n*n*n;
        public int numSteps = 1000000;
        public double temperature = 1;
        public double rho = 0.8;
        public double rc = 2.5;
        public boolean isGraphic = true;
        public double dbeta = 0.1;
        public double alpha = 1;

        public double kSine = 2110;
        public double lambda = 0.9998;
    }
}