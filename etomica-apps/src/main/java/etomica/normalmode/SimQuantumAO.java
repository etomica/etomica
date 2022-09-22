/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.ConformationGeneric;
import etomica.config.IConformation;
import etomica.data.*;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.potential.P1Anharmonic;
import etomica.potential.P1AnharmonicTIA;
import etomica.potential.P2Harmonic;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeField;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space1d.Space1D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
public class SimQuantumAO extends Simulation {
    public Box box;
    public IntegratorMC integrator;
    public MCMoveBox atomMove;
    public PotentialComputeField pcP1, pcP1EnTIA;
    public PotentialMasterBonding pmBonding;
    public PotentialCompute pm;
    public P1Anharmonic p1ah;
    public P1AnharmonicTIA p1ahUeff;
    public double betaN;
    public int nBeads;
    public double mass;
    public double k2_kin;
    public static final double kB = 1.0; // Constants.BOLTZMANN_K;
    public static final double hbar = 1.0;// Constants.PLANCK_H/(2.0*Math.PI);

    public SimQuantumAO(Space space, int nBeads, double temperature, double k2, double k4, double omega2, boolean isTIA, boolean moveReal) {
        super(space);
        Vector[] initCoords = new Vector[nBeads];
        for (int i = 0; i < nBeads; i++){
            initCoords[i] = space.makeVector();
        }
        IConformation conformation = new ConformationGeneric(initCoords);
        SpeciesGeneral species = new SpeciesBuilder(space)
                .addCount(AtomType.simpleFromSim(this), nBeads)
                .withConformation(conformation)
                .build();
        addSpecies(species);

        box = this.makeBox(new BoundaryRectangularNonperiodic(space));
        box.setNMolecules(species, 1);
        mass = species.getLeafType().getMass();

        //pm2 that uses the full PI potential, for data collection
        //spring P2 part (x_i-x_{i+1})^2
        pmBonding = new PotentialMasterBonding(getSpeciesManager(), box);

        double beta = 1.0/(kB*temperature);
        betaN = beta/nBeads;
        double omegaN = 1.0/(hbar*betaN);

        k2_kin = nBeads == 1 ? 0 : (mass*omegaN*omegaN/nBeads);

        P2Harmonic p2Bond = new P2Harmonic(k2_kin, 0);
        List<int[]> pairs = new ArrayList<>();
        for (int i=0; i<nBeads; i++) {
            int[] p = new int[]{i,(i+1)%nBeads};
            pairs.add(p);
        }
        pmBonding.setBondingPotentialPair(species, p2Bond, pairs);

        pcP1 = new PotentialComputeField(getSpeciesManager(), box);
        if (isTIA){
            double facUeff = 1.0;
            p1ahUeff = new P1AnharmonicTIA(space, k2, k4, nBeads, mass*omegaN*omegaN, facUeff);
            pcP1.setFieldPotential(species.getLeafType(), p1ahUeff);
        } else {
            p1ah = new P1Anharmonic(space, k2, k4, nBeads);
            pcP1.setFieldPotential(species.getLeafType(), p1ah);
        }
        pm = new PotentialComputeAggregate(pmBonding, pcP1);

        integrator = new IntegratorMC(pm, random, temperature, box);
        if (moveReal) {
            atomMove = new MCMoveHOReal(space, pm, random, temperature, omega2, box);
        }
        else {
            atomMove = new MCMoveHO(space, pm, random, temperature, omega2, box);
        }
        integrator.getMoveManager().addMCMove(atomMove);

        double facEn = 3.0;
        P1AnharmonicTIA p1ahEn = new P1AnharmonicTIA(space, k2, k4, nBeads, mass*omegaN*omegaN, facEn);
        pcP1EnTIA = new PotentialComputeField(getSpeciesManager(), box);
        pcP1EnTIA.setFieldPotential(species.getLeafType(), p1ahEn);
    }

    public Integrator getIntegrator() {
        return integrator;
    }

    public static void main(String[] args) {

        OctaneParams params = new OctaneParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            // custom parameters
            params.numSteps = 10000000;
            params.temperature = 1.0;
            params.nBeads = 33;
            params.k2 = 1.0;
            params.k4 = 24.0;
            params.moveReal = !true;
        }

        double temperature = params.temperature;
        double k2 = params.k2;
        double k4 = params.k4;
        int nBeads = params.nBeads;
        boolean graphics = params.graphics;
        long numSteps = params.numSteps;
        long numStepsEqu = numSteps/20;
        boolean isTIA = params.isTIA;
        boolean moveReal = params.moveReal;

        double omegaN = nBeads*kB*temperature/hbar; // 1/hbar*betan

        double massssss = 1.0;
        double omega2 = k2/massssss;
        if (isTIA){
            omega2 = omega2*(1.0 + omega2/12.0/omegaN/omegaN);
        }

        final SimQuantumAO sim = new SimQuantumAO(Space1D.getInstance(), nBeads, temperature, k2, k4, omega2, isTIA, moveReal);
        sim.integrator.reset();

        System.out.println(" numSteps: " +  numSteps + " numStepsEqu: " + numStepsEqu);
        System.out.println(" nBeads: " + nBeads);
        System.out.println(" temperature: " + temperature);
        System.out.println(" mass: " + sim.mass + " k2: " + k2 + " k4: " + k4);
        System.out.println(" hbar: " + hbar  + "  kB: " + kB);
        System.out.println(" k2_kin: " + sim.k2_kin);
        System.out.println(" isTIA: " + isTIA);

        MeterMSDHO meterMSDHO = new MeterMSDHO(nBeads, sim.box);

        MeterPIPrim meterPrim = null;
        MeterPIVir meterVir = null;
        MeterPICentVir meterCentVir = null;
        MeterPIHMAc meterHMAc = null;
        MeterPIHMA meterHMA = null;
        if (isTIA){
//            meterPrim = new MeterPIPrim(sim.pmBonding, sim.pcP1EnTIA, nBeads, sim.betaN);
//            meterVir = new MeterPIVirTIA(sim.pcP1EnTIA, sim.pcP1, sim.betaN, nBeads, sim.box);
//            meterCentVir = new MeterPICentVirTIA(sim.pcP1EnTIA, sim.pcP1, sim.betaN, nBeads, sim.box);
            meterHMAc = null;
//            meterHMA = new MeterPIHMATIA(sim.pmBonding, sim.pcP1EnTIA, sim.pcP1, sim.betaN, nBeads, omega2, sim.box);
//            meterHMAcent = null;
        } else {
            meterPrim = new MeterPIPrim(sim.pmBonding, sim.pcP1, nBeads, sim.betaN);
            meterVir = new MeterPIVir(sim.pcP1, sim.betaN, nBeads, sim.box);
            meterCentVir = new MeterPICentVir(sim.pcP1, sim.betaN, nBeads, sim.box);
            meterHMAc = new MeterPIHMAc(sim.pcP1, sim.betaN, nBeads, sim.box);
            meterHMA = new MeterPIHMA(sim.pmBonding, sim.pcP1, sim.betaN, nBeads, omega2, sim.box);
        }


        MeterPIVirMidPt meterCentVirBar = new MeterPIVirMidPt(sim.pcP1, sim.betaN, nBeads, sim.box); //Bad!!
        MeterPIHMAvir meterHMAvir = new MeterPIHMAvir(sim.pmBonding, sim.pcP1, sim.betaN, nBeads, omega2, sim.box);//Bad!!

        if (graphics) {
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
            DisplayTextBox timer = new DisplayTextBox();
            DataSourceCountSteps counter = new DataSourceCountSteps(sim.integrator);
            DataPumpListener counterPump = new DataPumpListener(counter, timer, 100);
            sim.integrator.getEventManager().addListener(counterPump);
            simGraphic.getPanel().controlPanel.add(timer.graphic());
            simGraphic.makeAndDisplayFrame(" PIMC ");
            return;
        }

        System.out.flush();
        final long startTime = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numStepsEqu));
        sim.integrator.getMoveManager().setEquilibrating(false);

        int numBlocks = 1000;
        int interval = 1;
        long blockSize = numSteps/numBlocks;
        if (blockSize == 0) blockSize = 1;
        System.out.println(" numBlocks: " + numBlocks + " blocksize: " + blockSize + " interval: " + interval);

        AccumulatorAverageFixed accumulatorMSD = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPumpMSD = new DataPumpListener(meterMSDHO, accumulatorMSD, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpMSD);

        //1 Primitive
        AccumulatorAverageCovariance accumulatorPrim = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener accumulatorPumpPrim = new DataPumpListener(meterPrim, accumulatorPrim, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpPrim);

        //2 Virial
        AccumulatorAverageCovariance accumulatorVir = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener accumulatorPumpVir = new DataPumpListener(meterVir, accumulatorVir, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpVir);

        //3 Centroid Virial
        AccumulatorAverageCovariance accumulatorCentVir = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener accumulatorPumpCentVir = new DataPumpListener(meterCentVir, accumulatorCentVir, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpCentVir);

         //4 HMAc (CLassical EC)
        AccumulatorAverageCovariance accumulatorHMAc = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener accumulatorPumpHMAc = new DataPumpListener(meterHMAc, accumulatorHMAc, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpHMAc);

        //5 HMAq (Quantum EC)
        AccumulatorAverageCovariance accumulatorHMA = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener accumulatorPumpHMA = new DataPumpListener(meterHMA, accumulatorHMA, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpHMA);

        long Nshort = numSteps/10;
        System.out.println(" N_short_sim = " + Nshort);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, Nshort));
        DataGroup dataHMAShort = (DataGroup) accumulatorHMA.getData();
        IData dataHMAShortAvg = dataHMAShort.getData(accumulatorHMA.AVERAGE.index);
        IData dataHMAShortErr = dataHMAShort.getData(accumulatorHMA.ERROR.index);
        double EnShift = dataHMAShortAvg.getValue(0);
        double errEnShift = dataHMAShortErr.getValue(0);
        System.out.println(" EnShift: " + EnShift + " +/- " + errEnShift);
        meterPrim.setShift(EnShift);
        meterVir.setShift(EnShift);
        meterCentVir.setShift(EnShift);
        meterHMAc.setShift(EnShift);
        meterHMA.setShift(EnShift);


        accumulatorPrim.reset();
        accumulatorVir.reset();
        accumulatorCentVir.reset();
        accumulatorHMAc.reset();
        accumulatorHMA.reset();




        //run
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));

        //MSD
        DataGroup dataMSD = (DataGroup)accumulatorMSD.getData();
        IData dataMSDAvg = dataMSD.getData(accumulatorMSD.AVERAGE.index);
        IData dataMSDErr = dataMSD.getData(accumulatorMSD.ERROR.index);
        IData dataMSDCorrelation = dataMSD.getData(accumulatorMSD.BLOCK_CORRELATION.index);
        double avgMSD = dataMSDAvg.getValue(0);
        double errMSD = dataMSDErr.getValue(0);
        double corMSD = dataMSDCorrelation.getValue(0);

        //Prim
        DataGroup dataPrim = (DataGroup)accumulatorPrim.getData();
        IData dataAvgPrim = dataPrim.getData(accumulatorPrim.AVERAGE.index);
        IData dataErrPrim = dataPrim.getData(accumulatorPrim.ERROR.index);
        IData dataCorPrim = dataPrim.getData(accumulatorPrim.BLOCK_CORRELATION.index);
        IData dataCovPrim  = dataPrim.getData(accumulatorPrim.COVARIANCE.index);

        double avgEnPrim = dataAvgPrim.getValue(0) + EnShift;
        double errEnPrim = dataErrPrim.getValue(0);
        double corEnPrim = dataCorPrim.getValue(0);
        System.out.println("\n En_prim: " + avgEnPrim  + " +/- " + errEnPrim + " cor: " + corEnPrim);


        //Vir
        DataGroup dataVir = (DataGroup)accumulatorVir.getData();
        IData dataAvgVir = dataVir.getData(accumulatorVir.AVERAGE.index);
        IData dataErrVir = dataVir.getData(accumulatorVir.ERROR.index);
        IData dataCorVir = dataVir.getData(accumulatorVir.BLOCK_CORRELATION.index);
        IData dataCovVir  = dataVir.getData(accumulatorVir.COVARIANCE.index);

        double avgEnVir = dataAvgVir.getValue(0) + EnShift;
        double errEnVir = dataErrVir.getValue(0);
        double corEnVir = dataCorVir.getValue(0);
        System.out.println(" En_vir:  " + avgEnVir  + " +/- " + errEnVir + " cor: " + corEnVir);

        //Cent-Vir
        DataGroup dataCentVir = (DataGroup)accumulatorCentVir.getData();
        IData dataAvgCentVir = dataCentVir.getData(accumulatorCentVir.AVERAGE.index);
        IData dataErrCentVir = dataCentVir.getData(accumulatorCentVir.ERROR.index);
        IData dataCorCentVir = dataCentVir.getData(accumulatorCentVir.BLOCK_CORRELATION.index);
        IData dataCovCentVir  = dataCentVir.getData(accumulatorCentVir.COVARIANCE.index);

        double avgEnCentVir = dataAvgCentVir.getValue(0) + EnShift;
        double errEnCentVir = dataErrCentVir.getValue(0);
        double corEnCentVir = dataCorCentVir.getValue(0);
        System.out.println(" En_cvir: " + avgEnCentVir  + " +/- " + errEnCentVir + " cor: " + corEnCentVir);

        //HMAc
        DataGroup dataHMAc = (DataGroup)accumulatorHMAc.getData();
        IData dataAvgHMAc = dataHMAc.getData(accumulatorHMAc.AVERAGE.index);
        IData dataErrHMAc = dataHMAc.getData(accumulatorHMAc.ERROR.index);
        IData dataCorHMAc = dataHMAc.getData(accumulatorHMAc.BLOCK_CORRELATION.index);
        IData dataCovHMAc  = dataHMAc.getData(accumulatorHMAc.COVARIANCE.index);

        double avgEnHMAc = dataAvgHMAc.getValue(0) + EnShift;
        double errEnHMAc = dataErrHMAc.getValue(0);
        double corEnHMAc = dataCorHMAc.getValue(0);
        System.out.println(" En_hmac: " + avgEnHMAc  + " +/- " + errEnHMAc + " cor: " + corEnHMAc);

        //HMA
        DataGroup dataHMA = (DataGroup)accumulatorHMA.getData();
        IData dataAvgHMA = dataHMA.getData(accumulatorHMA.AVERAGE.index);
        IData dataErrHMA = dataHMA.getData(accumulatorHMA.ERROR.index);
        IData dataCorHMA = dataHMA.getData(accumulatorHMA.BLOCK_CORRELATION.index);
        IData dataCovHMA = dataHMA.getData(accumulatorHMA.COVARIANCE.index);

        double avgEnHMA = dataAvgHMA.getValue(0) + EnShift;
        double errEnHMA = dataErrHMA.getValue(0);
        double corEnHMA = dataCorHMA.getValue(0);
        System.out.println(" En_hma:  " + avgEnHMA  + " +/- " + errEnHMA + " cor: " + corEnHMA);


        System.out.println("\n Quantum Harmonic Oscillator Theory");
        double omega = Math.sqrt(omega2);
        System.out.println(" ====================================");
        System.out.println(" MSD_sim: " + avgMSD + " +/- " + errMSD + " cor: " + corMSD);
        System.out.println(" MSDc: " + sim.kB*temperature/ sim.mass/omega2);
        System.out.println(" MSDq: " + hbar/sim.mass/omega*(0.5+1.0/(Math.exp(hbar*omega/temperature)-1.0)));

        double EnQ = hbar*omega*(0.5 + 1/(Math.exp(nBeads*sim.betaN*hbar*omega)-1.0));
        double EnC = temperature;
        System.out.println("\n EnC: " + EnC);
        System.out.println(" EnQ: " + EnQ);




        System.out.println("\n ********** Cvn ***********");
        double kB_beta2 = kB*sim.betaN*sim.betaN*nBeads*nBeads;
        double CvnPrim = kB_beta2*(dataAvgPrim.getValue(1) + dataCovPrim.getValue(0));
        double CvnVir = kB_beta2*(dataAvgVir.getValue(1) + dataCovVir.getValue(0));
        double CvnCentVir = kB_beta2*(dataAvgCentVir.getValue(1) + dataCovCentVir.getValue(0));
        double CvnHMAc = kB_beta2*(dataAvgHMAc.getValue(1) + dataCovHMAc.getValue(0));
        double CvnHMA  = kB_beta2*(dataAvgHMA.getValue(1) + dataCovHMA.getValue(0));


        System.out.println(" Cvn_prim: " + CvnPrim);
        System.out.println(" Cvn_vir: " + CvnVir);
        System.out.println(" Cvn_cvir: " + CvnCentVir);
        System.out.println(" Cvn_hmac: " + CvnHMAc);
        System.out.println(" Cvn_hma: " + CvnHMA);


        //Acceptance ratio
        System.out.println("\n acceptance %: " + 100*sim.atomMove.getTracker().acceptanceRatio());

        long endTime = System.currentTimeMillis();
        System.out.println("\n time (min): " + (endTime - startTime)/60.0/1000.0);
    }

    public static class OctaneParams extends ParameterBase {
        public double temperature = 1.0;
        public int nBeads = 11;
        public boolean graphics = false;
        public double k2 = 1.0;
        public double k4 = 24.0;
        public long numSteps = 1_000_000;
        public boolean isTIA = false;
        public boolean moveReal = !true;
    }
}
