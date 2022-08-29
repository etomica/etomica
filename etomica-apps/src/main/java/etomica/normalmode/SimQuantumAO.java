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
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.potential.*;
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

import java.awt.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
public class SimQuantumAO extends Simulation {

    public SpeciesGeneral species;
    public Box box;
    public IntegratorMC integrator;
    public MCMoveBox atomMove;
    public PotentialComputeField pcP1, pcP1EnTIA;
    public PotentialMasterBonding pmBonding;
    public PotentialCompute pm;
    public double betaN;
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
        box.getBoundary().setBoxSize(Vector.of(new double[]{1}));
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


        IPotential1 p1ah, p1ahUeff;
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

    public static void main(String[] args) {

        OctaneParams params = new OctaneParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            // custom parameters
            params.numSteps = 100000;
            params.temperature = 1;
            params.nBeads = 256;
            params.k4 = 0;
            params.k2 = 1;
            params.moveReal = true;
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
        DataSourceScalar meterPrim, meterVir, meterCentVir, meterCentVirMod, meterHMA, meterHMAcent;

        if (isTIA){
            meterPrim = new MeterPIPrim(sim.pmBonding, sim.pcP1EnTIA, sim.betaN);
            meterVir = new MeterPIVirTIA(sim.pcP1EnTIA, sim.pcP1, sim.betaN, nBeads, sim.box);
            meterCentVir = new MeterPICentVirTIA(sim.pcP1EnTIA, sim.pcP1, sim.betaN, nBeads, sim.box);
            meterCentVirMod = null;
            meterHMA = new MeterPIHMATIA(sim.pmBonding, sim.pcP1EnTIA, sim.pcP1, sim.betaN, nBeads, omega2, sim.box);
            meterHMAcent = null;
        } else {
            meterPrim = new MeterPIPrim(sim.pmBonding, sim.pcP1, sim.betaN);
            meterVir = new MeterPIVir(sim.pcP1, sim.betaN, nBeads, sim.box);
            meterCentVir = new MeterPICentVir(sim.pcP1, sim.betaN, nBeads, sim.box);
            meterCentVirMod = new MeterPICentVirMod(sim.pcP1, sim.betaN, nBeads, sim.box);
            meterHMA = new MeterPIHMA(sim.pmBonding, sim.pcP1, sim.betaN, nBeads, omega2, sim.box);
            meterHMAcent = null; // new MeterPIHMAcent(sim.pmBonding, sim.pcP1, sim.betaN, nBeads, omega2, sim.box);
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

        // Primitive
        AccumulatorAverageFixed accumulatorPrim = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPumpPrim = new DataPumpListener(meterPrim, accumulatorPrim, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpPrim);

        // Virial
        AccumulatorAverageFixed accumulatorVir = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPumpVir = new DataPumpListener(meterVir, accumulatorVir, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpVir);

        // Centroid Virial
        AccumulatorAverageFixed accumulatorCentVir = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPumpCentVir = new DataPumpListener(meterCentVir, accumulatorCentVir, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpCentVir);

        // Modified Centroid Virial
        AccumulatorAverageFixed accumulatorCentVirMod = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPumpCentVirMod = new DataPumpListener(meterCentVirMod, accumulatorCentVirMod, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpCentVirMod);

        // Virial-bar
        AccumulatorAverageFixed accumulatorVirBar = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPumpVirBar = new DataPumpListener(meterCentVirBar, accumulatorVirBar, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpVirBar);


        // HMA
        AccumulatorAverageFixed accumulatorHMA = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPumpHMA = new DataPumpListener(meterHMA, accumulatorHMA, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpHMA);

        // HMAvir
        AccumulatorAverageFixed accumulatorHMAvir = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPumpHMAvir = new DataPumpListener(meterHMAvir, accumulatorHMAvir, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpHMAvir);

        // HMA-cent
        AccumulatorAverageFixed accumulatorHMAcent = null;
        if (meterHMAcent != null) {
            accumulatorHMAcent =  new AccumulatorAverageFixed(blockSize);
            DataPumpListener accumulatorPumpHMAcent = new DataPumpListener(meterHMAcent, accumulatorHMAcent, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpHMAcent);
        }

        //run
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));

        //MSD
        DataGroup dataMSD = (DataGroup)accumulatorMSD.getData();
        IData dataMSDErr = dataMSD.getData(accumulatorMSD.ERROR.index);
        IData dataMSDAvg = dataMSD.getData(accumulatorMSD.AVERAGE.index);
        IData dataMSDCorrelation = dataMSD.getData(accumulatorMSD.BLOCK_CORRELATION.index);
        double avgMSD = dataMSDAvg.getValue(0);
        double errMSD = dataMSDErr.getValue(0);
        double corMSD = dataMSDCorrelation.getValue(0);

        //Prim
        DataGroup dataPrim = (DataGroup)accumulatorPrim.getData();
        IData dataErrPrim = dataPrim.getData(accumulatorPrim.ERROR.index);
        IData dataAvgPrim = dataPrim.getData(accumulatorPrim.AVERAGE.index);
        IData dataCorrelationPrim = dataPrim.getData(accumulatorPrim.BLOCK_CORRELATION.index);
        double avgEnPrim = dataAvgPrim.getValue(0);
        double errEnPrim = dataErrPrim.getValue(0);
        double corEnPrim = dataCorrelationPrim.getValue(0);
        System.out.println("\n En_prim: " + avgEnPrim  + " +/- " + errEnPrim + " cor: " + corEnPrim);

        //Vir
        DataGroup dataVir = (DataGroup)accumulatorVir.getData();
        IData dataErrVir = dataVir.getData(accumulatorVir.ERROR.index);
        IData dataAvgVir = dataVir.getData(accumulatorVir.AVERAGE.index);
        IData dataCorrelationVir = dataVir.getData(accumulatorVir.BLOCK_CORRELATION.index);
        double avgEnVir = dataAvgVir.getValue(0);
        double errEnVir = dataErrVir.getValue(0);
        double corEnVir = dataCorrelationVir.getValue(0);
        System.out.println(" En_vir:  " + avgEnVir  + " +/- " + errEnVir + " cor: " + corEnVir);

        //Cent Vir
        DataGroup dataCentVir = (DataGroup)accumulatorCentVir.getData();
        IData dataErrCentVir = dataCentVir.getData(accumulatorCentVir.ERROR.index);
        IData dataAvgCentVir = dataCentVir.getData(accumulatorCentVir.AVERAGE.index);
        IData dataCorrelationCentVir = dataCentVir.getData(accumulatorCentVir.BLOCK_CORRELATION.index);
        double avgEnCentVir = dataAvgCentVir.getValue(0);
        double errEnCentVir = dataErrCentVir.getValue(0);
        double corEnCentVir = dataCorrelationCentVir.getValue(0);
        System.out.println(" En_cvir: " + avgEnCentVir  + " +/- " + errEnCentVir + " cor: " + corEnCentVir);

        //Modified Cent Vir
        DataGroup dataCentVirMod = (DataGroup)accumulatorCentVirMod.getData();
        IData dataErrCentVirMod = dataCentVirMod.getData(accumulatorCentVirMod.ERROR.index);
        IData dataAvgCentVirMod = dataCentVirMod.getData(accumulatorCentVirMod.AVERAGE.index);
        IData dataCorrelationCentVirMod = dataCentVirMod.getData(accumulatorCentVirMod.BLOCK_CORRELATION.index);
        double avgEnCentVirMod = dataAvgCentVirMod.getValue(0);
        double errEnCentVirMod = dataErrCentVirMod.getValue(0);
        double corEnCentVirMod = dataCorrelationCentVirMod.getValue(0);
        System.out.println(" En_modcvir: " + avgEnCentVirMod  + " +/- " + errEnCentVirMod + " cor: " + corEnCentVirMod);


//        //Vir-bar
//        DataGroup dataVirBar = (DataGroup)accumulatorVirBar.getData();
//        IData dataErrVirBar = dataVirBar.getData(accumulatorVirBar.ERROR.index);
//        IData dataAvgVirBar = dataVirBar.getData(accumulatorVirBar.AVERAGE.index);
//        IData dataCorrelationVirBar = dataVirBar.getData(accumulatorVirBar.BLOCK_CORRELATION.index);
//        double avgEnVirBar = dataAvgVirBar.getValue(0);
//        double errEnVirBar = dataErrVirBar.getValue(0);
//        double corEnVirBar = dataCorrelationVirBar.getValue(0);
//        System.out.println(" En_vir_bar:  " + avgEnVirBar  + " +/- " + errEnVirBar + " cor: " + corEnVirBar);


        //HMA
        DataGroup dataHMA = (DataGroup)accumulatorHMA.getData();
        IData dataErrHMA = dataHMA.getData(accumulatorHMA.ERROR.index);
        IData dataAvgHMA = dataHMA.getData(accumulatorHMA.AVERAGE.index);
        IData dataCorrelationHMA = dataHMA.getData(accumulatorHMA.BLOCK_CORRELATION.index);
        double avgEnHMA = dataAvgHMA.getValue(0);
        double errEnHMA = dataErrHMA.getValue(0);
        double corEnHMA = dataCorrelationHMA.getValue(0);
        System.out.println(" En_hma:  " + avgEnHMA  + " +/- " + errEnHMA + " cor: " + corEnHMA);


//        //HMA-vir
//        DataGroup dataHMAvir = (DataGroup)accumulatorHMAvir.getData();
//        IData dataErrHMAvir = dataHMAvir.getData(accumulatorHMAvir.ERROR.index);
//        IData dataAvgHMAvir = dataHMAvir.getData(accumulatorHMAvir.AVERAGE.index);
//        IData dataCorrelationHMAvir = dataHMAvir.getData(accumulatorHMAvir.BLOCK_CORRELATION.index);
//        double avgEnHMAvir = dataAvgHMAvir.getValue(0);
//        double errEnHMAvir = dataErrHMAvir.getValue(0);
//        double corEnHMAvir = dataCorrelationHMAvir.getValue(0);
//        System.out.println(" En_hmavir:  " + avgEnHMAvir  + " +/- " + errEnHMAvir + " cor: " + corEnHMAvir);


        //HMA-cent
        if (accumulatorHMAcent != null) {
            DataGroup dataHMAcent = (DataGroup)accumulatorHMAcent.getData();
            IData dataErrHMAcent = dataHMAcent.getData(accumulatorHMAcent.ERROR.index);
            IData dataAvgHMAcent = dataHMAcent.getData(accumulatorHMAcent.AVERAGE.index);
            IData dataCorrelationHMAcent = dataHMAcent.getData(accumulatorHMAcent.BLOCK_CORRELATION.index);
            double avgEnHMAcent = dataAvgHMAcent.getValue(0);
            double errEnHMAcent = dataErrHMAcent.getValue(0);
            double corEnHMAcent = dataCorrelationHMAcent.getValue(0);
            System.out.println(" En_hmacent:  " + avgEnHMAcent  + " +/- " + errEnHMAcent + " cor: " + corEnHMAcent);
        }

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

        //Acceptance ratio
        System.out.println("\n acceptance %: " + 100*sim.atomMove.getTracker().acceptanceRatio());

        long endTime = System.currentTimeMillis();
        System.out.println("\n time (min): " + (endTime - startTime)/60.0/1000.0);
    }

    public static class OctaneParams extends ParameterBase {
        public double temperature = 1.0;
        public int nBeads = 11; //must be odd for now!
        public boolean graphics = false;
        public double k2 = 1.0;
        public double k4 = 24.0;
        public long numSteps = 1_000_000;
        public boolean isTIA = false;
        public boolean moveReal = false;
    }
}