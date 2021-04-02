/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.dielectric;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.ChargeSource;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.meter.MeterDipoleSumSquared;
import etomica.data.meter.MeterDipoleSumSquaredMappedAverageFasterer;
import etomica.data.meter.MeterPotentialEnergyFromIntegratorFasterer;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMCFasterer;
import etomica.integrator.mcmove.MCMoveMoleculeFasterer;
import etomica.integrator.mcmove.MCMoveMoleculeRotateFasterer;
import etomica.lattice.LatticeCubicFcc;
import etomica.models.water.P2WaterTIP4P;
import etomica.models.water.P2WaterTIP4PSoft;
import etomica.models.water.SpeciesWater4P;
import etomica.molecule.DipoleSourceMolecular;
import etomica.molecule.DipoleSourceMolecularGeneric;
import etomica.potential.*;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeEwaldFourier;
import etomica.potential.ewald.P2Ewald1Real;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomNumberGenerator;

import java.awt.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;

/**
 * Canonical ensemble Monte Carlo simulation (NVT)
 * dielectric constant (epsilon)
 * TIP4P water
 *
 * @author Weisong
 */
public class TIP4P extends Simulation {
    protected final PotentialCompute potentialMaster;
    protected final IntegratorMCFasterer integrator;
    protected final MCMoveMoleculeFasterer moveMolecule;//translation
    protected final MCMoveMoleculeRotateFasterer rotateMolecule;//rotation
    protected final Box box;
    protected SpeciesGeneral species;
    protected P2WaterTIP4PSoft pWater;
    private final static String APP_NAME = "TIP4P water";
    private static final int PIXEL_SIZE = 15;

    protected double sigmaLJ, epsilonLJ;
    protected double chargeM, chargeH;

     //************************************* constructor ********************************************//
     public TIP4P(Space space, int numberMolecules, double dielectricOutside, double boxSize, double temperature, double truncation) {
         super(space);

    	 setRandom(new RandomNumberGenerator(2));

         species = SpeciesWater4P.create();
         addSpecies(species);
         box = this.makeBox();
         box.setNMolecules(species, numberMolecules);
         box.getBoundary().setBoxSize(Vector.of(new double[]{boxSize, boxSize, boxSize}));
         //	 double mu=168.96979945736229;//in simulation unit

         //for potential truncated

         AtomType oType = species.getTypeByName("O");
         AtomType hType = species.getTypeByName("H");
         AtomType mType = species.getTypeByName("M");

         sigmaLJ = P2WaterTIP4P.s;
         epsilonLJ = P2WaterTIP4P.e;
         chargeM = -2*P2WaterTIP4P.qH;
         chargeH = P2WaterTIP4P.qH;

         PotentialComputeEwaldFourier ewaldFourier = new PotentialComputeEwaldFourier(getSpeciesManager(), box, BondingInfo.noBonding());
         PotentialComputeEwaldFourier.EwaldParams params = ewaldFourier.getOptimalParamsForCutoff(3, truncation);

         ewaldFourier.setkCut(params.kCut);
         ewaldFourier.setCharge(mType, chargeM);
         ewaldFourier.setCharge(hType, chargeH);
         ewaldFourier.setAlpha(params.alpha);

         PotentialMasterFasterer pm = new PotentialMasterFasterer(getSpeciesManager(), box, BondingInfo.noBonding());

         TruncationFactory tf = new TruncationFactorySimple(space, truncation);
         Potential2Soft p2OO = tf.make(new P2LennardJones(space, sigmaLJ, epsilonLJ));
         Potential2Soft p2MM = tf.make(new P2Ewald1Real(chargeM*chargeM, params.alpha));
         Potential2Soft p2HH = tf.make(new P2Ewald1Real(chargeH*chargeH, params.alpha));
         P2HardGeneric p2MHC = P2HardSphere.makePotential(0.1);
         Potential2Soft p2MH = tf.make(p2MHC, new P2Ewald1Real(chargeH*chargeM, params.alpha));
         pm.setPairPotential(oType, oType, p2OO);
         pm.setPairPotential(mType, mType, p2MM);
         pm.setPairPotential(hType, hType, p2HH);
         pm.setPairPotential(mType, hType, p2MH);

         potentialMaster = new PotentialComputeAggregate(pm, ewaldFourier);

         // integrator from potential master
         integrator = new IntegratorMCFasterer(potentialMaster, random, temperature, box);
         // add mc move
         moveMolecule = new MCMoveMoleculeFasterer(random, potentialMaster, box);
         rotateMolecule = new MCMoveMoleculeRotateFasterer(random, potentialMaster, box);
         this.getController().addActivity(new ActivityIntegrate(integrator));
         //******************************** periodic boundary condition ******************************** //
         BoxImposePbc imposePbc = new BoxImposePbc(box, space);
         imposePbc.setApplyToMolecules(true);
         //**************************** integrator ****************************** //
         integrator.setTemperature(temperature);
         integrator.getMoveManager().addMCMove(moveMolecule);
         integrator.getMoveManager().addMCMove(rotateMolecule);
         integrator.getEventManager().addListener(new IntegratorListenerAction(imposePbc));

         //******************************** initial configuration ******************************** //
         LatticeCubicFcc lattice = new LatticeCubicFcc(space);
         ConfigurationLattice configuration = new ConfigurationLattice(lattice, space);
         configuration.initializeCoordinates(box);
     }

     // **************************** simulation part **************************** //
     public static void main (String[] args){
//    	 System.out.println(Electron.UNIT.toSim( 0.20819434)*1.855);
//    	 System.exit(2);
        Param params = new Param();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
        }
        final long startTime = System.currentTimeMillis();
        DateFormat date = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
        Calendar cal = Calendar.getInstance();
        System.out.println("startTime : " + date.format(cal.getTime()));
        Space space = Space3D.getInstance();
        int steps = params.steps;
        boolean isGraphic = params.isGraphic;
        boolean difInterval = params.difInterval;
        boolean mSquare = params.mSquare;
        boolean aEE = params.aEE;
        double temperature = Kelvin.UNIT.toSim(params.temperatureK);// convert Kelvin temperature to T(sim), essentially kT
//    	 System.out.println(temperature +" "+Kelvin.UNIT.toSim(1)+" "+Kelvin.UNIT.fromSim(1));
//    	 System.exit(2);
        int numberMolecules = params.numberMolecules;
        double density = params.density;//mol/L
        // double sigmaLJ=3.1535779419764953;
        // double epsilonLJ=64.8694333333333;
        // double chargeM = Electron.UNIT.toSim(-1.04);
        // double chargeH = Electron.UNIT.toSim(+0.52);
        // double hardCore = 0.7735839376326238;

        //double dipoleStrength=168.96979945736229;// in sim unit
        double dielectricOutside = params.dielectricOutside;
        double densitySim = density / (2 * Hydrogen.INSTANCE.getMass() + Oxygen.INSTANCE.getMass()) * Constants.AVOGADRO * 1e-24;  // convert to sim unit; in 1/(A)^3
//         System.out.println("Constants.AVOGADRO * 1e-27: "+Constants.AVOGADRO * 1e-27);
//         System.exit(2);
        double boxSize = Math.pow(numberMolecules / densitySim, (1.0 / 3.0));
        double truncation = boxSize * 0.49;
//         System.out.println("******************* TIP4P water, dielectric constant, NVT********************");
        System.out.println("steps = " + steps);
        System.out.println("numberMolecules = " + numberMolecules);
//         System.out.println("density= "+density+" mol/L");
//         System.out.println("denisty(sim)= "+densitySim + "1/(A)^3");
        System.out.println("density = " + density + " g/cm^3");
        System.out.println("temperature= " + params.temperatureK + " K");
//         System.out.println("temperature in sim unit = "+temperature);
//         System.out.println("box size= "+boxSize);
//         System.out.println("truncation= "+truncation);
//         System.out.println("dielectric constant outside = "+dielectricOutside);

        final TIP4P sim = new TIP4P(space, numberMolecules, dielectricOutside, boxSize, temperature, truncation);

//         System.out.println("*********** potential parameters:");

         double sigmaLJ=sim.sigmaLJ;
         double epsilonLJ=sim.epsilonLJ;
         double chargeM=sim.chargeM;
         double chargeH=sim.chargeH;


         double dipoleStrength = 168.96979945736229;

         if (isGraphic){
        	  SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
        	  simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(PIXEL_SIZE));
        	  simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

        	  ((DiameterHashByType)((DisplayBox)simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species.getTypeByName("H"),1);
        	  ((DiameterHashByType)((DisplayBox)simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species.getTypeByName("O"),1);
        	  ((DiameterHashByType)((DisplayBox)simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species.getTypeByName("M"),1);

            ColorSchemeByType colorScheme = (ColorSchemeByType) simGraphic.getDisplayBox(sim.box).getColorScheme();
            colorScheme.setColor(sim.species.getTypeByName("H"), Color.red);
            colorScheme.setColor(sim.species.getTypeByName("O"), Color.green);
            colorScheme.setColor(sim.species.getTypeByName("M"), Color.blue);

            simGraphic.makeAndDisplayFrame(APP_NAME);
            simGraphic.getDisplayBox(sim.box).repaint();
            return;
        }

        ////////////////////////////////////////////////////////////////////
        int blockNumber = 1000;
        int sampleAtInterval = numberMolecules;
        int samplePerBlock = steps / sampleAtInterval / blockNumber;
//         System.out.println("number of blocks is : "+blockNumber);
//         System.out.println("sample per block is : "+samplePerBlock);
        ////////////////////////////////////////////////////////////////////
         sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps / 5));//

         sim.integrator.getMoveManager().setEquilibrating(false);

//         MeterPotentialEnergyFasterer meterPE = new MeterPotentialEnergyFasterer(sim.potentialMaster);

         MeterPotentialEnergyFromIntegratorFasterer meterPE = new MeterPotentialEnergyFromIntegratorFasterer(sim.integrator);
         AccumulatorAverage accPE = new AccumulatorAverageFixed(samplePerBlock);
         DataPumpListener pumpPE = new DataPumpListener(meterPE, accPE, sampleAtInterval);
         sim.integrator.getEventManager().addListener(pumpPE);

         ChargeSource chargeSource = new ChargeSource(sim.getSpeciesManager());
         chargeSource.setCharge(sim.species.getTypeByName("M"), -2*P2WaterTIP4P.qH);
         chargeSource.setCharge(sim.species.getTypeByName("H"), P2WaterTIP4P.qH);
         DipoleSourceMolecular dipoleSource = new DipoleSourceMolecularGeneric(sim.box, chargeSource, null);

//         System.out.println("equilibration finished");
        // dipoleSumSquared
        MeterDipoleSumSquared dipoleMeter = null;
        AccumulatorAverage dipoleAccumulator = null;
        if (mSquare) {
            if (difInterval) {
                sampleAtInterval = numberMolecules / 40;
                samplePerBlock = steps / sampleAtInterval / blockNumber;
            }


            System.out.println("sampleperBlock " + samplePerBlock);
            System.out.println("SampleArInterval " + sampleAtInterval);
            System.out.println("BlockNumber " + (steps / sampleAtInterval / samplePerBlock));


            dipoleMeter = new MeterDipoleSumSquared(sim.box, dipoleSource);
            dipoleAccumulator = new AccumulatorAverageFixed(samplePerBlock);
            DataPump dipolePump = new DataPump(dipoleMeter, dipoleAccumulator);
            IntegratorListenerAction dipoleListener = new IntegratorListenerAction(dipolePump);
            dipoleListener.setInterval(sampleAtInterval);
            sim.integrator.getEventManager().addListener(dipoleListener);
        }
        // energy
//         MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
//         AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed(10);
//         DataPump energyPump = new DataPump(energyMeter, energyAccumulator);
//         energyAccumulator.setBlockSize(50);
//         IntegratorListenerAction energyListener = new IntegratorListenerAction(energyPump);
//         sim.integrator.getEventManager().addListener(energyListener);
        //externalField
//		MeterExternalFieldPerturbationWater meterExternalfiled =
//		new MeterExternalFieldPerturbationWater(space, sim.box,dipoleStrength,temperature, sim.potentialMaster);
//		AccumulatorAverage externalFieldAccumlator = new AccumulatorAverageFixed(samplePerBlock);
//		DataPumpListener externalFieldPumpListener = new DataPumpListener(meterExternalfiled,externalFieldAccumlator,sampleAtInterval);
//		sim.integrator.getEventManager().addListener(externalFieldPumpListener);

        //AEE

         MeterDipoleSumSquaredMappedAverageFasterer AEEMeter = new MeterDipoleSumSquaredMappedAverageFasterer(sim.box, dipoleStrength, temperature,
                 sim.potentialMaster, dipoleSource);
         AccumulatorAverageCovariance AEEAccumulator = new AccumulatorAverageCovariance(samplePerBlock, true);

        if (aEE) {
            if (difInterval) {
                sampleAtInterval = numberMolecules * 4;
                samplePerBlock = steps / sampleAtInterval / blockNumber;
            }
            System.out.println("sampleperBlock " + samplePerBlock);
            System.out.println("SampleArInterval " + sampleAtInterval);
            System.out.println("BlockNumber " + (steps / sampleAtInterval / samplePerBlock));

            AEEAccumulator = new AccumulatorAverageCovariance(samplePerBlock, true);
            DataPump AEEPump = new DataPump(AEEMeter, AEEAccumulator);
            IntegratorListenerAction AEEListener = new IntegratorListenerAction(AEEPump);
            AEEListener.setInterval(sampleAtInterval);
            //AEEListener.setInterval(1);//debug only to have more test samples
            sim.integrator.getEventManager().addListener(AEEListener);
        }
         sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));

        double avgPE = accPE.getData(accPE.AVERAGE).getValue(0);
        double errPE = accPE.getData(accPE.ERROR).getValue(0);
        System.out.println("PE: "+avgPE+" "+errPE);

        //calculate dipoleSumSquared average
        double dipoleSumSquared = 0;
        double dipoleSumSquaredERR = 0;
        double dipoleSumCor = 0;
        if (mSquare) {
            dipoleSumSquared = ((DataDouble) ((DataGroup) dipoleAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index)).x;
            dipoleSumSquaredERR = ((DataDouble) ((DataGroup) dipoleAccumulator.getData()).getData(AccumulatorAverage.ERROR.index)).x;
            dipoleSumCor = ((DataDouble) ((DataGroup) dipoleAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index)).x;
        }
//         //externalField
//         double UE = (((DataGroup)externalFieldAccumlator.getData()).getData(externalFieldAccumlator.AVERAGE.index)).getValue(0);
//         double UEERR =  (((DataGroup)externalFieldAccumlator.getData()).getData(externalFieldAccumlator.ERROR.index)).getValue(0);
//
//         double UEE = (((DataGroup)externalFieldAccumlator.getData()).getData(externalFieldAccumlator.AVERAGE.index)).getValue(1);
//         double UEEERR =  (((DataGroup)externalFieldAccumlator.getData()).getData(externalFieldAccumlator.ERROR.index)).getValue(1);
//
//         double UE2 = (((DataGroup)externalFieldAccumlator.getData()).getData(externalFieldAccumlator.AVERAGE.index)).getValue(2);
//         double UE2ERR =  (((DataGroup)externalFieldAccumlator.getData()).getData(externalFieldAccumlator.ERROR.index)).getValue(2);
//
//         double JE = (((DataGroup)externalFieldAccumlator.getData()).getData(externalFieldAccumlator.AVERAGE.index)).getValue(3);
//         double JEERR =  (((DataGroup)externalFieldAccumlator.getData()).getData(externalFieldAccumlator.ERROR.index)).getValue(3);
//
//         double JEE = (((DataGroup)externalFieldAccumlator.getData()).getData(externalFieldAccumlator.AVERAGE.index)).getValue(4);
//         double JEEERR =  (((DataGroup)externalFieldAccumlator.getData()).getData(externalFieldAccumlator.ERROR.index)).getValue(4);
//
//
//         System.out.println("UE = \t" + UE + " UEERR =  \t" + UEERR );
//         System.out.println("UEE = \t" + UEE + " UEEERR =  \t" + UEEERR);
//         System.out.println("UE2 = \t" + UE2 + " UE2RR =  \t" + UE2ERR );
//         System.out.println("JE = \t" + JE + " JEERR =  \t" + JEERR );
//         System.out.println("JEE = \t" + JEE + " JEERR =  \t" + JEEERR );

        //AEE
        double AEE = 0;
        double AEEER = 0;
        double AEECor = 0;
        if (aEE) {
            double sum0 = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index).getValue(0);
            double ERsum0 = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.ERROR.index).getValue(0);
            AEECor = ((DataGroup) AEEAccumulator.getData()).getData(AccumulatorAverage.BLOCK_CORRELATION.index).getValue(0);
            AEE = sum0;
            AEEER = ERsum0;
        }

        long endTime = System.currentTimeMillis();

        double totalTime = (endTime - startTime) / (1000.0 * 60.0);
        if (mSquare) {
            System.out.println("-<M^2>*bt*bt:\t" + (-dipoleSumSquared / temperature / temperature)
                    + " mSquareErr:\t" + (dipoleSumSquaredERR / temperature / temperature)
                    + " mSquareDifficulty:\t" + (dipoleSumSquaredERR / temperature / temperature) * Math.sqrt(totalTime)
                    + " dipolesumcor = " + dipoleSumCor);
            System.out.println("mSquare_Time: " + (endTime - startTime) / (1000.0 * 60.0));
        }
        if (aEE) {
            System.out.println("AEE_new:\t" + (AEE)
                    + " AEEErr:\t" + AEEER
                    + " AEEDifficulty:\t" + AEEER * Math.sqrt(totalTime)
                    + " AEECor = " + AEECor);
            System.out.println("AEE_Time: " + (endTime - startTime) / (1000.0 * 60.0));
        }
    }

    // ******************* parameters **********************//
    public static class Param extends ParameterBase {
        public boolean isGraphic = false;
        public boolean difInterval = !true;//TODO should be true when actually run with big system
        public boolean mSquare = true;
        public boolean aEE = false;
        public double temperatureK = 256;
        public int numberMolecules = 2;
        public double density = 0.1;//g/cm^3
        public double dielectricOutside = 1.0E11;
        public int steps = 100000;
    }
}
