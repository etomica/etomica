/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedDensity;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterProfileByVolume;
import etomica.data.types.DataFunction;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.math.function.FunctionDifferentiable;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.util.List;

/**
 * Simple Simulation for Monte Carlo simulation of Lennard-Jonesium.  Interactions are
 * tracked with cell lists.
 *
 * @author Andrew Schultz
 */

public class MappedSingletDensity extends Simulation {

    public final PotentialMasterCell potentialMaster;
    public final ActivityIntegrate activityIntegrate;
    public final IntegratorMC integrator;

    /**
     * Creates simulation with default parameters from {@link SimParams}
     */
    public MappedSingletDensity() {
        this(new SimParams());
    }

    /**
     * Creates simulation with the given parameters
     */
    public MappedSingletDensity(SimParams params) {
        super(Space3D.getInstance());

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        Box box = new Box(space);
        addBox(box);
        box.setNMolecules(species, params.numAtoms);
        BoxInflate inflater = new BoxInflate(box, space, params.density);
        inflater.actionPerformed();

        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box);

        double rc = 3;
        potentialMaster = new PotentialMasterCell(this, rc, space);

        P2LennardJones p2lj = new P2LennardJones(space);
        P2SoftSphericalTruncated p2 = new P2SoftSphericalTruncated(space, p2lj, rc);
        AtomType atomType = species.getLeafType();
        potentialMaster.addPotential(p2, new AtomType[]{atomType, atomType});
//SINE IS ACTUALLY LN(2+SIN)
         if (params.field==Field.SINE || params.field==Field.PHISINEPSINESUM) {       //AS FOR PHISINEPSINESUM EXTERNAL POT IS LN(2+SIN)-SAME AS FIELD.SINE
            P1Sine p1 = new P1Sine(space, 5, params.temperature);
            potentialMaster.addPotential(p1, new AtomType[]{atomType});
        } else if (params.field == Field.PARABOLIC || params.field == Field.PHIPARABOLICPSUMOFGAUSSIANS || params.field==Field.UNIFORM) {
            P1Parabolic p1 = new P1Parabolic(space);
            potentialMaster.addPotential(p1, new AtomType[]{atomType});
        }

        integrator = new IntegratorMC(this, potentialMaster, box);
        integrator.setTemperature(params.temperature);

        MCMoveAtom mcMoveAtom = new MCMoveAtom(random, potentialMaster, space);
        integrator.getMoveManager().addMCMove(mcMoveAtom);

        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        potentialMaster.setCellRange(2);
        potentialMaster.reset();
        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());
    }

    public static void main(String[] args) {
        SimParams params = new SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
      //      params.field = Field.SINE;
      //          params.field = Field.PARABOLIC;
            params.field = Field.PHIPARABOLICPSUMOFGAUSSIANS;
       //        params.field = Field.PHISINEPSINESUM;
        //    params.field = Field.UNIFORM;
            // modify parameters here for interactive testing
        }

        double aa=params.density;
        double bb=0.0;
        double cc=0.0;
        if(params.temperature==5.0 && params.density==0.125){   bb=0.04835;  cc=0.002934;}
        if(params.temperature==5.0 && params.density==0.25){  bb=0.07508;  cc=0.01126;}
        if(params.temperature==5.0 && params.density==0.375){  bb=0.08929;  cc=0.02323;}
        if(params.temperature==5.0 && params.density==0.5){  bb=0.09321;  cc=0.03947;}
        if(params.temperature==5.0 && params.density==0.625){  bb=0.08904;  cc=0.05789;}
        if(params.temperature==5.0 && params.density==0.75){  bb=0.08278;  cc=0.07673;}
        if(params.temperature==5.0 && params.density==0.875){   bb=0.06952;  cc=0.09625;}
        if(params.temperature==5.0 && params.density==1.0){   bb=0.06052;  cc=0.1058;}

        double a1=0.0;
        double a2=0.0;
        double a3=0.0;
        double a4=0.0;
        double b1=0.0;
        double b2=0.0;
        double b3=0.0;
        double b4=0.0;
        double c1=1.0;
        double c2=1.0;
        double c3=1.0;
        double c4=1.0;
        if(params.temperature==5.0 && params.density==0.125){a1=1.0; a2=-0.4529; b1=0.809; b2=0.8364; c1=0.8498; c2=0.6394;}
        if(params.temperature==5.0 && params.density==0.25){ a1	=	0.2884;a2	=	0.686;b1	=	1.795;b2	=	-0.7539;c1	=	0.6091;c2	=	-1.041;}
        if(params.temperature==5.0 && params.density==0.375){ a1	=	0.8073;a2	=	0.3485;b1	=	0.8895;b2	=	2.079;c1	=	-1.178;c2	=	0.6238;}
        if(params.temperature==5.0 && params.density==0.5){ a1	=	-19.71;a2	=	0.3829;a3	=	-0.006994;a4	=	20.64;b1	=	1.069;b2	=	2.489;b3	=	-0.9805;b4	=	1.069;c1	=	0.902;c2	=	0.5318;c3	=	0.0214;c4	=	0.9197;}
        if(params.temperature==5.0 && params.density==0.625){ a1	=	0.6381;a2	=	0.2972;a3	=	-2.397;a4	=	2.915;b1	=	1.921;b2	=	-2.74;b3	=	0.006159;b4	=	0;c1	=	1.064;c2	=	0.4744;c3	=	2.086;c4	=	1.994;}
        if(params.temperature==5.0 && params.density==0.75){ a1	=	0.5239;a2	=	0.9187;a3	=	0.1986;a4	=	0.8576;b1	=	2.97;b2	=	0.6462;b3	=	1.253;b4	=	2.053;c1	=	0.5098;c2	=	0.9231;c3	=	0.4761;c4	=	0.7447;}
        if(params.temperature==5.0 && params.density==0.875){ a1	=	-0.3673;a2	=	1.452;a3	=	-8.242;a4	=	8.938;b1	=	-0.9993;b2	=	1.055;b3	=	2.766;b4	=	2.77;c1	=	0.5828;c2	=	1.123;c3	=	0.5688;c4	=	0.5953;}
        if(params.temperature==5.0 && params.density==1.0){ a1	=	0.6222;a2	=	1.174;a3	=	0.9962;a4	=	-0.283;b1	=	3.391;b2	=	2.247;b3	=	0.7525;b4	=	2.079;c1	=	0.4002;c2	=	0.9006;c3	=	1.061;c4	=	0.4084;}


        MappedSingletDensity sim = new MappedSingletDensity(params);
        long steps = params.steps;
        int interval = 5* params.numAtoms;
        int blocks = 100;
        long blockSize = steps / (interval * blocks);

        System.out.println("Lennard-Jones Monte Carlo simulation");
        System.out.println("N: " + params.numAtoms);
        System.out.println("T: " + params.temperature);
        System.out.println("density: " + params.density);
        System.out.println("steps: " + params.steps);

        // equilibration
        long t1 = System.currentTimeMillis();
        sim.activityIntegrate.setMaxSteps(steps / 10);
        sim.activityIntegrate.actionPerformed();
        System.out.println("equilibration finished");

        // data collection
        MeterProfileByVolume densityMeter = new MeterProfileByVolume(sim.space);
        densityMeter.setBox(sim.box());
        densityMeter.setProfileDim(2);
        densityMeter.getXDataSource().setNValues(50);  //conventional bins=100
        densityMeter.reset();

        MeterNMolecules meterNMolecules = new MeterNMolecules();
        densityMeter.setDataSource(meterNMolecules);
        AccumulatorAverageFixed acc = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pump = new DataPumpListener(densityMeter, acc, interval);
        sim.getIntegrator().getEventManager().addListener(pump);

        MeterProfileForceSum densityMeterForce = new MeterProfileForceSum(sim.box(), sim.potentialMaster, params.temperature);
        densityMeterForce.setProfileDim(2);
        densityMeterForce.getXDataSource().setNValues(50);  //pap bins=1000
        densityMeterForce.reset();
        AccumulatorAverageFixed accForce = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pumpForce = new DataPumpListener(densityMeterForce, accForce, interval);
        sim.getIntegrator().getEventManager().addListener(pumpForce);

        FunctionDifferentiable f;
        double L = sim.box().getBoundary().getBoxSize().getX(2);
        switch (params.field) {
            case PARABOLIC:
                f = new FunctionParabolic(L, params.temperature);
                break;
            case SINE:
                f = new FunctionSine(5, L);
                break;
            case UNIFORM:
                f = new FunctionUniform(L);
                break;
            case PHIPARABOLICPSUMOFGAUSSIANS:
                f = new FunctionPhiparabolicpsumofgaussians(L,a1,b1,c1,a2,b2,c2,a3,b3,c3,a4,b4,c4);
                break;
            case PHISINEPSINESUM:
                f = new FunctionPhisinepsinesum(L,aa,bb,cc);
                break;
            default:
                throw new RuntimeException("not yet");
        }
        MeterProfileMappedAvg densityMeterMappedAvg = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
        densityMeterMappedAvg.setProfileDim(2);
        densityMeterMappedAvg.getXDataSource().setNValues(50);  //mappedavg bins=1000
        densityMeterMappedAvg.reset();
        AccumulatorAverageFixed accMappedAvg = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pumpMappedAvg = new DataPumpListener(densityMeterMappedAvg, accMappedAvg, 5 * params.numAtoms);
        sim.getIntegrator().getEventManager().addListener(pumpMappedAvg);

        sim.activityIntegrate.setMaxSteps(steps);
        sim.getIntegrator().resetStepCount();
        sim.integrator.getMoveManager().setEquilibrating(false);
        sim.activityIntegrate.actionPerformed();

        long t2 = System.currentTimeMillis();
        IData data =  acc.getData(acc.AVERAGE);
        IData dataunc =  acc.getData(acc.ERROR);
        IData dataForce =  accForce.getData(accForce.AVERAGE);
        IData dataForceunc =  accForce.getData(accForce.ERROR);
        IData dataMappedAvg =  accMappedAvg.getData(accMappedAvg.AVERAGE);
        IData dataMappedAvgunc =  accMappedAvg.getData(accMappedAvg.ERROR);

       IData zdata= ((DataFunction.DataInfoFunction)((DataGroup.DataInfoGroup)acc.getDataInfo()).getSubDataInfo(0)).getXDataSource().getIndependentData(0);
       for (int i=0;i<zdata.getLength();i++){
           System.out.println(zdata.getValue(i)+" "+data.getValue(i)+" "+dataForce.getValue(i)+" "+dataMappedAvg.getValue(i)+" "+dataunc.getValue(i)+" "+dataForceunc.getValue(i)+" "+dataMappedAvgunc.getValue(i));
       }
        System.out.println("time: " + (t2 - t1) * 0.001);
    }

    public static class Graphic {
        public static void main(String[] args) {
            SimParams params = new SimParams();

            if (args.length > 0) {
                ParseArgs.doParseArgs(params, args);
            }
            else {
            //    params.field = Field.SINE;
           //         params.field = Field.PARABOLIC;
                params.field = Field.PHIPARABOLICPSUMOFGAUSSIANS;
           //      params.field = Field.PHISINEPSINESUM;
             //   params.field = Field.UNIFORM;
                // modify parameters here for interactive testing
            }


            double aa=params.density;
            double bb=0.0;
            double cc=0.0;
            if(params.temperature==5.0 && params.density==0.125){   bb=0.04835;  cc=0.002934;}
            if(params.temperature==5.0 && params.density==0.25){  bb=0.07508;  cc=0.01126;}
            if(params.temperature==5.0 && params.density==0.375){  bb=0.08929;  cc=0.02323;}
            if(params.temperature==5.0 && params.density==0.5){  bb=0.09321;  cc=0.03947;}
            if(params.temperature==5.0 && params.density==0.625){  bb=0.08904;  cc=0.05789;}
            if(params.temperature==5.0 && params.density==0.75){  bb=0.08278;  cc=0.07673;}
            if(params.temperature==5.0 && params.density==0.875){   bb=0.06952;  cc=0.09625;}
            if(params.temperature==5.0 && params.density==1.0){   bb=0.06052;  cc=0.1058;}

            double a1=0.0;
            double a2=0.0;
            double a3=0.0;
            double a4=0.0;
            double b1=0.0;
            double b2=0.0;
            double b3=0.0;
            double b4=0.0;
            double c1=1.0;
            double c2=1.0;
            double c3=1.0;
            double c4=1.0;
            if(params.temperature==5.0 && params.density==0.125){a1=1.0; a2=-0.4529; b1=0.809; b2=0.8364; c1=0.8498; c2=0.6394;}
            if(params.temperature==5.0 && params.density==0.25){ a1	=	0.2884;a2	=	0.686;b1	=	1.795;b2	=	-0.7539;c1	=	0.6091;c2	=	-1.041;}
            if(params.temperature==5.0 && params.density==0.375){ a1	=	0.8073;a2	=	0.3485;b1	=	0.8895;b2	=	2.079;c1	=	-1.178;c2	=	0.6238;}
            if(params.temperature==5.0 && params.density==0.5){ a1	=	-19.71;a2	=	0.3829;a3	=	-0.006994;a4	=	20.64;b1	=	1.069;b2	=	2.489;b3	=	-0.9805;b4	=	1.069;c1	=	0.902;c2	=	0.5318;c3	=	0.0214;c4	=	0.9197;}
            if(params.temperature==5.0 && params.density==0.625){ a1	=	0.6381;a2	=	0.2972;a3	=	-2.397;a4	=	2.915;b1	=	1.921;b2	=	-2.74;b3	=	0.006159;b4	=	0;c1	=	1.064;c2	=	0.4744;c3	=	2.086;c4	=	1.994;}
            if(params.temperature==5.0 && params.density==0.75){ a1	=	0.5239;a2	=	0.9187;a3	=	0.1986;a4	=	0.8576;b1	=	2.97;b2	=	0.6462;b3	=	1.253;b4	=	2.053;c1	=	0.5098;c2	=	0.9231;c3	=	0.4761;c4	=	0.7447;}
            if(params.temperature==5.0 && params.density==0.875){ a1	=	-0.3673;a2	=	1.452;a3	=	-8.242;a4	=	8.938;b1	=	-0.9993;b2	=	1.055;b3	=	2.766;b4	=	2.77;c1	=	0.5828;c2	=	1.123;c3	=	0.5688;c4	=	0.5953;}
            if(params.temperature==5.0 && params.density==1.0){ a1	=	0.6222;a2	=	1.174;a3	=	0.9962;a4	=	-0.283;b1	=	3.391;b2	=	2.247;b3	=	0.7525;b4	=	2.079;c1	=	0.4002;c2	=	0.9006;c3	=	1.061;c4	=	0.4084;}


            MappedSingletDensity sim = new MappedSingletDensity(params);
            SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);

            int blockSize = 100;
            MeterProfileByVolume densityMeter = new MeterProfileByVolume(sim.space);
            densityMeter.setBox(sim.box());
            densityMeter.setProfileDim(2);
            MeterNMolecules meterNMolecules = new MeterNMolecules();
            densityMeter.setDataSource(meterNMolecules);
            AccumulatorAverageFixed acc = new AccumulatorAverageFixed(5 * blockSize);
            DataPumpListener pump = new DataPumpListener(densityMeter, acc, params.numAtoms);
            sim.getIntegrator().getEventManager().addListener(pump);

            MeterProfileForceSum densityMeterForce = new MeterProfileForceSum(sim.box(), sim.potentialMaster, params.temperature);
            densityMeterForce.setProfileDim(2);
            AccumulatorAverageFixed accForce = new AccumulatorAverageFixed(blockSize);
            DataPumpListener pumpForce = new DataPumpListener(densityMeterForce, accForce, 5 * params.numAtoms);
            sim.getIntegrator().getEventManager().addListener(pumpForce);

            FunctionDifferentiable f;
            double L = sim.box().getBoundary().getBoxSize().getX(2);
            switch (params.field) {
                case PARABOLIC:
                    f = new FunctionParabolic(L, params.temperature);
                    break;
                case SINE:
                    f = new FunctionSine(5, L);
                    break;
                case UNIFORM:
                    f = new FunctionUniform(L);
                    break;
                case PHISINEPSINESUM:
                    f = new FunctionPhisinepsinesum(L,aa,bb,cc);
                    break;
                case PHIPARABOLICPSUMOFGAUSSIANS:
                    f = new FunctionPhiparabolicpsumofgaussians(L,a1,b1,c1,a2,b2,c2,a3,b3,c3,a4,b4,c4);
                    break;
                default:
                    throw new RuntimeException("not yet");
            }
            MeterProfileMappedAvg densityMeterMappedAvg = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterMappedAvg.setProfileDim(2);
            AccumulatorAverageFixed accMappedAvg = new AccumulatorAverageFixed(blockSize);
            DataPumpListener pumpMappedAvg = new DataPumpListener(densityMeterMappedAvg, accMappedAvg, 5 * params.numAtoms);
            sim.getIntegrator().getEventManager().addListener(pumpMappedAvg);

            MeterProfileMappedAvg densityMeterP = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterP.getXDataSource().setNValues(104);
            densityMeterP.setProfileDim(2);
            densityMeterP.setBehavior(MeterProfileMappedAvg.Behavior.P);
            MeterProfileMappedAvg densityMeterZidot0 = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterZidot0.getXDataSource().setNValues(104);
            densityMeterZidot0.setProfileDim(2);
            densityMeterZidot0.setBehavior(MeterProfileMappedAvg.Behavior.ZIDOT);
            densityMeterZidot0.setZidotZ(0);
            MeterProfileMappedAvg densityMeterZidot1 = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterZidot1.getXDataSource().setNValues(104);
            densityMeterZidot1.setProfileDim(2);
            densityMeterZidot1.setBehavior(MeterProfileMappedAvg.Behavior.ZIDOT);
            densityMeterZidot1.setZidotZ(L / 8);
            MeterProfileMappedAvg densityMeterZidot2 = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterZidot2.getXDataSource().setNValues(104);
            densityMeterZidot2.setProfileDim(2);
            densityMeterZidot2.setBehavior(MeterProfileMappedAvg.Behavior.ZIDOT);
            densityMeterZidot2.setZidotZ(L / 4);

            MeterProfileMappedAvg densityMeterDZidot0 = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterDZidot0.getXDataSource().setNValues(104);
            densityMeterDZidot0.setProfileDim(2);
            densityMeterDZidot0.setBehavior(MeterProfileMappedAvg.Behavior.DZIDOT);
            densityMeterDZidot0.setZidotZ(0);
            MeterProfileMappedAvg densityMeterDZidot1 = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterDZidot1.getXDataSource().setNValues(104);
            densityMeterDZidot1.setProfileDim(2);
            densityMeterDZidot1.setBehavior(MeterProfileMappedAvg.Behavior.DZIDOT);
            densityMeterDZidot1.setZidotZ(L / 8);
            MeterProfileMappedAvg densityMeterDZidot2 = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterDZidot2.getXDataSource().setNValues(104);
            densityMeterDZidot2.setProfileDim(2);
            densityMeterDZidot2.setBehavior(MeterProfileMappedAvg.Behavior.DZIDOT);
            densityMeterDZidot2.setZidotZ(L / 4);

            DisplayPlot densityPlot = new DisplayPlot();
            acc.addDataSink(densityPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.AVERAGE});
            accForce.addDataSink(densityPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.AVERAGE});
            accMappedAvg.addDataSink(densityPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.AVERAGE});

            new DataPump(densityMeterP, densityPlot.getDataSet().makeDataSink()).actionPerformed();

            densityPlot.setLabel("density");
            graphic.add(densityPlot);

            DisplayPlot errorPlot = new DisplayPlot();
            acc.addDataSink(errorPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.ERROR});
            accForce.addDataSink(errorPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.ERROR});
            accMappedAvg.addDataSink(errorPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.ERROR});
            errorPlot.setLabel("error");
            graphic.add(errorPlot);

            DisplayPlot corPlot = new DisplayPlot();
            acc.addDataSink(corPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.BLOCK_CORRELATION});
            accForce.addDataSink(corPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.BLOCK_CORRELATION});
            accMappedAvg.addDataSink(corPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.BLOCK_CORRELATION});
            corPlot.setLabel("correlation");
            graphic.add(corPlot);

            DisplayPlot zidotPlot = new DisplayPlot();
            new DataPump(densityMeterZidot0, zidotPlot.getDataSet().makeDataSink()).actionPerformed();
            new DataPump(densityMeterZidot1, zidotPlot.getDataSet().makeDataSink()).actionPerformed();
            new DataPump(densityMeterZidot2, zidotPlot.getDataSet().makeDataSink()).actionPerformed();
            zidotPlot.setLabel("zidot");
            graphic.add(zidotPlot);

            DisplayPlot dzidotPlot = new DisplayPlot();
            new DataPump(densityMeterDZidot0, dzidotPlot.getDataSet().makeDataSink()).actionPerformed();
            new DataPump(densityMeterDZidot1, dzidotPlot.getDataSet().makeDataSink()).actionPerformed();
            new DataPump(densityMeterDZidot2, dzidotPlot.getDataSet().makeDataSink()).actionPerformed();
            dzidotPlot.setLabel("-dzidot");
            graphic.add(dzidotPlot);

            List<DataPump> pumps = graphic.getController().getDataStreamPumps();
            pumps.add(pump);
            pumps.add(pumpForce);
            pumps.add(pumpMappedAvg);

            graphic.makeAndDisplayFrame();

        }
    }

    enum Field {
        SINE, UNIFORM, PARABOLIC, PHISINEPSINESUM, PHIPARABOLICPSUMOFGAUSSIANS
    }

    public static class SimParams extends ParameterBase {
        public long steps = 250000;
        public double density = 1.0;
        public double temperature = 5.0;
        public int numAtoms = 500;
        public Field field = Field.PHIPARABOLICPSUMOFGAUSSIANS;
    }
}
