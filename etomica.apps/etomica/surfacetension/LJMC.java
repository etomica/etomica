/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.surfacetension;
import java.util.List;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataDump;
import etomica.data.DataProcessor;
import etomica.data.DataPump;
import etomica.data.DataPumpListener;
import etomica.data.DataTag;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterPressureTensor;
import etomica.data.meter.MeterProfileByVolume;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.listener.IntegratorListenerAction;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

public class LJMC extends Simulation {
    
    public final PotentialMasterCell potentialMaster;
    public final SpeciesSpheresMono species;
    public final IBox box;
    public final ActivityIntegrate activityIntegrate;
    public final IntegratorMC integrator;
    public final MCMoveAtom mcMoveAtom, mcMoveAtomBigStep;

    public LJMC(ISpace _space, int numAtoms, double temperature, double aspectRatio) {
        super(_space);
        double rc = 4.0;
        potentialMaster = new PotentialMasterCell(this, rc, space);
        potentialMaster.setCellRange(2);
        potentialMaster.lrcMaster().setEnabled(false);
        
        //controller and integrator
	    integrator = new IntegratorMC(potentialMaster, random, 1.0);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

	    //species and potentials
	    species = new SpeciesSpheresMono(this, space);//index 1
	    species.setIsDynamic(true);
        addSpecies(species);
        
        //instantiate several potentials for selection in combo-box
	    P2LennardJones potential = new P2LennardJones(space);
        P2SoftSphericalTruncated p2Truncated = new P2SoftSphericalTruncated(space, potential, rc);
	    potentialMaster.addPotential(p2Truncated, new IAtomType[]{species.getLeafType(), species.getLeafType()});
	    
        //construct box
	    box = new Box(space);
        addBox(box);
        IVectorMutable dim = space.makeVector();
        box.getBoundary().setBoxSize(dim);
        box.setNMolecules(species, numAtoms);
        integrator.setBox(box);
        potentialMaster.getNbrCellManager(box).assignCellAll();
        
        mcMoveAtom = new MCMoveAtom(random, potentialMaster, space);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        mcMoveAtomBigStep = new MCMoveAtom(random, potentialMaster, space);
        mcMoveAtomBigStep.setStepSize(6);
        mcMoveAtomBigStep.setStepSizeMin(5);
        integrator.getMoveManager().addMCMove(mcMoveAtomBigStep);
        integrator.getMoveManager().setFrequency(mcMoveAtomBigStep, 0.2);
        
        //new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        double Tc = 1.31;
        double rhoc = 0.314;
        double T = temperature;
        double rhoV = rhoc - 0.477*Math.pow(Tc-T,1.0/3.0) + 0.05333*(Tc-T) + 0.1261*Math.pow(Tc-T, 1.5);
        double rhoL = rhoc + 0.477*Math.pow(Tc-T,1.0/3.0) + 0.2124*(Tc-T) - 0.01151*Math.pow(Tc-T, 1.5);
        double xL = 1.0/(1.0 + Math.pow(rhoL/rhoV,1.0/3.0));
        // V = NV/rhoV + NL/rhoL
        // rho = 0.75 rhoV + 0.25 rhoL
        double rho = (1-xL)*rhoV + xL*rhoL;
        double V = numAtoms / rho;

        // Lx * Lyz * Lyz = V
        // Lx = Lyz * 4
        // 4 * Lyz^3 = V
        // Lyz = (V/4)^(1/3)
        
        double Lyz = Math.pow(V/aspectRatio, 1.0/3.0);
        if (Lyz < 2*rc) {
            throw new RuntimeException("cutoff("+rc+") too large for box size ("+Lyz+")");
        }
        double Lx = Lyz*aspectRatio;
        
        //initialize box to just liquid size, put atoms on FCC lattice
        box.getBoundary().setBoxSize(space.makeVector(new double[]{Lx*xL,Lyz,Lyz}));
        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        // then expand box size
        
        box.getBoundary().setBoxSize(space.makeVector(new double[]{Lx,Lyz,Lyz}));

        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());
    }
    
    public static void main(String[] args) {
        LJMCParams params = new LJMCParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
        }
        ISpace space = Space3D.getInstance();
        int numAtoms = params.numAtoms;
        double temperature = params.temperature;
        double aspectRatio = params.aspectRatio;
        

        final LJMC sim = new LJMC(space, numAtoms, temperature, aspectRatio);
        
        MeterPressureTensor meterPTensor = new MeterPressureTensor(sim.potentialMaster, space);
        meterPTensor.setBox(sim.box);
        meterPTensor.setTemperature(temperature);
        
        if (true) {
            SimulationGraphic ljmcGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
            SimulationGraphic.makeAndDisplayFrame(ljmcGraphic.getPanel(), "LJ surface tension");
            
            IAction recenterAction = new ActionRecenter(sim);
            IntegratorListenerAction recenterActionListener = new IntegratorListenerAction(recenterAction);
            sim.integrator.getEventManager().addListener(recenterActionListener);
            recenterActionListener.setInterval(100);
            
            
            List<DataPump> dataStreamPumps = ljmcGraphic.getController().getDataStreamPumps();
            
            MeterProfileByVolume densityProfileMeter = new MeterProfileByVolume(space);
            densityProfileMeter.setBox(sim.box);
            MeterNMolecules meterNMolecules = new MeterNMolecules();
            meterNMolecules.setSpecies(sim.species);
            densityProfileMeter.setDataSource(meterNMolecules);
            AccumulatorAverageFixed densityProfileAvg = new AccumulatorAverageFixed(10);
            densityProfileAvg.setPushInterval(10);
            DataPump profilePump = new DataPump(densityProfileMeter, densityProfileAvg);
            DataDump profileDump = new DataDump();
            densityProfileAvg.addDataSink(profileDump, new AccumulatorAverage.StatType[]{densityProfileAvg.AVERAGE});
            IntegratorListenerAction profilePumpListener = new IntegratorListenerAction(profilePump);
            sim.integrator.getEventManager().addListener(profilePumpListener);
            profilePumpListener.setInterval(10);
            dataStreamPumps.add(profilePump);

            final FitTanh fitTanh = new FitTanh();
            densityProfileAvg.addDataSink(fitTanh, new AccumulatorAverage.StatType[]{densityProfileAvg.AVERAGE});

            DisplayPlot profilePlot = new DisplayPlot();
            densityProfileAvg.addDataSink(profilePlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{densityProfileAvg.AVERAGE});
            fitTanh.setDataSink(profilePlot.getDataSet().makeDataSink());
            profilePlot.setLegend(new DataTag[]{densityProfileAvg.getTag()}, "density");
            profilePlot.setLegend(new DataTag[]{densityProfileAvg.getTag(), fitTanh.getTag()}, "fit");
            profilePlot.setDoLegend(true);
            profilePlot.setLabel("Density");
            ljmcGraphic.add(profilePlot);
            
            DataProcessor surfaceTension = new DataProcessorSurfaceTension();
            
            DataPumpListener tensionPump = new DataPumpListener(meterPTensor, surfaceTension, numAtoms);
            AccumulatorAverageCollapsing surfaceTensionAvg = new AccumulatorAverageCollapsing();
            surfaceTension.setDataSink(surfaceTensionAvg);
            DisplayTextBoxesCAE tensionBox = new DisplayTextBoxesCAE();
            tensionBox.setAccumulator(surfaceTensionAvg);
            tensionBox.setLabel("Surface Tension");
            tensionBox.setPrecision(4);
            dataStreamPumps.add(tensionPump);
            sim.integrator.getEventManager().addListener(tensionPump);
            ljmcGraphic.add(tensionBox);
            
            
            
            
            return;
        }
        
        sim.getController().actionPerformed();

    }
    
    public static class LJMCParams extends ParameterBase {
        public int numAtoms = 1000;
        public double temperature = 1.1;
        public double aspectRatio = 6;
    }
}
