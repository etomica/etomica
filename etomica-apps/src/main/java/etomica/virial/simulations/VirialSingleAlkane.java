package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.meter.MeterGyrationTensor;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graph.model.Graph;
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graphics.*;
import etomica.integrator.*;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveAtomFasterer;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.parser.ParserAMBER;
import etomica.potential.IPotential;
import etomica.potential.Potential2Soft;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.PotentialMoleculePair;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;
import etomica.units.*;
import etomica.units.dimensions.*;
import etomica.units.dimensions.Dimension;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.collections.IntArrayList;
import etomica.virial.MayerHardSphere;
import etomica.virial.PotentialComputeIntramolecular;
import etomica.virial.mcmove.*;

import javax.swing.*;
import java.awt.*;

/**
 * Shape descriptors calculation for alkanes using the AMBER force field.
 */

public class VirialSingleAlkane {
    public static void main(String[] args) {
        VirialSingleAlkaneParam params = new VirialSingleAlkaneParam();
        boolean isCommandline = args.length > 0;
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.temperature = 400.0;// Kelvin
            params.numSteps = 1000000;
            params.file = "CCCCC";
        }

        double temperature = params.temperature;
        long steps = params.numSteps;
        String file = params.file;
        ParserAMBER.Stuff stuff= ParserAMBER.makeStuff(file, new ParserAMBER.Options());
        SpeciesManager speciesManager = stuff.speciesManager;
        ISpecies alkane = speciesManager.getSpecies(0);
        AtomType C = alkane.getAtomType(0);
        AtomType H = alkane.getAtomType(1);
        PotentialMasterBonding.FullBondingInfo fullBondingInfo = stuff.fullBondingInfo;
        Potential2Soft[][] potential2Soft = stuff.potential2Soft;

        IntArrayList[] bonding = new IntArrayList[speciesManager.getSpecies(0).getLeafAtomCount()];
        for(int atomIndex = 0; atomIndex < bonding.length; atomIndex++){
            bonding[atomIndex] = new IntArrayList();
        }

        for(int[][][] bondType : fullBondingInfo.bondedPairPartners[0].values()){
            for(int atomIndex = 0; atomIndex < bondType.length; atomIndex++){
                for(int[] bonds: bondType[atomIndex]){
                    int second = bonds[0] + bonds[1] - atomIndex;
                    bonding[atomIndex].add(second);
                }
            }
        }

        Space space = Space3D.getInstance();

        PotentialMoleculePair potentialMoleculePair = new PotentialMoleculePair(space, speciesManager);
        potentialMoleculePair.setAtomPotential(C, C, potential2Soft[0][0]);
        potentialMoleculePair.setAtomPotential(C, H, potential2Soft[0][1]);
        potentialMoleculePair.setAtomPotential(H, H, potential2Soft[1][1]);

        System.out.println("numSteps " + steps);
        System.out.println("Shape Descriptors at " + temperature +"K");
        temperature = Kelvin.UNIT.toSim(temperature);

        long t1 = System.currentTimeMillis();

        final Simulation sim = new Simulation(space, speciesManager);
        Box box = new Box(new BoundaryRectangularNonperiodic(space), space);
        sim.addBox(box);
        box.setNMolecules(alkane, 1);

        PotentialComputeIntramolecular potentialComputeIntramolecular0 = new PotentialComputeIntramolecular(space, sim.box(), speciesManager, fullBondingInfo, potential2Soft);
        PotentialMasterBonding potentialMasterBonding0 = new PotentialMasterBonding(speciesManager, sim.box(), fullBondingInfo);
        PotentialComputeAggregate potentialComputeAggregate0 = new PotentialComputeAggregate(potentialComputeIntramolecular0, potentialMasterBonding0);

        MeterGyrationTensor meterGyrationTensor = new MeterGyrationTensor(space);
        meterGyrationTensor.setBox(sim.box());

        IntegratorMCFasterer integrator = new IntegratorMCFasterer(potentialComputeAggregate0, sim.getRandom(), temperature, sim.box());

        MCMoveAtomFasterer mcMoveAtomFasterer = new MCMoveAtomFasterer(sim.getRandom(), potentialComputeAggregate0, sim.box());
        integrator.getMoveManager().addMCMove(mcMoveAtomFasterer);

        MCMoveAngle mcMoveAngle0 = new MCMoveAngle(potentialComputeAggregate0, space, bonding, sim.getRandom(), 1);
        integrator.getMoveManager().addMCMove(mcMoveAngle0);

        MCMoveTorsion mcMoveTorsion0 = new MCMoveTorsion(potentialComputeAggregate0, space,bonding, sim.getRandom(), 2*Math.PI);
        if(fullBondingInfo.bondedQuads[0].size() > 0) {
            integrator.getMoveManager().addMCMove(mcMoveTorsion0);
        }

        if(false) {
            double size = (bonding.length/3 + 5) * 1.5;
            sim.getController().addActivity(new ActivityIntegrate(integrator));
            sim.box().getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box());
            displayBox0.setPixelUnit(new Pixel(300.0 / size));
            displayBox0.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);

            DiameterHashByType diameterManager = (DiameterHashByType) displayBox0.getDiameterHash();
            diameterManager.setDiameter(speciesManager.getSpecies(0).getAtomType(0), 1);
            diameterManager.setDiameter(speciesManager.getSpecies(0).getAtomType(1), 0.5);
            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim.getSpeciesManager(), sim.box(), sim.getRandom());
            displayBox0.setColorScheme(colorScheme);
            simGraphic.makeAndDisplayFrame();
            return;
        }

        sim.getController().runActivityBlocking(new ActivityIntegrate(integrator, steps/10));

        ActivityIntegrate ai = new ActivityIntegrate(integrator, steps);

        System.out.println("equilibration finished");
        integrator.getMoveManager().setEquilibrating(false);
        AccumulatorAverageFixed accumulator = new AccumulatorAverageFixed((steps/100)/speciesManager.getSpecies(0).getLeafAtomCount());
        DataPumpListener pump = new DataPumpListener(meterGyrationTensor, accumulator, speciesManager.getSpecies(0).getLeafAtomCount());
        integrator.getEventManager().addListener(pump);
        sim.getController().runActivityBlocking(ai);

        System.out.println("Acceptance of Torsion Move "+ mcMoveTorsion0.getTracker().acceptanceProbability());

        double RGAvg = accumulator.getData(AccumulatorAverage.AVERAGE).getValue(0);
        double RGErr = accumulator.getData(AccumulatorAverage.ERROR).getValue(0);
        double RGCor = accumulator.getData(AccumulatorAverage.BLOCK_CORRELATION).getValue(0);
        double RGSD = accumulator.getData(AccumulatorAverage.STANDARD_DEVIATION).getValue(0);

        double AsphericityAvg = accumulator.getData(AccumulatorAverage.AVERAGE).getValue(1);
        double AsphericityErr = accumulator.getData(AccumulatorAverage.ERROR).getValue(1);
        double AsphericityCor = accumulator.getData(AccumulatorAverage.BLOCK_CORRELATION).getValue(1);
        double AsphericitySD = accumulator.getData(AccumulatorAverage.STANDARD_DEVIATION).getValue(1);

        double AcylindricityAvg = accumulator.getData(AccumulatorAverage.AVERAGE).getValue(2);
        double AcylindricityErr = accumulator.getData(AccumulatorAverage.ERROR).getValue(2);
        double AcylindricityCor = accumulator.getData(AccumulatorAverage.BLOCK_CORRELATION).getValue(2);
        double AcylindricitySD = accumulator.getData(AccumulatorAverage.STANDARD_DEVIATION).getValue(2);

        double AnisotropyAvg = accumulator.getData(AccumulatorAverage.AVERAGE).getValue(3);
        double AnisotropyErr = accumulator.getData(AccumulatorAverage.ERROR).getValue(3);
        double AnisotropyCor = accumulator.getData(AccumulatorAverage.BLOCK_CORRELATION).getValue(3);
        double AnisotropySD = accumulator.getData(AccumulatorAverage.STANDARD_DEVIATION).getValue(3);

        System.out.print(String.format("RadiusOfGyration avg: %13.6e  sd: %13.6e  err: %11.4e   cor: % 4.2f\n", RGAvg, RGSD, RGErr, RGCor));
        System.out.print(String.format("Asphericity avg: %13.6e  sd: %13.6e  err: %11.4e   cor: % 4.2f\n", AsphericityAvg, AsphericitySD, AsphericityErr, AsphericityCor));
        System.out.print(String.format("Acylindricity avg: %13.6e  sd: %13.6e  err: %11.4e   cor: % 4.2f\n", AcylindricityAvg, AcylindricitySD, AcylindricityErr, AcylindricityCor));
        System.out.print(String.format("Anisotropy avg: %13.6e  sd: %13.6e  err: %11.4e   cor: % 4.2f\n", AnisotropyAvg, AnisotropySD, AnisotropyErr, AnisotropyCor));

        long t2 = System.currentTimeMillis();
        System.out.println("time: "+(t2-t1)*0.001);
    }

    /**
     * Inner class for parameters
     */
    public static class VirialSingleAlkaneParam extends ParameterBase {
        public double temperature = 300.0;// Kelvin
        public long numSteps = 100000000;
        public String file = "CCCCCC";
    }
}
