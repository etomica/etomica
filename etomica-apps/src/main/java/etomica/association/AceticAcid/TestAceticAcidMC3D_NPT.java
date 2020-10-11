/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association.AceticAcid;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;

import etomica.action.activity.ActivityIntegrate;
import etomica.association.*;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.atom.iterator.Atomset4IteratorIndexList;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.integrator.mcmove.MCMoveTorsionAceticAcid;
import etomica.integrator.mcmove.MCMoveWiggleAceticAcid;
import etomica.lattice.LatticeCubicFcc;
import etomica.models.OPLS.AceticAcidPotentialHelper;
import etomica.models.OPLS.DipoleSourceAceticAcid;
import etomica.models.OPLS.SpeciesAceticAcid;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.nbr.cell.molecule.BoxAgentSourceCellManagerMolecular;
import etomica.nbr.cell.molecule.NeighborCellManagerMolecular;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesGeneral;
import etomica.units.*;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;
import java.io.File;

/**
 * OPLS Acetic Acid Monte Carlo NPT simulation in 3D.
 * average density = N*<1/V>
 * Translation/Rotation/Torsion/Wiggle/Volume change + BiasUB
 * 
 * @author Hye Min Kim
 */
public class TestAceticAcidMC3D_NPT extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    public MCMoveMoleculeMonomer mcMoveMoleculeMonomer;
    public MCMoveMoleculeSmer mcMoveMoleculeSmer;
    public MCMoveMoleculeRotateAssociated mcMoveMoleculeRotate;
    public MCMoveVolumeAssociatedMolecule mcMoveVolume;
    public MCMoveTorsionAceticAcid mcMoveTorsion;
    public MCMoveWiggleAceticAcid mcMoveWiggle;
    public MCMoveBiasUBMolecule mcMoveBiasUB;
    public SpeciesGeneral species;
    public Box box;
    public PotentialGroup potential;
    public P2ReactionFieldDipole reactionField;

    public AssociationManagerMolecule associationManager;
    public AssociationHelperMolecule associationHelper;
    public BiasVolumeAceticAcid bv;
        
    public TestAceticAcidMC3D_NPT(int numAtoms, double pressureBar, double densityMolLiter, double temperatureK, long numSteps) {
        super(Space3D.getInstance());

        species = SpeciesAceticAcid.create();
        addSpecies(species);

        PotentialMaster potentialMaster = new PotentialMaster();
        //setRandom(new RandomNumberGenerator(3));
        IMoleculePositionDefinition positionDefinition = new IMoleculePositionDefinition() {//anonymous class
            public Vector position(IMolecule molecule) {
                return molecule.getChildList().get(SpeciesAceticAcid.indexC).getPosition();
            }
        };
        BoxAgentSourceCellManagerMolecular bASCellManagerMolecular = new BoxAgentSourceCellManagerMolecular(this, positionDefinition, space);//tracking neighbors
        bASCellManagerMolecular.setRange(4.2);//association is made within 4.2A of C-C distance
        BoxAgentManager<NeighborCellManagerMolecular> cellAgentManager = new BoxAgentManager<NeighborCellManagerMolecular>(bASCellManagerMolecular, this);
        System.out.println("pressure = " + pressureBar + "bar");
        System.out.println("initial density = " + densityMolLiter + "mol/L");
        System.out.println("temperature = " + temperatureK + "K");
        System.out.println("numSteps = " + numSteps + "steps");
        CompoundUnit rhoUnit = new CompoundUnit(new Unit[]{Mole.UNIT, Liter.UNIT}, new double[]{1, -1});
        double density = rhoUnit.toSim(densityMolLiter);
        double pressure = Bar.UNIT.toSim(pressureBar);
        double temperature = Kelvin.UNIT.toSim(temperatureK);
        double volume = 1 / (density / numAtoms);
        double boxLength = Math.pow(volume, 1.0 / 3.0);

        box = this.makeBox();
        integrator = new IntegratorMC(this, potentialMaster, box);
        integrator.setTemperature(temperature);
        mcMoveMoleculeMonomer = new MCMoveMoleculeMonomer(this, potentialMaster, space);//Standard Monte Carlo atom-displacement trial move
        mcMoveMoleculeSmer = new MCMoveMoleculeSmer(this, potentialMaster, space);
        mcMoveMoleculeRotate = new MCMoveMoleculeRotateAssociated(potentialMaster, random, space);//Performs a rotation of an atom (not a molecule) that has an orientation coordinate;
        bv = new BiasVolumeAceticAcid(space, random, box);
        bv.setBox(box);

        associationManager = new AssociationManagerMolecule(this, box, cellAgentManager, bv, 4.2);//check the association within the range 4.2A, maximum distance of association = 4.2A
        associationHelper = new AssociationHelperMolecule(space, box, associationManager);
        mcMoveBiasUB = new MCMoveBiasUBMolecule(potentialMaster, bv, random, space);
        mcMoveMoleculeMonomer.setAssociationManager(associationManager, associationHelper);
        mcMoveMoleculeSmer.setAssociationManager(associationManager, associationHelper);
        mcMoveMoleculeRotate.setAssociationManager(associationManager, associationHelper);
        mcMoveBiasUB.setAssociationManager(associationManager, associationHelper);
        mcMoveVolume = new MCMoveVolumeAssociatedMolecule(this, potentialMaster, space);//volume change move
        mcMoveVolume.setPressure(pressure);
        mcMoveVolume.setAssociationManager(associationManager, associationHelper);

        mcMoveTorsion = new MCMoveTorsionAceticAcid(potentialMaster, space, random);
        mcMoveWiggle = new MCMoveWiggleAceticAcid(potentialMaster, random, 0.1, space);

        integrator.getMoveEventManager().addListener(associationManager);

        ((MCMoveStepTracker) mcMoveMoleculeMonomer.getTracker()).setNoisyAdjustment(true);
        ((MCMoveStepTracker) mcMoveMoleculeSmer.getTracker()).setNoisyAdjustment(true);
        ((MCMoveStepTracker) mcMoveMoleculeRotate.getTracker()).setNoisyAdjustment(true);
        ((MCMoveStepTracker) mcMoveVolume.getTracker()).setNoisyAdjustment(true);
        ((MCMoveStepTracker) mcMoveWiggle.getTracker()).setNoisyAdjustment(true);
        integrator.getMoveManager().addMCMove(mcMoveMoleculeMonomer);
        integrator.getMoveManager().addMCMove(mcMoveMoleculeSmer);
        integrator.getMoveManager().addMCMove(mcMoveMoleculeRotate);
        integrator.getMoveManager().addMCMove(mcMoveTorsion);
        integrator.getMoveManager().addMCMove(mcMoveWiggle);
        integrator.getMoveManager().addMCMove(mcMoveVolume);
        integrator.getMoveManager().addMCMove(mcMoveBiasUB);
        integrator.getMoveManager().setEquilibrating(true);
        this.getController().addActivity(new ActivityIntegrate(integrator), numSteps);
        //actionIntegrate.setSleepPeriod(1);
        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space);//Performs actions that cause volume of system to expand or contract
        inflater.setTargetDensity(density);
        inflater.actionPerformed();

        potential = new PotentialGroup(2);
        potential.setBox(box);
        System.out.println("cut-off Distance " + boxLength * 0.49 + " A");
        AceticAcidPotentialHelper.initPotential(space, species, potential);

        double thetaEqCCDBO = 126 * Math.PI / 180;//equilibrium bond angle [=] degree, C-C=O
        double thetaEqDBOCSBO = 123 * Math.PI / 180;//O=C-O
        double thetaEqCSBOH = 107 * Math.PI / 180;//C-O-H
        double thetaEqCCSBO = 111 * Math.PI / 180;//C-C-O
        double kThetaCCDBO = Kelvin.UNIT.toSim(80600); //force constant [=] K, C-C=O
        double kThetaDBOCSBO = Kelvin.UNIT.toSim(80600);//O=C-O
        double kThetaCSBOH = Kelvin.UNIT.toSim(35200);//C-O-H
        double kThetaCCSBO = Kelvin.UNIT.toSim(70600);//C-C-O
        PotentialGroup pIntra = potentialMaster.makePotentialGroup(1);

        P3BondAngle uBendingCCDBO = new P3BondAngle(space);//C-C=O
        P3BondAngle uBendingDBOCSBO = new P3BondAngle(space);//O=C-O
        P3BondAngle uBendingCSBOH = new P3BondAngle(space);//C-O-H
        P3BondAngle uBendingCCSBO = new P3BondAngle(space);//C-C-O

        uBendingCCDBO.setAngle(thetaEqCCDBO);
        uBendingCCDBO.setEpsilon(kThetaCCDBO);
        uBendingDBOCSBO.setAngle(thetaEqDBOCSBO);
        uBendingDBOCSBO.setEpsilon(kThetaDBOCSBO);
        uBendingCSBOH.setAngle(thetaEqCSBOH);
        uBendingCSBOH.setEpsilon(kThetaCSBOH);
        uBendingCCSBO.setAngle(thetaEqCCSBO);
        uBendingCCSBO.setEpsilon(kThetaCCSBO);

        pIntra.addPotential(uBendingCCDBO, new Atomset3IteratorIndexList(new int[][]{{0, 1, 2}}));//CH3:0, C:1, dBO:2, sBO:3, H:4
        pIntra.addPotential(uBendingDBOCSBO, new Atomset3IteratorIndexList(new int[][]{{2, 1, 3}}));
        pIntra.addPotential(uBendingCSBOH, new Atomset3IteratorIndexList(new int[][]{{1, 3, 4}}));
        pIntra.addPotential(uBendingCCSBO, new Atomset3IteratorIndexList(new int[][]{{0, 1, 3}}));
        potentialMaster.addPotential(pIntra, new ISpecies[]{species});

        //utorsional = c1*[1+cos(phi+f1)]+c2*[1-cos(2*phi)], c1/kB = 630K, c2/kB = 1562.4K, f1=180°(for OCOH), 0°(for CCOH)

        P4BondTorsion p4OCOH = new P4BondTorsion(space, 2 * Kelvin.UNIT.toSim(630.0), Kelvin.UNIT.toSim(-630.0), Kelvin.UNIT.toSim(1562.4), Kelvin.UNIT.toSim(0.0));
        P4BondTorsion p4CCOH = new P4BondTorsion(space, 0, Kelvin.UNIT.toSim(630.0), Kelvin.UNIT.toSim(1562.4), Kelvin.UNIT.toSim(0.0));
        pIntra.addPotential(p4OCOH, new Atomset4IteratorIndexList(new int[][]{{2, 1, 3, 4}}));
        pIntra.addPotential(p4CCOH, new Atomset4IteratorIndexList(new int[][]{{0, 1, 3, 4}}));

        DipoleSourceAceticAcid dipoleSource = new DipoleSourceAceticAcid(space);

        System.out.println("number of molecules " + box.getMoleculeList().size());
        System.out.println("volume " + box.getBoundary().volume());
        System.out.println("Rho " + box.getMoleculeList().size() / box.getBoundary().volume());

        BoxImposePbc pbc = new BoxImposePbc(box, space);
        pbc.setApplyToMolecules(true);
        integrator.getEventManager().addListener(new IntegratorListenerAction(pbc));

        potentialMaster.addPotential(potential, new ISpecies[]{species, species});//two-body
        String configFile = "AceticAcid_NPT_" + numAtoms + "atoms" + temperatureK + "T" + pressureBar + "Bar" + numSteps + "steps";
        if (false && new File(configFile + ".pos").exists()) {
            ConfigurationFile config = new ConfigurationFile(configFile);
            config.initializeCoordinates(box);
        } else {
            ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
            config.initializeCoordinates(box);
        }

    }
 
    public static void main(String[] args) {  	

    	AceticAcidMCParam params = new AceticAcidMCParam();
    	ParseArgs.doParseArgs(params, args);
    	int numAtoms = params.numAtoms;
    	double pressure = params.pressure;
    	double density = params.density;
        double temperature = params.temperature;
        long numSteps = params.numSteps;

        final TestAceticAcidMC3D_NPT sim = new TestAceticAcidMC3D_NPT(numAtoms, pressure, density, temperature, numSteps);
        System.out.println("equilibrium period = " +numSteps/10);//equilibrium period
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps / 10));
System.out.println("equilibrium finished");

MeterDensity rhoMeter = new MeterDensity(sim.box);
        AccumulatorAverage rhoAccumulator = new AccumulatorAverageFixed(1000);//Accumulator that keeps statistics for averaging and error analysis
        DataPump rhoPump = new DataPump(rhoMeter,rhoAccumulator);
        IntegratorListenerAction listener = new IntegratorListenerAction(rhoPump);
        listener.setInterval(1);
        sim.integrator.getEventManager().addListener(listener);
        final MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);//E from integrator
        AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed(10);
        DataPump energyManager = new DataPump(energyMeter, energyAccumulator);
        energyAccumulator.setBlockSize(50);
        IntegratorListenerAction energyListener = new IntegratorListenerAction(energyManager);
        energyListener.setInterval(100);
        sim.integrator.getEventManager().addListener(energyListener);

        //AccumulatorAverage smerAccumulator = new AccumulatorAverageFixed(10);

        AccumulatorAverage energy2Accumulator = new AccumulatorAverageFixed(10);
        final MeterPotentialEnergy energy2Meter = new MeterPotentialEnergy(sim.integrator.getPotentialMaster());//true energy
        energy2Meter.setBox(sim.box);
        DataPump energy2Pump = new DataPump(energy2Meter,energy2Accumulator);
        IntegratorListenerAction energy2Listener = new IntegratorListenerAction(energy2Pump);
        energy2Listener.setInterval(100);
        sim.integrator.getEventManager().addListener(energy2Listener);

        if (false) {
        	SimulationGraphic graphic = new SimulationGraphic(sim,SimulationGraphic.TABBED_PANE,"acetic acid", 1);
        	ISpecies species = sim.getSpecies(0);
            AtomType typeCH3 = species.getTypeByName("CH3");
            AtomType typeC = species.getTypeByName("C");
            AtomType typeDBO = species.getTypeByName("DBO");
            AtomType typeSBO = species.getTypeByName("SBO");
            AtomType typeH = species.getTypeByName("H");
            ((ColorSchemeByType)graphic.getDisplayBox(sim.box).getColorScheme()).setColor(typeCH3, Color.GREEN);
            ((DiameterHashByType)graphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(typeCH3, 2*1.7);
            ((ColorSchemeByType)graphic.getDisplayBox(sim.box).getColorScheme()).setColor(typeC, Color.BLUE);
            ((ColorSchemeByType)graphic.getDisplayBox(sim.box).getColorScheme()).setColor(typeDBO, Color.RED);
            ((ColorSchemeByType)graphic.getDisplayBox(sim.box).getColorScheme()).setColor(typeSBO, Color.YELLOW);
            ((ColorSchemeByType)graphic.getDisplayBox(sim.box).getColorScheme()).setColor(typeH, Color.WHITE);
        	AccumulatorHistory densityHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
            rhoAccumulator.addDataSink(densityHistory, new StatType[]{rhoAccumulator.MOST_RECENT});
            DisplayPlot rhoPlot = new DisplayPlot();
        	densityHistory.setDataSink(rhoPlot.getDataSet().makeDataSink());
        	rhoPlot.setLabel("density");
        	graphic.add(rhoPlot);
        	//AccumulatorHistory smerHistory1 = new AccumulatorHistory(new HistoryCollapsingAverage());
        	//smerAccumulator.addDataSink(smerHistory, new StatType[]{StatType.MOST_RECENT});
        	//DisplayPlot smerPlot1 = new DisplayPlot();
        	//smerHistory1.setDataSink(smerPlot1.getDataSet().makeDataSink());
        	//smerPlot1.setLabel("smer fraction "+associationEnergy+"association Energy");
        	//graphic.add(smerPlot1);
        	DataSourceCountSteps stepCounter = new DataSourceCountSteps(sim.integrator);
        	AccumulatorHistory energyHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
            energyAccumulator.addDataSink(energyHistory, new StatType[]{energyAccumulator.MOST_RECENT});
            DisplayPlot energyPlot = new DisplayPlot();
        	AccumulatorHistory energy2History = new AccumulatorHistory(new HistoryCollapsingAverage());
        	energyHistory.setTimeDataSource(stepCounter);
        	energy2History.setTimeDataSource(stepCounter);
            energy2Accumulator.addDataSink(energy2History, new StatType[]{energy2Accumulator.MOST_RECENT});
            energyHistory.setDataSink(energyPlot.getDataSet().makeDataSink());
        	energyPlot.setLabel("energy");
        	graphic.add(energyPlot);
        	energy2History.setDataSink(energyPlot.getDataSet().makeDataSink());
        	graphic.makeAndDisplayFrame();
        	return;
        }
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));
        
    	CompoundUnit rhoUnit = new CompoundUnit(new Unit[]{Mole.UNIT,Liter.UNIT},new double[]{1,-1});
        double finalDensity = rhoMeter.getDataAsScalar();
        finalDensity = rhoUnit.fromSim(finalDensity);
        System.out.println("next initial Density "+finalDensity);
        System.out.println("numAtom=" +numAtoms);
        double avgDensity = ((DataDouble) ((DataGroup) rhoAccumulator.getData()).getData(rhoAccumulator.AVERAGE.index)).x;//average density
        double errDensity = ((DataDouble) ((DataGroup) rhoAccumulator.getData()).getData(rhoAccumulator.ERROR.index)).x;
        double correlationBlock = ((DataDouble) ((DataGroup) rhoAccumulator.getData()).getData(rhoAccumulator.BLOCK_CORRELATION.index)).x;
        System.out.println("err "+errDensity);
        System.out.println("correlationBlock "+correlationBlock);
        //double avgSmerFraction = ((DataDouble)((DataGroup)smerAccumulator1.getData()).getData(smerAccumulator1.AVERAGE.index)).x;//average density
        System.out.println("average density(sim)= " +avgDensity);
        avgDensity = rhoUnit.fromSim(avgDensity);
        System.out.println("average density(mol/liter)= " +avgDensity);
        //System.out.println("smer fraction of "+associationEnergy+" association energy "+avgSmerFraction);
        
        double Z = pressure/(avgDensity*sim.integrator.getTemperature());
        double avgPE = ((DataDouble) ((DataGroup) energyAccumulator.getData()).getData(energyAccumulator.AVERAGE.index)).x;
        System.out.println("average energy= "+avgPE);
        double finalEnergy = energyMeter.getDataAsScalar();
        double finalEnergy2 = energy2Meter.getDataAsScalar();
        System.out.println("final energy with energymeter "+finalEnergy);
        System.out.println("final energy with energymeter2 "+finalEnergy2);
        avgPE /= numAtoms;
        System.out.println("Z="+Z);
        System.out.println("PE/epsilon="+avgPE+"(sim)");
    	Unit calPerMoles = new CompoundUnit(new Unit[]{Calorie.UNIT,Mole.UNIT},new double[]{1.0,-1.0});
    	System.out.println("PE/epsilon="+calPerMoles.fromSim(avgPE)+"cal/mole");
        double temp = sim.integrator.getTemperature();
        double Cv = ((DataDouble) ((DataGroup) energyAccumulator.getData()).getData(energyAccumulator.STANDARD_DEVIATION.index)).x;
        Cv /= temp;
        Cv *= Cv/numAtoms;
        System.out.println("Cv/k="+Cv);
        
        if (Double.isNaN(Z) || Math.abs(Z+0.25) > 0.15) {
            System.exit(1);
        }
        if (Double.isNaN(avgPE) || Math.abs(avgPE+4.56) > 0.03) {
            System.exit(1);
        }
        if (Double.isNaN(Cv) || Math.abs(Cv-0.61) > 0.45) {  // actual average seems to be 0.51
            System.exit(1);
        }
          
    }
    public static class AceticAcidMCParam extends ParameterBase {
		public int numAtoms = 20;
		public double pressure = 10;//bar
		public double density = 0.2;//1g/cm3=1000/60.052mol/L
		public double temperature = 574.08;//Kelvin
		public long numSteps = 1000000;
	}

}
