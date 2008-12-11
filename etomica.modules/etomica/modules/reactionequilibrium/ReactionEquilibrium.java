package etomica.modules.reactionequilibrium;

import javax.swing.JPanel;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IAtomTypeSphere;
import etomica.api.IBox;
import etomica.api.IController;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.box.Box;
import etomica.data.meter.MeterTemperature;
import etomica.integrator.IntegratorHard;
import etomica.potential.P1HardPeriodic;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space2d.Space2D;
import etomica.species.SpeciesSpheresMono;

public class ReactionEquilibrium extends Simulation implements AgentSource {

    public IController controller1;
    public JPanel panel = new JPanel(new java.awt.BorderLayout());
    public IntegratorHard integratorHard1;
    public java.awt.Component display;
    public IBox box;
    public etomica.action.SimulationRestart restartAction;
    public boolean initializing = true;
    public MeterTemperature thermometer;
    public SpeciesSpheresMono speciesA;
    public SpeciesSpheresMono speciesB;
    public P2SquareWellBonded AAbonded;
    public P2SquareWellBonded ABbonded;
    public P2SquareWellBonded BBbonded;
    public MeterDimerFraction meterDimerFraction;
    public ActivityIntegrate activityIntegrate;
    private AtomLeafAgentManager agentManager = null;
    public IAtomLeaf[] agents;
    
    public ReactionEquilibrium() {
        super(Space2D.getInstance());
        PotentialMaster potentialMaster = new PotentialMasterMonatomic(this);
        controller1 = getController();

        double diameter = 1.0;

        //controller and integrator
        integratorHard1 = new IntegratorHard(this, potentialMaster, space);
        integratorHard1.setIsothermal(true);

        //construct box
        box = new Box(this, space);
        addBox(box);
        box.setBoundary(new BoundaryRectangularPeriodic(space, random, 30.0));
        integratorHard1.setBox(box);
        speciesA = new SpeciesSpheresMono(this, space);
        speciesB = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(speciesA);
        getSpeciesManager().addSpecies(speciesB);
        ((IAtomTypeSphere)speciesA.getLeafType()).setDiameter(diameter);
        ((IAtomTypeSphere)speciesB.getLeafType()).setDiameter(diameter);
        box.setNMolecules(speciesA, 30);
        box.setNMolecules(speciesB, 30);
        P1HardPeriodic nullPotential = new P1HardPeriodic(space, diameter);
        integratorHard1.setNullPotential(nullPotential, speciesA.getLeafType());
        integratorHard1.setNullPotential(nullPotential, speciesB.getLeafType());

        agentManager = new AtomLeafAgentManager(this,box);

        //potentials
        AAbonded = new P2SquareWellBonded(space, agentManager, 0.5 * diameter, //core
                2.0, //well multiplier
                1.0, true);
        ABbonded = new P2SquareWellBonded(space, agentManager, 0.5 * diameter, //core
                2.0, //well multiplier
                1.0, true);
        BBbonded = new P2SquareWellBonded(space, agentManager, 0.5 * diameter, //core
                2.0, //well multiplier
                1.0, true);
        potentialMaster.addPotential(AAbonded,
                new IAtomTypeLeaf[] { speciesA.getLeafType(), speciesA.getLeafType() });
        potentialMaster.addPotential(ABbonded,
                new IAtomTypeLeaf[] { speciesA.getLeafType(), speciesB.getLeafType() });
        potentialMaster.addPotential(BBbonded,
                new IAtomTypeLeaf[] { speciesB.getLeafType(), speciesB.getLeafType() });

        meterDimerFraction = new MeterDimerFraction(agentManager);
        meterDimerFraction.setSpeciesA(speciesA);
        meterDimerFraction.setBox(box);
        thermometer = new MeterTemperature(box, space.D());
        
        activityIntegrate = new ActivityIntegrate(integratorHard1);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);
        integratorHard1.addIntervalAction(new BoxImposePbc(box, space));
	}
    
    public Class getAgentClass() {
        return IAtomLeaf.class;
    }

    public AtomLeafAgentManager getAgentManager() {
    	return agentManager;
    }

    /**
     * Implementation of Atom.AgentSource interface. Agent is the 
     * bonding partner.
     * 
     * @param a  ignored
     * @return Object always null
     */
    public Object makeAgent(IAtomLeaf a) {
        return null;
    }
    
    public void releaseAgent(Object agent, IAtomLeaf atom) {}


}

