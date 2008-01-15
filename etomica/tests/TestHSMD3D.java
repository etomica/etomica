package etomica.tests;

import etomica.action.ActionIntegrate;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeMolecule;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.data.meter.MeterPressureHard;
import etomica.integrator.IntegratorHard;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardSphere;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

/**
 * Simple hard-sphere molecular dynamics simulation in 3D.
 * Initial configurations at http://rheneas.eng.buffalo.edu/etomica/tests/
 * @author David Kofke
 */
 
public class TestHSMD3D extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public IntegratorHard integrator;
    public SpeciesSpheresMono species, species2;
    public Box box;

    public TestHSMD3D(Space space, int numAtoms) {
        super(space, true);
        PotentialMasterList potentialMaster = new PotentialMasterList(this);
        
        double neighborRangeFac = 1.6;
        double sigma = 1.0;
        // makes eta = 0.35
        double l = 14.4573*Math.pow((numAtoms/2000.0),1.0/3.0);
        potentialMaster.setCellRange(1);
        potentialMaster.setRange(neighborRangeFac*sigma);
        integrator = new IntegratorHard(this, potentialMaster);
        integrator.setTimeStep(0.01);
        integrator.setIsothermal(true);
        ActionIntegrate actionIntegrate = new ActionIntegrate(integrator,false);
        getController().addAction(actionIntegrate);
        actionIntegrate.setMaxSteps(20000000/numAtoms);
        species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);
        getSpeciesManager().removeSpecies(species);
        species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);
        species2 = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species2);
        AtomType type1 = ((AtomTypeMolecule)species.getMoleculeType()).getChildTypes()[0];
        AtomType type2 = ((AtomTypeMolecule)species2.getMoleculeType()).getChildTypes()[0];

        potentialMaster.addPotential(new P2HardSphere(space, sigma, false),new AtomType[]{type1, type1});

        potentialMaster.addPotential(new P2HardSphere(space, sigma, false),new AtomType[]{type1, type2});

        potentialMaster.addPotential(new P2HardSphere(space, sigma, false),new AtomType[]{type2, type2});
        
        box = new Box(this);
        addBox(box);
        box.setNMolecules(species, numAtoms);
        box.setNMolecules(species2, numAtoms/100);
        box.setDimensions(Space.makeVector(new double[]{l,l,l}));
        NeighborListManager nbrManager = potentialMaster.getNeighborManager(box);
        integrator.addIntervalAction(nbrManager);
        integrator.addNonintervalListener(nbrManager);
        integrator.setBox(box);
        ConfigurationFile config = new ConfigurationFile("HSMD3D"+Integer.toString(numAtoms));
        config.initializeCoordinates(box);
        
//        WriteConfiguration writeConfig = new WriteConfiguration("foo",box,1);
//        integrator.addIntervalListener(writeConfig);
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        int numAtoms = 500;
        if (args.length > 0) {
            numAtoms = Integer.valueOf(args[0]).intValue();
        }
        TestHSMD3D sim = new TestHSMD3D(Space3D.getInstance(), numAtoms);

        MeterPressureHard pMeter = new MeterPressureHard(sim.space);
        pMeter.setIntegrator(sim.integrator);
        
        sim.getController().actionPerformed();
        
        double Z = pMeter.getDataAsScalar()*sim.box.volume()/(sim.box.moleculeCount()*sim.integrator.getTemperature());
        System.out.println("Z="+Z);
        
        // compressibility factor for this system should be 5.22
        if (Double.isNaN(Z) || Math.abs(Z-5.22) > 0.03) {
            System.exit(1);
        }
    }
}
