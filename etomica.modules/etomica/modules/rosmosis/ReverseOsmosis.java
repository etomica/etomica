package etomica.modules.rosmosis;
import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Dalton;
import etomica.units.Kelvin;

/**
 * Reverse osmosis simulation, based on simulation code by Sohail Murad.
 * 
 * @author Andrew Schultz
 */
public class ReverseOsmosis extends Simulation {

    private static final long serialVersionUID = 1L;
    public SpeciesSpheresMono speciesSolvent, speciesSolute, speciesMembrane;
    public IBox box;
    public IntegratorVelocityVerlet integrator;
    public P2LennardJones potential11, potential12, potential22;
    public P2LennardJones potentialMM, potentialM1, potentialM2;
    public ActivityIntegrate activityIntegrate;
    public ConfigurationMembrane configMembrane;
    public P1Tether potentialTether;
    
    public ReverseOsmosis(Space _space) {
        super(_space);
        PotentialMaster potentialMaster = new PotentialMaster(space); //List(this, 2.0);
        
        int N = 360;  //number of atoms, originally 768
        
        //controller and integrator
	    integrator = new IntegratorVelocityVerlet(potentialMaster, getRandom(), 0.01, Kelvin.UNIT.toSim(125), space);
	    integrator.setIsothermal(false);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatInterval(1);
        integrator.setTimeStep(0.02);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        //solute (1)
        speciesSolute = new SpeciesSpheresMono(this, space);
        ((ElementSimple)speciesSolute.getLeafType().getElement()).setMass(Dalton.UNIT.toSim(40));
        getSpeciesManager().addSpecies(speciesSolute);
        
        //solvent (2)
        speciesSolvent = new SpeciesSpheresMono(this, space);
        ((ElementSimple)speciesSolvent.getLeafType().getElement()).setMass(Dalton.UNIT.toSim(40));
        getSpeciesManager().addSpecies(speciesSolvent);

        //membrane
        speciesMembrane = new SpeciesSpheresMono(this, space);
        ((ElementSimple)speciesMembrane.getLeafType().getElement()).setMass(Dalton.UNIT.toSim(80));
        getSpeciesManager().addSpecies(speciesMembrane);

        double epsSolute = Kelvin.UNIT.toSim(125.0);
        double sigSolute = 3.5;
        double epsSolvent = Kelvin.UNIT.toSim(125.0);
        double sigSolvent = 0.5*3.5;
        double epsMembrane = Kelvin.UNIT.toSim(12.5);
        double sigMembrane = 0.988 * 3.5; // ???
        
        double xSize = 66+2.0/3.0;// 80 originally
        double yzSize = 21;       // 28 originally
        double rCut = 0.5*yzSize;
        
        //instantiate several potentials for selection in combo-box
	    potential11 = new P2LennardJones(space, sigSolute, epsSolute);
	    P2SoftSphericalTruncatedShifted pTrunc = new P2SoftSphericalTruncatedShifted(potential11, rCut);
        potentialMaster.addPotential(pTrunc,new AtomType[]{speciesSolute.getLeafType(),speciesSolute.getLeafType()});
	    
        potential22 = new P2LennardJones(space, sigSolvent, epsSolvent);
        pTrunc = new P2SoftSphericalTruncatedShifted(potential22, rCut);
        potentialMaster.addPotential(pTrunc,new AtomType[]{speciesSolvent.getLeafType(),speciesSolvent.getLeafType()});
        
        potential12 = new P2LennardJones(space, 0.5*(sigSolvent+sigSolute), Math.sqrt(epsSolvent*epsSolute));
        pTrunc = new P2SoftSphericalTruncatedShifted(potential12, rCut);
        potentialMaster.addPotential(pTrunc,new AtomType[]{speciesSolvent.getLeafType(),speciesSolute.getLeafType()});
        
        potentialMM = new P2LennardJones(space, sigMembrane, epsMembrane);
        pTrunc = new P2SoftSphericalTruncatedShifted(potentialMM, rCut);
        potentialMaster.addPotential(pTrunc,new AtomType[]{speciesMembrane.getLeafType(),speciesMembrane.getLeafType()});
        
        potentialM1 = new P2LennardJones(space, 0.5*(sigMembrane+sigSolute), Math.sqrt(epsMembrane*epsSolute));
        pTrunc = new P2SoftSphericalTruncatedShifted(potentialM1, rCut);
        potentialMaster.addPotential(pTrunc,new AtomType[]{speciesMembrane.getLeafType(),speciesSolute.getLeafType()});
        
        potentialM2 = new P2LennardJones(space, 0.5*(sigMembrane+sigSolvent), Math.sqrt(epsMembrane*epsSolvent));
        pTrunc = new P2SoftSphericalTruncatedShifted(potentialM2, rCut);
        potentialMaster.addPotential(pTrunc,new AtomType[]{speciesMembrane.getLeafType(),speciesSolvent.getLeafType()});


        ((AtomTypeSphere)speciesSolute.getLeafType()).setDiameter(sigSolute);
        ((AtomTypeSphere)speciesSolvent.getLeafType()).setDiameter(sigSolvent);
        ((AtomTypeSphere)speciesMembrane.getLeafType()).setDiameter(sigMembrane);

        //construct box
	    box = new Box(this, space);
        addBox(box);
        IVector dim = space.makeVector();
        dim.E(new double[]{xSize, yzSize, yzSize});
        box.setDimensions(dim);
        configMembrane = new ConfigurationMembrane(this, space);
        configMembrane.setMembraneDim(0);
        configMembrane.setMembraneThicknessPerLayer(10.0/3.0);
        configMembrane.setNumMembraneLayers(2);
        configMembrane.setMembraneWidth(3);
        double density = 0.525/Math.pow(sigSolute, 3);
        configMembrane.setSolutionChamberDensity(density);
        configMembrane.setSolventChamberDensity(density);
        configMembrane.setSpeciesMembrane(speciesMembrane);
        configMembrane.setSpeciesSolute(speciesSolute);
        configMembrane.setSpeciesSolvent(speciesSolvent);
        configMembrane.initializeCoordinates(box);
        
        potentialTether = new P1Tether(box, speciesMembrane, space);
        potentialTether.setEpsilon(20000);
        potentialMaster.addPotential(potentialTether, new AtomType[]{speciesMembrane.getLeafType()});
        
        integrator.setBox(box);

//        integrator.addIntervalAction(potentialMaster.getNeighborManager(box));
//        integrator.addNonintervalListener(potentialMaster.getNeighborManager(box));
        integrator.addIntervalAction(new BoxImposePbc(box, space));
    }
    
    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        if(args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    space = Space3D.getInstance();
                }
            } catch(NumberFormatException e) {}
        }
            
        ReverseOsmosis sim = new ReverseOsmosis(space);
        sim.getController().actionPerformed();
    }//end of main
    
}
