package etomica.modules.rosmosis;
import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtomPositioned;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.IVector;
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
    public Box box;
    public IntegratorVelocityVerlet integrator;
    public P2LennardJones potential11, potential12, potential22;
    public P2LennardJones potentialMM, potentialM1, potentialM2;
    public ActivityIntegrate activityIntegrate;
    
    public ReverseOsmosis(Space space) {
        super(space);
        PotentialMaster potentialMaster = new PotentialMaster(space); //List(this, 2.0);
        
        int N = 360;  //number of atoms, originally 768
        
        //controller and integrator
	    integrator = new IntegratorVelocityVerlet(potentialMaster, getRandom(), 0.01, Kelvin.UNIT.toSim(125));
	    integrator.setIsothermal(false);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatInterval(1);
        integrator.setTimeStep(0.02);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        //solute (1)
        speciesSolute = new SpeciesSpheresMono(this);
        ((ElementSimple)((AtomTypeLeaf)speciesSolute.getMoleculeType()).getElement()).setMass(Dalton.UNIT.toSim(40));
        getSpeciesManager().addSpecies(speciesSolute);
        
        //solvent (2)
        speciesSolvent = new SpeciesSpheresMono(this);
        ((ElementSimple)((AtomTypeLeaf)speciesSolvent.getMoleculeType()).getElement()).setMass(Dalton.UNIT.toSim(40));
        getSpeciesManager().addSpecies(speciesSolvent);

        //membrane
        speciesMembrane = new SpeciesSpheresMono(this);
        ((ElementSimple)((AtomTypeLeaf)speciesMembrane.getMoleculeType()).getElement()).setMass(Dalton.UNIT.toSim(80));
        getSpeciesManager().addSpecies(speciesMembrane);
        
        double epsSolute = Kelvin.UNIT.toSim(125.0);
        double sigSolute = 3.5;
        double epsSolvent = Kelvin.UNIT.toSim(125.0);
        double sigSolvent = 1.75;
        double epsMembrane = Kelvin.UNIT.toSim(12.5);
        double sigMembrane = 0.988 * 3.5; // ???
        
        double xSize = 66+2.0/3.0;// 80 originally
        double yzSize = 21;       // 28 originally
        double rCut = 0.5*yzSize;
        
        //instantiate several potentials for selection in combo-box
	    potential11 = new P2LennardJones(space, sigSolute, epsSolute);
	    P2SoftSphericalTruncatedShifted pTrunc = new P2SoftSphericalTruncatedShifted(potential11, rCut);
        potentialMaster.addPotential(pTrunc,new Species[]{speciesSolute,speciesSolute});
	    
        potential22 = new P2LennardJones(space, sigSolvent, epsSolvent);
        pTrunc = new P2SoftSphericalTruncatedShifted(potential22, rCut);
        potentialMaster.addPotential(pTrunc,new Species[]{speciesSolvent,speciesSolvent});
        
        potential12 = new P2LennardJones(space, 0.5*(sigSolvent+sigSolute), Math.sqrt(epsSolvent*epsSolute));
        pTrunc = new P2SoftSphericalTruncatedShifted(potential12, rCut);
        potentialMaster.addPotential(pTrunc,new Species[]{speciesSolvent,speciesSolute});
        
        potentialMM = new P2LennardJones(space, sigMembrane, epsMembrane);
        pTrunc = new P2SoftSphericalTruncatedShifted(potentialMM, rCut);
        potentialMaster.addPotential(pTrunc,new Species[]{speciesMembrane,speciesMembrane});
        
        potentialM1 = new P2LennardJones(space, 0.5*(sigMembrane+sigSolute), Math.sqrt(epsMembrane*epsSolute));
        pTrunc = new P2SoftSphericalTruncatedShifted(potentialM1, rCut);
        potentialMaster.addPotential(pTrunc,new Species[]{speciesMembrane,speciesSolute});
        
        potentialM2 = new P2LennardJones(space, 0.5*(sigMembrane+sigSolvent), Math.sqrt(epsMembrane*epsSolvent));
        pTrunc = new P2SoftSphericalTruncatedShifted(potentialM2, rCut);
        potentialMaster.addPotential(pTrunc,new Species[]{speciesMembrane,speciesSolvent});


        ((AtomTypeSphere)speciesSolute.getMoleculeType()).setDiameter(sigSolute);
        ((AtomTypeSphere)speciesSolvent.getMoleculeType()).setDiameter(sigSolvent);
        ((AtomTypeSphere)speciesMembrane.getMoleculeType()).setDiameter(sigMembrane);

        //construct box
	    box = new Box(this);
        addBox(box);
        IVector dim = space.makeVector();
        dim.E(new double[]{xSize, yzSize, yzSize});
        box.setDimensions(dim);
        box.setNMolecules(speciesSolute, N);
        new ConfigurationLattice(space.D() == 3 ? new LatticeCubicFcc() : new LatticeOrthorhombicHexagonal()).initializeCoordinates(box);
        AtomSet molecules = box.getMoleculeList(speciesSolute);
        double membraneSize = 5.0;
        double wallX1 = -0.25*xSize;
        double wallX2 =  0.25*xSize;
        IVector pos = space.makeVector();
        for (int i=0; i<molecules.getAtomCount(); ) {
            IAtomPositioned atom = (IAtomPositioned)molecules.getAtom(i);
            pos.E(atom.getPosition());
            double x = pos.x(0);
            if (Math.abs(x - wallX1) < membraneSize*0.5 ||
                Math.abs(x - wallX2) < membraneSize*0.5) {
                // membrane atom
                box.removeMolecule(atom);
                IAtomPositioned newAtom = (IAtomPositioned)speciesMembrane.getMoleculeFactory().makeAtom();
                newAtom.getPosition().E(pos);
                box.addMolecule(newAtom);
            }
            else if (Math.abs(x) < wallX2) {
                // solvent atom
                box.removeMolecule(atom);
                IAtomPositioned newAtom = (IAtomPositioned)speciesSolvent.getMoleculeFactory().makeAtom();
                newAtom.getPosition().E(pos);
                box.addMolecule(newAtom);
            }
            else {
                // atom is OK (solute).  advance to next one
                i++;
            }
        }
        
        P1Tether membraneTether = new P1Tether(box, speciesMembrane);
        membraneTether.setEpsilon(20000);
        potentialMaster.addPotential(membraneTether, new Species[]{speciesMembrane});
        
        integrator.setBox(box);

//        integrator.addIntervalAction(potentialMaster.getNeighborManager(box));
//        integrator.addNonintervalListener(potentialMaster.getNeighborManager(box));
        integrator.addIntervalAction(new BoxImposePbc(box));
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
