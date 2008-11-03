package etomica.kmc;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IAtomTypeSphere;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.atom.AtomArrayList;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.Configuration;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationLattice;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.SimulationGraphic;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

/**
 * Simulation using Henkelman's Dimer method to find a saddle point for
 * an adatom on a surface, modeled with LJ.
 * 
 * @author msellers
 *
 */

public class SimKMCLJadatom extends Simulation{

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "DimerLJadatom";
    public final PotentialMaster potentialMaster;
    public IntegratorKMC integratorKMC;
    public IBox box;
    public SpeciesSpheresMono fixed, movable;
    public ActivityIntegrate activityIntegrateKMC;
    public IAtomSet movableSet;
    public IVector adAtomPos;
    

    public SimKMCLJadatom() {
    	super(Space3D.getInstance(), true);
    	potentialMaster = new PotentialMasterMonatomic(this, space);
    	
    //SIMULATION BOX
        box = new Box(new BoundaryRectangularSlit(random, 0, 5, space), space);
        addBox(box);
        
    //SPECIES
    	double sigma = 1.0;
    	fixed = new SpeciesSpheresMono(this, space, new ElementSimple("A", Double.POSITIVE_INFINITY));
        movable = new SpeciesSpheresMono(this, space);  
        getSpeciesManager().addSpecies(fixed);
        getSpeciesManager().addSpecies(movable);
        ((IAtomTypeSphere)fixed.getLeafType()).setDiameter(sigma);
        ((IAtomTypeSphere)movable.getLeafType()).setDiameter(sigma);
    	
        // Must be in same order as the respective species is added to SpeciesManager
        box.setNMolecules(fixed, 256);    	
    	
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(1);
        inflater.actionPerformed();
    	
		potentialMaster.addPotential(new P2LennardJones(space, sigma, 1.0), new IAtomTypeLeaf[]{movable.getLeafType(), fixed.getLeafType()});
		potentialMaster.addPotential(new P2LennardJones(space, sigma, 1.0), new IAtomTypeLeaf[]{movable.getLeafType(), movable.getLeafType()});

		
    //CRYSTAL
        Configuration config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box); 
       
        //ADATOM CREATION AND PLACEMENT
        
        IMolecule iMolecule = movable.makeMolecule();
        box.addMolecule(iMolecule);
        adAtomPos = ((IAtomPositioned)iMolecule.getChildList().getAtom(0)).getPosition();
        //adAtomPos = getSpace().makeVector();
        adAtomPos.setX(0, 3.5);
        adAtomPos.setX(1, -0.30);
        adAtomPos.setX(2, -0.30);
        IVector newBoxLength = space.makeVector();
        newBoxLength.E(box.getBoundary().getDimensions());
        newBoxLength.setX(0, 2.0*adAtomPos.x(0)+2.0);
        box.getBoundary().setDimensions(newBoxLength);

    }
    
    public void setMovableAtoms(double distance, IVector center){
        //distance = distance*distance;
        IVector rij = space.makeVector();
        AtomArrayList movableList = new AtomArrayList();
        IAtomSet loopSet = box.getMoleculeList();
        for (int i=0; i<loopSet.getAtomCount(); i++){
            rij.Ev1Mv2(center,((IAtomPositioned)((IMolecule)loopSet.getAtom(i)).getChildList().getAtom(0)).getPosition());
            if(rij.x(0) > (box.getBoundary().getDimensions().x(0) - 3.0)){continue;}
            //box.getBoundary().nearestImage(rij);
            if(rij.x(0)< distance){
               movableList.add(loopSet.getAtom(i));
            } 
        }
        for (int i=0; i<movableList.getAtomCount(); i++){
            IMolecule newMolecule = movable.makeMolecule();
            box.addMolecule(newMolecule);
            ((IAtomPositioned)newMolecule.getChildList().getAtom(0)).getPosition().E(((IAtomPositioned)((IMolecule)movableList.getAtom(i)).getChildList().getAtom(0)).getPosition());
            box.removeMolecule((IMolecule)movableList.getAtom(i));
        }
        movableSet = box.getMoleculeList(movable);
    }

    
    //Must be run after setMovableAtoms
    public void removeAtoms(double distance, IVector center){
        distance = distance*distance;
        IVector rij = space.makeVector();
        
        IAtomSet loopSet = box.getMoleculeList(movable);
        for (int i=0; i<loopSet.getAtomCount(); i++){
            rij.Ev1Mv2(center,((IAtomPositioned)((IMolecule)loopSet.getAtom(i)).getChildList().getAtom(0)).getPosition());
            box.getBoundary().nearestImage(rij);
            if(rij.squared() < distance){
               box.removeMolecule((IMolecule)loopSet.getAtom(i));
            } 
        }   
    }
    
    public void initializeConfiguration(String fileName){
        ConfigurationFile config = new ConfigurationFile(fileName);
        config.initializeCoordinates(box);
    }
    
    public void integratorKMC(){
        integratorKMC = new IntegratorKMC(this, potentialMaster, 273.15, this.getRandom(), new ISpecies[]{movable}, this.getSpace());
        integratorKMC.setBox(box);
        activityIntegrateKMC = new ActivityIntegrate(integratorKMC);
        getController().addAction(activityIntegrateKMC);
    }
    
    public static void main(String[] args){
       
        final SimKMCLJadatom sim = new SimKMCLJadatom();
        IVector vect = sim.getSpace().makeVector();
        vect.setX(0, 3.5);
        vect.setX(1, 0.0);
        vect.setX(2, 0.0);
        
        sim.initializeConfiguration("0");
        sim.setMovableAtoms(2.0, vect);

        sim.integratorKMC();
        sim.integratorKMC.setInitialStateConditions(-0.06976750944145352, 1.7236382371736393E90);
        sim.integratorKMC.createIntegrators();
        sim.integratorKMC.setSearchLimit(10);
        
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME,1, sim.getSpace(), sim.getController());
        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
        
        sim.integratorKMC.addIntervalAction(simGraphic.getPaintAction(sim.box));
        sim.integratorKMC.integratorDimer.addIntervalAction(simGraphic.getPaintAction(sim.box));
        sim.integratorKMC.integratorMin1.addIntervalAction(simGraphic.getPaintAction(sim.box));
        sim.integratorKMC.integratorMin2.addIntervalAction(simGraphic.getPaintAction(sim.box));
        
        ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
        
        colorScheme.setColor(sim.fixed.getLeafType(),java.awt.Color.gray);
        colorScheme.setColor(sim.movable.getLeafType(),java.awt.Color.red);

        simGraphic.makeAndDisplayFrame(APP_NAME);
    }

}
