package etomica.dimer;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.box.Box;
import etomica.chem.elements.Copper;
import etomica.chem.elements.Tin;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.SimulationGraphic;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.BasisBetaSnA5;
import etomica.lattice.crystal.PrimitiveTetragonal;
import etomica.meam.ParameterSetMEAM;
import etomica.meam.PotentialMEAM;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.IVector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

/**
 * Simulation using Henkelman's Dimer method to find a saddle point for
 * an adatom of Sn on a surface, modeled with MEAM.
 * 
 * @author msellers
 *
 */

public class EnergyMapMEAMadatomSn extends Simulation{

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "MEAM Md3D";
    public final PotentialMaster potentialMaster;
    public IntegratorEnergyMap integratorMAP;
    public Box box;
    public SpeciesSpheresMono sn, snFix, snAdatom, movable;
    public PotentialMEAM potential;
    public ActivityIntegrate activityIntegrateMAP;
    
    public static void main(String[] args){
        double height1 = 10.0;
        String fileTail1 = ""+height1;
    	final String APP_NAME = "DimerMEAMadatomCu";
    	final EnergyMapMEAMadatomSn sim = new EnergyMapMEAMadatomSn(height1, fileTail1);
    	
    	sim.activityIntegrateMAP.setMaxSteps(1); 
    	SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME);
    	simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

    	sim.integratorMAP.addIntervalAction(simGraphic.getPaintAction(sim.box));
    	
    	ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());

    	//Sn
    
    	colorScheme.setColor(sim.sn.getMoleculeType(),java.awt.Color.gray);
        colorScheme.setColor(sim.snFix.getMoleculeType(),java.awt.Color.blue);
        colorScheme.setColor(sim.snAdatom.getMoleculeType(),java.awt.Color.red);
        colorScheme.setColor(sim.movable.getMoleculeType(),java.awt.Color.PINK);
      
    	
    	simGraphic.makeAndDisplayFrame(APP_NAME);
    	
    }
    
    
    public EnergyMapMEAMadatomSn(double height, String fileTail) {
    	super(Space3D.getInstance(), true);
    	
    	potentialMaster = new PotentialMaster(space);
        // Sn
        Tin tinFixed = new Tin("SnFix", Double.POSITIVE_INFINITY);
        
        snFix = new SpeciesSpheresMono(this, tinFixed);
        sn = new SpeciesSpheresMono(this, Tin.INSTANCE);
        snAdatom = new SpeciesSpheresMono(this, Tin.INSTANCE);
        movable = new SpeciesSpheresMono(this, Tin.INSTANCE);
        
        getSpeciesManager().addSpecies(snFix);
        getSpeciesManager().addSpecies(sn);
        getSpeciesManager().addSpecies(snAdatom);
        getSpeciesManager().addSpecies(movable);
        
        ((AtomTypeSphere)snFix.getMoleculeType()).setDiameter(3.022); 
        ((AtomTypeSphere)sn.getMoleculeType()).setDiameter(3.022); 
        ((AtomTypeSphere)snAdatom.getMoleculeType()).setDiameter(3.022);
        ((AtomTypeSphere)movable.getMoleculeType()).setDiameter(3.022);
        
        box = new Box(new BoundaryRectangularSlit(space, random, 0, 5));
        addBox(box);
            	
		// Sn
        box.setNMolecules(snFix, 72);
        box.setNMolecules(sn, 144);
        box.setNMolecules(snAdatom, 0);
        
        potential = new PotentialMEAM(space);
        
    	potential.setParameters(snFix, ParameterSetMEAM.Sn);
        potential.setParameters(sn, ParameterSetMEAM.Sn);
        potential.setParameters(snAdatom, ParameterSetMEAM.Sn);
        potential.setParameters(movable, ParameterSetMEAM.Sn);
        
        this.potentialMaster.addPotential(potential, new Species[]{sn, snFix, snAdatom, movable});
		//*/		
    	
    	

        // Sn box
     
        box.setDimensions(new Vector3D(5.8314*3, 5.8314*3, 3.1815*6));
        PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.8318, 3.1819);
        BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisBetaSnA5());
        Configuration config = new ConfigurationLattice(crystal);
        config.initializeCoordinates(box); 

        // Sn
        IAtom iAtom = snAdatom.getMoleculeFactory().makeAtom();
        box.getAgent(snAdatom).addChildAtom(iAtom);

         ((IAtomPositioned)iAtom).getPosition().setX(0, height);
         ((IAtomPositioned)iAtom).getPosition().setX(1, -4.37);
         ((IAtomPositioned)iAtom).getPosition().setX(2, -1.6);
         
         /**
         
         IVector rij = space.makeVector();
         AtomArrayList movableList = new AtomArrayList();
         AtomSet loopSet = box.getMoleculeList(sn);
         
         for (int i=0; i<loopSet.getAtomCount(); i++){
             rij.Ev1Mv2(((IAtomPositioned)iAtom).getPosition(),((IAtomPositioned)loopSet.getAtom(i)).getPosition());
           
             if(rij.squared()<32.49){
                movableList.add(loopSet.getAtom(i));
             } 
         }
        
        for (int i=0; i<movableList.getAtomCount(); i++){
            ((IAtomPositioned)box.addNewMolecule(movable)).getPosition().E(((IAtomPositioned)movableList.getAtom(i)).getPosition());
            box.removeMolecule(movableList.getAtom(i));
        }
        
         */
         
         integratorMAP = new IntegratorEnergyMap(this, potentialMaster, iAtom, fileTail);
         integratorMAP.setBox(box);
         activityIntegrateMAP = new ActivityIntegrate(integratorMAP);

         getController().addAction(activityIntegrateMAP);
       
    }
}
