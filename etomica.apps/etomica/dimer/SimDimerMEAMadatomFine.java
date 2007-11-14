package etomica.dimer;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.box.Box;
import etomica.chem.elements.Silver;
import etomica.chem.elements.Tin;
import etomica.config.Configuration;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationLattice;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.SimulationGraphic;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.BasisBetaSnA5;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.PrimitiveCubic;
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

public class SimDimerMEAMadatomFine extends Simulation{

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "DimerMEAMadatomSn";
    public final PotentialMaster potentialMaster;
    public IntegratorDimerRT integratorDimer;
    public Box box;
    public IVector [] saddle;
    public SpeciesSpheresMono sn, snFix, snAdatom, ag, agFix, agAdatom, movable;
    public PotentialMEAM potential;
    public ActivityIntegrate activityIntegrateDimer;
    
    public static void main(String[] args){
    	final String APP_NAME = "DimerFineMEAMadatomSn";
    	final SimDimerMEAMadatomFine sim = new SimDimerMEAMadatomFine();

    	sim.activityIntegrateDimer.setMaxSteps(700);
    	
    	SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME);
    	simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
  	
    	ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
    	
    	//Sn

    	colorScheme.setColor(sim.sn.getMoleculeType(),java.awt.Color.gray);
    	colorScheme.setColor(sim.snFix.getMoleculeType(),java.awt.Color.blue);
    	colorScheme.setColor(sim.snAdatom.getMoleculeType(),java.awt.Color.red);
        colorScheme.setColor(sim.movable.getMoleculeType(),java.awt.Color.PINK);
        
        /**
        //Ag
        colorScheme.setColor(sim.ag.getMoleculeType(),java.awt.Color.darkGray);
        colorScheme.setColor(sim.agFix.getMoleculeType(),java.awt.Color.green);
        colorScheme.setColor(sim.agAdatom.getMoleculeType(),java.awt.Color.red);
        colorScheme.setColor(sim.movable.getMoleculeType(),java.awt.Color.PINK);
    	*/
        
    	simGraphic.makeAndDisplayFrame(APP_NAME);
    	
    }
    
    
    public SimDimerMEAMadatomFine() {
    	super(Space3D.getInstance(), true);
    	
    	potentialMaster = new PotentialMaster(space);
    	
    	// Sn
    	Tin tinFixed = new Tin("SnFix", Double.POSITIVE_INFINITY);
    	
        snFix = new SpeciesSpheresMono(this, tinFixed);
        sn = new SpeciesSpheresMono(this, tinFixed);
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
        
        /**
        //Ag
 
        Silver silverFixed = new Silver("AgFix", Double.POSITIVE_INFINITY);
         
        agFix = new SpeciesSpheresMono(this, silverFixed);
        ag = new SpeciesSpheresMono(this, silverFixed);
        agAdatom = new SpeciesSpheresMono(this, Silver.INSTANCE);
        movable = new SpeciesSpheresMono(this, Silver.INSTANCE);
 
        getSpeciesManager().addSpecies(agFix);
        getSpeciesManager().addSpecies(ag);
        getSpeciesManager().addSpecies(agAdatom);
        getSpeciesManager().addSpecies(movable);
        
        ((AtomTypeSphere)agFix.getMoleculeType()).setDiameter(2.8895); 
        ((AtomTypeSphere)ag.getMoleculeType()).setDiameter(2.8895); 
        ((AtomTypeSphere)agAdatom.getMoleculeType()).setDiameter(2.8895);
        ((AtomTypeSphere)movable.getMoleculeType()).setDiameter(2.8895);
        //*/ 
        
        box = new Box(new BoundaryRectangularSlit(space, random, 0, 5));
        addBox(box);
        
        integratorDimer = new IntegratorDimerRT(this, potentialMaster, new Species[]{snAdatom}, "SnAdatom-Fine");     
    	/**
    	//Ag
    	integratorDimer = new IntegratorDimerRT(this, potentialMaster, new Species[]{agAdatom}, "AgAdatom-FINE");
    	 */
        
        activityIntegrateDimer = new ActivityIntegrate(integratorDimer);


    	getController().addAction(activityIntegrateDimer);
    	
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
		
		/**
        //Ag
        box.setNMolecules(agFix, 64);
        box.setNMolecules(ag, 128);
        box.setNMolecules(agAdatom, 0);
        
        potential = new PotentialMEAM(space);
        
        potential.setParameters(agFix, ParameterSetMEAM.Ag);
        potential.setParameters(ag, ParameterSetMEAM.Ag);
        potential.setParameters(agAdatom, ParameterSetMEAM.Ag);
        potential.setParameters(movable, ParameterSetMEAM.Ag);
        
        this.potentialMaster.addPotential(potential, new Species[]{ag, agFix, agAdatom, movable});
        */        
        
    	
    	integratorDimer.setBox(box);	

    	//beta-Sn box
        //The dimensions of the simulation box must be proportional to those of
        //the unit cell to prevent distortion of the lattice.  The values for the 
        //lattice parameters for tin's beta box (a = 5.8314 angstroms, c = 3.1815 
        //angstroms) are taken from the ASM Handbook. 
        
    	//box.setDimensions(new Vector3D(5.8314*3, 5.8314*3, 3.1815*6));
    	//PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.8318, 3.1819);
    	//BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisBetaSnA5());
        
        //Alternatively, using the parameters calculated in Ravelo & Baskes (1997)
        box.setDimensions(new Vector3D(5.92*3, 5.92*3, 3.23*6));
        PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.92, 3.23);
        BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisBetaSnA5());
     
        /**
        //Ag
        box.setDimensions(new Vector3D(4.0863*3, 4.0863*4, 4.0863*4));
        PrimitiveCubic primitive = new PrimitiveCubic(space, 3.6148);
        BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisCubicFcc());
        */
        
        Configuration config = new ConfigurationLattice(crystal);
        config.initializeCoordinates(box); 
        
        
        
        // Sn
        IAtom iAtom = snAdatom.getMoleculeFactory().makeAtom();
        box.getAgent(snAdatom).addChildAtom(iAtom);
        ((IAtomPositioned)iAtom).getPosition().setX(0, 10.0);
        ((IAtomPositioned)iAtom).getPosition().setX(1, 0.1);
        ((IAtomPositioned)iAtom).getPosition().setX(2, -0.1);
        
        /**
        //Ag
        IAtom iAtom = agAdatom.getMoleculeFactory().makeAtom();
        box.getAgent(agAdatom).addChildAtom(iAtom);
        ((IAtomPositioned)iAtom).getPosition().setX(0, 7.0);
        ((IAtomPositioned)iAtom).getPosition().setX(1, 0.9477016722828758);
        ((IAtomPositioned)iAtom).getPosition().setX(2, 1.0709520701043456);
        */
        
        integratorDimer.deltaR = 0.0005;
        integratorDimer.dXl = 10E-5;       
        integratorDimer.deltaXmax = 0.005;
        integratorDimer.dFsq = 0.0001*0.0001;
        integratorDimer.dFrot = 0.01;
        

        String fTail = "-6121967.366434635_saddle_38cut";
        
        ConfigurationFile configFile = new ConfigurationFile(fTail);
        configFile.initializeCoordinates(box);
       
        
         
        IVector rij = space.makeVector();
        AtomArrayList movableList = new AtomArrayList();
        AtomSet loopSet = box.getMoleculeList(sn);
        
        for (int i=0; i<loopSet.getAtomCount(); i++){
            rij.Ev1Mv2(((IAtomPositioned)iAtom).getPosition(),((IAtomPositioned)loopSet.getAtom(i)).getPosition());
          
            if(rij.squared()<38.0){
               movableList.add(loopSet.getAtom(i));
            } 
        }
       
       for (int i=0; i<movableList.getAtomCount(); i++){
           ((IAtomPositioned)box.addNewMolecule(movable)).getPosition().E(((IAtomPositioned)movableList.getAtom(i)).getPosition());
           box.removeMolecule(movableList.getAtom(i));
       }
       
       
       
       
    }
}
