package etomica;

import etomica.atom.iterator.AtomIteratorCompound;
import etomica.space.Vector;

/**
 * Fills phase with molecules on a lattice, taking each molecule in successive order
 * from the linked list of molecules.  Takes no special action when list moves from
 * one species to the next.
 * Wall "molecules" are ignored becuase the (super) add method will not add them.
 *
 * Need to improve this to handle different dimensions more elegantly
 */

/* History
 * 01/04/03 (SKK/DAK) added hexagonal lattice for 2D configuration
 * 01/14/03 (DAK) fixed typo in name of getSquareConfig method
 */
 
public class ConfigurationSequential extends Configuration {

	private boolean fill;
	private Vector dimensions;
	private boolean squareConfig;
    
	public ConfigurationSequential(Space space) {
		super(space);
		setFillVertical(true);
		setSquareConfig(false); // hexagonalLattice is Default!!
		dimensions = space.makeVector();
		dimensions.E(Default.BOX_SIZE);
	}
    
	public void setDimensions(Vector dimensions) {this.dimensions.E(dimensions);}
    
	public void setFillVertical(boolean b) {fill = b;}
	public boolean getFillVertical() {return fill;}
    
	public void setSquareConfig(boolean b){ squareConfig = b;}
	public boolean getSquareConfig() {return squareConfig;}
    
	public void initializePositions(AtomIterator[] iterators) {

		AtomIteratorCompound iterator = new AtomIteratorCompound(iterators);//lump 'em all together

		double Lx = dimensions.x(0);
		double Ly = 0.0;
		double Lz = 0.0;
		if(dimensions.length()>1)  Ly = dimensions.x(1);
		if(dimensions.length()>2)  Lz = dimensions.x(2);

		int sumOfMolecules = iterator.size();
         
		if(sumOfMolecules == 0) return;
 //       System.out.println("ConfigurationSequential sumOfMolecules = "+sumOfMolecules);
        
		Vector[] rLat;
        
		switch(space.D()) {
			case 1:
				rLat = lineLattice(sumOfMolecules, Lx);
				break;
			default:
			case 2:
//skkwak
				if(squareConfig){rLat = squareLattice(sumOfMolecules, Lx, Ly, fill);
				} else {rLat = hexagonalLattice(sumOfMolecules,Lx,Ly,fill);}
         
				break;
			case 3:
				rLat = null;
///                rLat = new etomica.lattice.LatticeFCC(sumOfMolecules, Default.BOX_SIZE).positions();//ConfigurationFcc.lattice(sumOfMolecules);
				break;
		}
        
   // Place molecules     
		int i = 0;
		iterator.reset();
		while(iterator.hasNext()) {
			Atom a = iterator.nextAtom();
			if(a.node.parentSpecies() instanceof SpeciesWalls) continue;
			//initialize coordinates of child atoms
			try {//may get null pointer exception when beginning simulation
				a.creator().getConfiguration().initializeCoordinates(a);
			} catch(NullPointerException e) {}
//            a.coord.translateTo(rLat[i]);
            a.coord.position().E(rLat[i]);
 //           System.out.println("configurationsequential: "+rLat[i].toString());
			i++;
		}
   //     initializeMomenta(phase.speciesMaster());
	}
}
