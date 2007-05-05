package etomica.paracetamol;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtomPositioned;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.config.Conformation;
import etomica.space.Space;

/*
 *  Geometry of Published Paracetamol Molecule
 * 
 * @author Tai Tan
 */

public class ConformationParacetamolOrthorhombic extends Conformation{
	
	private final AtomIteratorArrayListSimple iterator;
	
	public ConformationParacetamolOrthorhombic(Space space) {
		super(space);
		iterator = new AtomIteratorArrayListSimple();	
	}
	
	public void initializePositions(AtomArrayList List){
		
		iterator.setList(List);
		double x = 0.0;
		double y = 0.0;
		double z = 0.0;
		
		iterator.reset();
		
        IAtomPositioned c1 = (IAtomPositioned)iterator.nextAtom();
		c1.getPosition().E(new double [] {x, y, z});
		
        IAtomPositioned c2 = (IAtomPositioned)iterator.nextAtom();
		x = - 0.00881;
		y = - 1.27191;
		z = - 0.56806;
		c2.getPosition().E(new double [] {x, y, z});
		
        IAtomPositioned c3 = (IAtomPositioned)iterator.nextAtom();
		x =   1.14793;
		y = - 2.02553;
		z = - 0.59664;
		c3.getPosition().E(new double [] {x, y, z});
		
        IAtomPositioned c4 = (IAtomPositioned)iterator.nextAtom();
		x =   2.33229;
		y = - 1.53226;
		z = - 0.07189;
		c4.getPosition().E(new double [] {x, y, z});
		
        IAtomPositioned o1 = (IAtomPositioned)iterator.nextAtom();
		x =   3.43035;
		y = - 2.32151;
		z = - 0.13878;
		o1.getPosition().E(new double [] {x, y, z});
		
        IAtomPositioned c5 = (IAtomPositioned)iterator.nextAtom();
		x =   2.34693;
		y = - 0.27036;
		z =   0.49261;
		c5.getPosition().E(new double [] {x, y, z});
		
        IAtomPositioned c6 = (IAtomPositioned)iterator.nextAtom();
		x =   1.18543;
		y =   0.48363;
		z =   0.52922;
		c6.getPosition().E(new double [] {x, y, z});
		
        IAtomPositioned n1 = (IAtomPositioned)iterator.nextAtom();
		x = - 1.13216;
		y =   0.84584;
		z =   0.04563;
		n1.getPosition().E(new double [] {x, y, z});
		
        IAtomPositioned c7 = (IAtomPositioned)iterator.nextAtom();
		x = - 2.43631;
		y =   0.52961;
		z = - 0.17712;
		c7.getPosition().E(new double [] {x, y, z});
		
        IAtomPositioned c8 = (IAtomPositioned)iterator.nextAtom();
		x = - 3.38720;
		y =   1.70470;
		z = - 0.08041;
		c8.getPosition().E(new double [] {x, y, z});
		
        IAtomPositioned o2 = (IAtomPositioned)iterator.nextAtom();
		x = - 2.82585;
		y = - 0.56758;
		z = - 0.46120;
		o2.getPosition().E(new double [] {x, y, z});
	
	}

	private static final long serialVersionUID = 1L;
}
