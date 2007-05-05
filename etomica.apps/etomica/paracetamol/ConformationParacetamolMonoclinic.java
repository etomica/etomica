package etomica.paracetamol;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtomPositioned;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.config.Conformation;
import etomica.space.Space;

/*
 *  Geometry of Published Paracetamol Molecule (Monoclinic)
 * 
 * @author Tai Tan
 */

public class ConformationParacetamolMonoclinic extends Conformation{
	
	private final AtomIteratorArrayListSimple iterator;
	
	public ConformationParacetamolMonoclinic(Space space) {
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
		x =   0.87932;
		y = - 0.01900;
		z = - 1.24329;
		c2.getPosition().E(new double [] {x, y, z});
		
        IAtomPositioned c3 = (IAtomPositioned)iterator.nextAtom();
		x =   2.04562;
		y = - 0.78423;
		z = - 1.40610;
		c3.getPosition().E(new double [] {x, y, z});
		
        IAtomPositioned c4 = (IAtomPositioned)iterator.nextAtom();
		x =   2.36686;
		y = - 1.53602;
		z = - 0.34044;
		c4.getPosition().E(new double [] {x, y, z});
		
        IAtomPositioned o1 = (IAtomPositioned)iterator.nextAtom();
		x =   3.52947;
		y = - 2.25738;
		z = - 0.56693;
		o1.getPosition().E(new double [] {x, y, z});
		
        IAtomPositioned c5 = (IAtomPositioned)iterator.nextAtom();
		x =   1.49876;
		y = - 1.52024;
		z =   0.89670;
		c5.getPosition().E(new double [] {x, y, z});
		
        IAtomPositioned c6 = (IAtomPositioned)iterator.nextAtom();
		x =   0.32387;
		y = - 0.75997;
		z =   1.05799;
		c6.getPosition().E(new double [] {x, y, z});
		
        IAtomPositioned n1 = (IAtomPositioned)iterator.nextAtom();
		x = - 1.19895;
		y =   0.77570;
		z =   0.25787;
		n1.getPosition().E(new double [] {x, y, z});
		
        IAtomPositioned c7 = (IAtomPositioned)iterator.nextAtom();
		x = - 1.85739;
		y =   1.40014;
		z = - 0.64694;
		c7.getPosition().E(new double [] {x, y, z});
		
        IAtomPositioned c8 = (IAtomPositioned)iterator.nextAtom();
		x = - 3.10135;
		y =   2.18139;
		z = - 0.00476;
		c8.getPosition().E(new double [] {x, y, z});
		
        IAtomPositioned o2 = (IAtomPositioned)iterator.nextAtom();
		x = - 1.48999;
		y =   1.38080;
		z = - 1.85311;
		o2.getPosition().E(new double [] {x, y, z});
	
	}

	private static final long serialVersionUID = 1L;
}
