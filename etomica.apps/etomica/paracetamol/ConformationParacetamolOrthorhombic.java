package etomica.paracetamol;

import etomica.atom.AtomSet;
import etomica.atom.IAtomLeaf;
import etomica.atom.IAtomPositioned;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.config.Conformation;
import etomica.space.Space;

/*
 *  Geometry of Published Paracetamol Molecule (Orthorhombic)
 * 
 * @author Tai Tan
 */

public class ConformationParacetamolOrthorhombic extends Conformation{
	
	private final AtomIteratorArrayListSimple iterator;
	
	public ConformationParacetamolOrthorhombic(Space space) {
		super(space);
		iterator = new AtomIteratorArrayListSimple();	
	}
	
	public void initializePositions(AtomSet List){
		AtomParacetamol molecule = (AtomParacetamol)((IAtomLeaf)List.getAtom(0)).getParentGroup();
		iterator.setList(List);
		double x = 0.0;
		double y = 0.0;
		double z = 0.0;
		
		iterator.reset();
		
		IAtomPositioned c1 = molecule.C1;
		x =   0.11531;
		y =   0.08253;
		z =   0.04163;
		c1.getPosition().E(new double [] {x, y, z});
		
		IAtomPositioned c2 = molecule.C2;
		x =   0.10650;
		y = - 1.18938;
		z = - 0.52643;
		c2.getPosition().E(new double [] {x, y, z});
		
		IAtomPositioned h1 = molecule.H1;
		x = - 0.79998;
		y = - 1.58600;
		z = - 0.93135;
		h1.getPosition().E(new double [] {x, y, z});
		
		IAtomPositioned c3 =  molecule.C3;
		x =   1.26324;
		y = - 1.94300;
		z = - 0.55501;
		c3.getPosition().E(new double [] {x, y, z});
		
		IAtomPositioned h2 = molecule.H2;
		x =   1.25933;
		y = - 2.92500;
		z = - 0.99131;
		h2.getPosition().E(new double [] {x, y, z});
		
		IAtomPositioned c4 = molecule.C4;
		x =   2.44760;
		y = - 1.44973;
		z = - 0.03026;
		c4.getPosition().E(new double [] {x, y, z});
		
		IAtomPositioned o1 = molecule.O1;
		x =   3.54566;
		y = - 2.23898;
		z = - 0.09715;
		o1.getPosition().E(new double [] {x, y, z});
		
		IAtomPositioned h5 = molecule.H5;
		x =   4.29280;
		y = - 1.80970;
		z =   0.28439;
		h5.getPosition().E(new double [] {x, y, z});
		
		IAtomPositioned c5 = molecule.C5;
		x =   2.46224;
		y = - 0.18782;
		z =   0.53424;
		c5.getPosition().E(new double [] {x, y, z});
		
		IAtomPositioned h3 = molecule.H3;
		x =   3.37076;
		y =   0.21390;
		z =   0.95074;
		h3.getPosition().E(new double [] {x, y, z});
		
		IAtomPositioned c6 = molecule.C6;
		x =   1.30074;
		y =   0.56616;
		z =   0.57085;
		c6.getPosition().E(new double [] {x, y, z});
		
		IAtomPositioned h4 = molecule.H4;
		x =   1.32879;
		y =   1.54416;
		z =   1.02133;
		h4.getPosition().E(new double [] {x, y, z});
		
		IAtomPositioned n1 = molecule.N1;
		x = - 1.01685;
		y =   0.92837;
		z =   0.08726;
		n1.getPosition().E(new double [] {x, y, z});
		
		IAtomPositioned h6 = molecule.H6;
		x = - 0.82779;
		y =   1.87097;
		z =   0.33401;
		h6.getPosition().E(new double [] {x, y, z});
		
		IAtomPositioned c7 = molecule.C7;
		x = - 2.32100;
		y =   0.61215;
		z = - 0.13549;
		c7.getPosition().E(new double [] {x, y, z});
		
		IAtomPositioned o2 = molecule.O2;
		x = - 2.71054;
		y = - 0.48505;
		z = - 0.41957;
		o2.getPosition().E(new double [] {x, y, z});
		
		IAtomPositioned c8 = molecule.C8;
		x = - 3.27190;
		y =   1.78724;
		z = - 0.03877;
		c8.getPosition().E(new double [] {x, y, z});
		
		IAtomPositioned h7 = molecule.H7;
		x = - 2.95649;
		y =   2.52790;
		z =   0.68746;
		h7.getPosition().E(new double [] {x, y, z});
		
		IAtomPositioned h8 = molecule.H8;
		x = - 4.25132;
		y =   1.41465;
		z =   0.22431;
		h8.getPosition().E(new double [] {x, y, z});
	
		IAtomPositioned h9 = molecule.H9;
		x = - 3.33710;
		y =   2.26662;
		z = - 1.01090;
		h9.getPosition().E(new double [] {x, y, z});
		
	}

	private static final long serialVersionUID = 1L;
}
