package etomica.paracetamol;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomSet;
import etomica.atom.IAtomLeaf;
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
	
	public void initializePositions(AtomSet List){
		AtomParacetamol molecule = (AtomParacetamol)((IAtomLeaf)List.getAtom(0)).getParentGroup();
		iterator.setList(List);
		
		double x = 0.0;
		double y = 0.0;
		double z = 0.0;
		
		iterator.reset();
		
		AtomLeaf c1 = molecule.C1;
		x =   0.02112;
		y = - 0.05209;
		z =   0.13276;
		c1.getPosition().E(new double [] {x, y, z});
		
		AtomLeaf c2 = molecule.C2;
		x =   0.90043;
		y = - 0.07109;
		z = - 1.11054;
		c2.getPosition().E(new double [] {x, y, z});
		
		AtomLeaf h1 = molecule.H1;
		x =   0.66461;
		y =   0.49952;
		z = - 1.94367;
		h1.getPosition().E(new double [] {x, y, z});
		
		AtomLeaf c3 = molecule.C3;
		x =   2.06673;
		y = - 0.83632;
		z = - 1.27335;
		c3.getPosition().E(new double [] {x, y, z});
		
		AtomLeaf h2 = molecule.H2;
		x =   2.74509;
		y = - 0.85319;
		z = - 2.23162;
		h2.getPosition().E(new double [] {x, y, z});
		
		AtomLeaf c4 = molecule.C4;
		x =   2.38797;
		y = - 1.58811;
		z = - 0.20768;
		c4.getPosition().E(new double [] {x, y, z});
		
		AtomLeaf o1 = molecule.O1;
		x =   3.55058;
		y = - 2.30947;
		z = - 0.43417;
		o1.getPosition().E(new double [] {x, y, z});
		
		AtomLeaf h5 = molecule.H5;
		x =   3.67175;
		y = - 2.77967;
		z =   0.35244;
		h5.getPosition().E(new double [] {x, y, z});
		
		AtomLeaf c5 = molecule.C5;
		x =   1.51987;
		y = - 1.57233;
		z =   1.02946;
		c5.getPosition().E(new double [] {x, y, z});
		
		AtomLeaf h3 = molecule.H3;
		x =   1.74805;
		y = - 2.15126;
		z =   1.86983;
		h3.getPosition().E(new double [] {x, y, z});
		
		AtomLeaf c6 = molecule.C6;
		x =   0.34498;
		y = - 0.81206;
		z =   1.19075;
		c6.getPosition().E(new double [] {x, y, z});
		
		AtomLeaf h4 = molecule.H4;
		x = - 0.32286;
		y = - 0.81763;
		z =   2.15972;
		h4.getPosition().E(new double [] {x, y, z});
		
		AtomLeaf n1 = molecule.N1;
		x = - 1.17784;
		y =   0.72361;
		z =   0.39062;
		n1.getPosition().E(new double [] {x, y, z});
		
		AtomLeaf h6 = molecule.H6;
		x = - 1.59813;
		y =   0.80096;
		z =   1.36267;
		h6.getPosition().E(new double [] {x, y, z});
		
		AtomLeaf c7 = molecule.C7;
		x = - 1.83628;
		y =   1.34805;
		z = - 0.51418;
		c7.getPosition().E(new double [] {x, y, z});
		
		AtomLeaf o2 = molecule.O2;
		x = - 1.46888;
		y =   1.32871;
		z = - 1.72035;
		o2.getPosition().E(new double [] {x, y, z});
		
		AtomLeaf c8 = molecule.C8;
		x = - 3.08023;
		y =   2.12930;
		z =   0.12800;
		c8.getPosition().E(new double [] {x, y, z});
		
		AtomLeaf h7 = molecule.H7;
		x = - 3.60051;
		y =   1.68308;
		z =   1.06275;
		h7.getPosition().E(new double [] {x, y, z});
		
		AtomLeaf h8 = molecule.H8;
		x = - 3.76528;
		y =   2.19509;
		z = - 0.59402;
		h8.getPosition().E(new double [] {x, y, z});
	
		AtomLeaf h9 = molecule.H9;
		x = - 2.77116;
		y =   3.13489;
		z =   0.35058;
		h9.getPosition().E(new double [] {x, y, z});
	
	}

	private static final long serialVersionUID = 1L;
}
