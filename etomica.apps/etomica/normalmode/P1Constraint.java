package etomica.normalmode;

import etomica.api.IAtomList;
import etomica.api.IAtomPositioned;
import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.lattice.crystal.Primitive;
import etomica.potential.Potential1;
import etomica.space.ISpace;

public class P1Constraint extends Potential1{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	public P1Constraint(ISpace space, Primitive primitive, IBox box, CoordinateDefinition coordinateDefinition) {
		super(space);
		
		siteIndex = box.getLeafList().getAtomCount();
		latticeSite = space.makeVectorArray(siteIndex);
				
		double neighborRadius = primitive.getSize()[0];
		radiusInner = neighborRadius*neighborRadius/8;
		radiusOuter = neighborRadius*neighborRadius/6;
		
		//Lattice Site Assignment
		for(int i=0; i<siteIndex; i++){
			latticeSite[i] = coordinateDefinition.getLatticePosition(box.getLeafList().getAtom(i));
		}
		
		//Check for neighboring sites
		IVectorMutable distance = space.makeVector();
		neighborSite = new IVectorMutable[siteIndex][siteIndex];
		
		for (int i=0; i<siteIndex; i++){
			for (int j=0; j<siteIndex; j++){
				
				neighborSite[i][j]= space.makeVector();

			}
		}
		//Determine which are the neighboring sites
		for (int i=0; i<siteIndex; i++){
			
			int k=0;
			for(int j=0; j<siteIndex; j++){
			
				distance.Ev1Mv2(latticeSite[j], latticeSite[i]);
				double distCheck = Math.sqrt(distance.squared());
				
				if (distCheck >0.00000001 && distCheck < neighborRadius-0.000001){
					
					neighborSite[i][k].E(latticeSite[j]);
					k++;
				}
			
			}
		}
		
		
	}//End of Constructor

	@Override
	public double energy(IAtomList atoms) {
		
		IVectorMutable posAtom = ((IAtomPositioned)atoms.getAtom(0)).getPosition();
		
		int atomIndex = atoms.getAtom(0).getLeafIndex();
		double d = posAtom.Mv1Squared(latticeSite[atomIndex]);
				
		if (d < radiusInner){
			return 0;
		} 
		
		if (d > radiusOuter){
			return Double.POSITIVE_INFINITY;
		} 
		
		for (int i=0; i<neighborSite[atomIndex].length; i++){
			double dNeighbor = posAtom.Mv1Squared(neighborSite[atomIndex][i]);
			
			if (d > dNeighbor){
				return Double.POSITIVE_INFINITY;
			} 
		}
		
		return 0;
	}

	
	private IVectorMutable[] latticeSite; 
	private IVectorMutable[][] neighborSite;
	private int siteIndex;
	private final double radiusInner, radiusOuter;
	



}
