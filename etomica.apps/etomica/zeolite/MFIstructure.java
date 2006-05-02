package etomica.apps.zeolite;
/*
 * Created on Feb 2, 2006
 */
import etomica.lattice.Basis;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
/**
 * 	Class that creates the model for MFI zeolites, of which ZSM-5 is a member.
 * 
 * 
 * 
 * 
 * 
 */
public class MFIstructure implements Basis {
	
	public MFIstructure(){
		//need to determine size (number of atoms)
		size = PerBU.length;
		coordinates = new Vector[size];
		for(int i=0;i<size;i++){
			coordinates[i] = (Vector)PerBU[i].clone();
		}
	}
	public int size(){
		return size;
	}
	/**
	 *Definition of Periodic Building Unit suggested by Building scheme for MEL and MFI
	 *as taken from http://www.iza-structure.org/databases/ModelBuilding/MFI.pdf
	 * 
	 * 
	 */
	private static final Vector3D[] PerBU = new Vector3D[]{
		new Vector3D(0.4224,0.0565,-0.336),
		new Vector3D(0.3072,0.0277,-0.1893),
		new Vector3D(0.2791,0.0613,0.0312),
		new Vector3D(0.1221,0.063,0.0267),
		new Vector3D(0.0713,0.0272,-0.1855),
		new Vector3D(0.1864,0.059,-0.3282),
		new Vector3D(0.4227,-0.1725,-0.3272),
		new Vector3D(0.3078,-0.1302,-0.1855),
		new Vector3D(0.2755,-0.1728,0.0311),
		new Vector3D(0.1206,-0.1731,0.0298),
		new Vector3D(0.0704,-0.1304,-0.182),
		new Vector3D(0.1871,-0.1733,-0.3193)
	};
	
	public Vector[] positions(){
		
		return coordinates;
	}
	
	private int size;
	private Vector[] coordinates;
	private double oldAB;
	private double oldC;
	
}
