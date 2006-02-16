package testing;

import etomica.simulation.Simulation;
import etomica.math.geometry.Circle;
import etomica.math.geometry.Polytope;
import etomica.space.Vector;
import etomica.space.Space;

public class BoundaryCircular extends etomica.space.Boundary {
	public BoundaryCircular(Simulation sim){
		super(sim.space,new CirclePolytope(sim.space));
		this.space = sim.space;
		v = space.makeVector();
	}
	public Vector centralImage(Vector r){
		return v;
	}
	public void nearestImage(Vector dr){}
	public Vector getDimensions(){
		return v;
	}
	public void setDimensions(Vector r){}
	public Vector randomPosition(){
		return v;
	}
	public float[][] getOverflowShifts(Vector r, double distance){
		return new float[10][10];
	}
	public double[][] imageOrigins(int nShells){
		return new double[1][1];
	}
	public Vector getBoundingBox(){
		return v;
	}
	public Vector getCenter(){
		return v;
	}
	
	
	private final Vector v;
	private final Space space;
}
