package etomica.lattice;
import etomica.*;

/**
 * Primitive group for a monoclinic system.  One primitive-vector
 * angle not normal, and vectors not necessarily of equal length.
 * a != b != c; alpha = gamma = 90deg, beta >= 90deg
 */
public class PrimitiveMonoclinic extends Primitive implements Primitive3D {
    
    private double[] sizeCopy;
    private double a, b, c, beta, sinBeta, cosBeta;
    
    public PrimitiveMonoclinic(Simulation sim) {
        this(sim, 1.0, 1.0, 1.0, rightAngle);
    }
    public PrimitiveMonoclinic(Simulation sim, double a, double b, double c, double beta) {
        super(sim);
        size = new double[sim.space.D()];
        sizeCopy = new double[sim.space.D()];
        setA(a);
        setB(b);
        setC(c);
        setBeta(beta);
    }
    
    public void setA(double a) {
        this.a = a;
        latticeVectors[0].setX(0,a);
        if(lattice != null) lattice.update();
    }
    public double getA() {return a;}
    
    public void setB(double b) {
        this.b = b;
        latticeVectors[1].setX(1,b);
        if(lattice != null) lattice.update();
    }
    public double getB() {return b;}
        
    public void setC(double c) {
        this.c = c;
        latticeVectors[2].setX(0,c*cosBeta);
        latticeVectors[2].setX(2,c*sinBeta);
        if(lattice != null) lattice.update();
    }
    public double getC() {return c;}
    
    public void setAlpha(double t) {}//no adjustment of angle permitted
    public double getAlpha() {return rightAngle;}
    
    public void setBeta(double t) {
        if(t < rightAngle) t = rightAngle;
        if(t > Math.PI) t = Math.PI;
        beta = t;
        cosBeta = Math.cos(beta);
        sinBeta = Math.sin(beta);
        setC(c);
    }
    public double getBeta() {return beta;}
    
    public void setGamma(double t) {}
    public double getGamma() {return rightAngle;}
    

    /**
     * Returns a new, identical instance of this primitive.
     */
    public Primitive copy() {
        return new PrimitiveMonoclinic(simulation, a, b, c, beta);
    }
    
        
    /**
     * Returns a copy of the array of primitive-vector sizes.
     */
     //used by IteratorFactoryCell
    public double[] getSize() {
        for(int i=0; i<D; i++) sizeCopy[i] = size[i];
        return sizeCopy;
    }
    
    public void scaleSize(double scale) {
        setA(a*scale);
        setB(b*scale);
        setC(c*scale);
        if(lattice != null) lattice.update();
    }        
    
    public int[] latticeIndex(Space.Vector q) {
        for(int i=0; i<D; i++) {
            double x = q.x(i)/size[i];
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
        }
        return idx;
    }

    public int[] latticeIndex(Space.Vector q, int[] dimensions) {
        for(int i=0; i<D; i++) {
            double x = q.x(i)/size[i];
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
            while(idx[i] >= dimensions[i]) idx[i] -= dimensions[i];
            while(idx[i] < 0)              idx[i] += dimensions[i];
        }
        return idx;
    }
    
    public Primitive reciprocal() {
        throw new RuntimeException("method PrimitiveOrthorhombic.reciprocal not yet implemented");
    }
    
    public AtomFactory wignerSeitzCellFactory() {
        throw new RuntimeException("method PrimitiveOrthorhombic.wignerSeitzCell not yet implemented");
    }
    
    public AtomFactory unitCellFactory() {
        return new UnitCellFactory(simulation);
    }
    
    private double[] size;
    
    public String toString() {return "Monoclinic";}

///////////////////////////////////////////////////////////////////////////////////////////

public class UnitCellFactory extends AtomFactory {

    AtomType atomType;
    
    public UnitCellFactory(Simulation sim) {
        super(sim);
        setType(new AtomType(this));//default
    }
    
    public boolean isGroupFactory() {return false;}
    
    public void setType(AtomType t) {atomType = t;}
    public AtomType type() {return atomType;}

    /**
     * Builds a single unit cell.
     */
    protected Atom build(AtomTreeNodeGroup parent) {
        return new UnitCell(space, atomType, parent);
    }

    public Atom build(Atom atom) {
        if(!(atom instanceof UnitCell)) throw new IllegalArgumentException("PrimitiveOrthorhombic.UnitCellFactory.build(Atom) attempted using an atom that is not an instance of UnitCell");
        return atom;
    }
    
}//end of UnitCellFactory

///////////////////////////////////////////////////////////////////////////////////////////

/**
 * A orthorhombic unit cell.  Position of the cell is given by the vertex
 * in which each coordinate is minimized.
 */
public class UnitCell extends AbstractCell {
    
    private final Space.Vector delta;

    public UnitCell(Space space, AtomType type, AtomTreeNodeGroup parent) {
        super(space, type, parent);
        delta = space.makeVector();
    }
    /**
     * Dimension of the space occupied by the cell
     */
     public int D() {return space.D();}
     
    /**
     * Returns the volume of the orthorhombic cell.
     */
    public double volume() {
        double sizeN = 1.0;
        for(int i=D()-1; i>=0; i--) sizeN *= size[i];
        return sizeN;
    }
    /**
     * Returns the positions of the vertices relative to the cell position.
     * Absolute positions are obtained by adding the coordinate.position vector.
     * Note that vertices might be computed on-the-fly, with each call of the method, rather than
     * computed once and stored; thus it may be worthwhile to store the values if using them often, 
     * but if doing so be careful to update them if any transformations are done to the lattice.
     */
    public Space.Vector[] vertex() {
        Space.Vector[] vertices = new Space.Vector[space.powerD(2)];//number of vertices is 2^D
        for(int i=0; i<vertices.length; i++) {
            vertices[i] = space.makeVector();
            int mask = 1;
            for(int j=0; j<D; j++) {
                //the bits of i indicate whether corresponding primitive vector is added
                //do bitwise "and" with mask (which has exactly 1 nonzero bit) to see if
                //each bit is 1 or 0
                if((i & mask) != 0) vertices[i].PE(latticeVectors[j]);//actual lattice vectors (not scaled)
                mask *= 2;
            }
            vertices[i].PE(coord.position());//translate to position
        }
        return vertices;
    }
    
    /**
     * Returns <code>true</code> if the given vector lies inside (or on the surface of)
     * this cell, <code>false</code> otherwise.
     */
    public boolean inCell(Space.Vector v) {
        delta.Ev1Mv2(v, coord.position());
        double x;
        switch(D()) {
            case 3: x = delta.x(2);
                    if(x < 0.0 || x > size[2]) return false;
            case 2: x = delta.x(1);//fall through to check all dimensions
                    if(x < 0.0 || x > size[1]) return false;
            case 1: x = delta.x(0);
                    if(x < 0.0 || x > size[0]) return false;
                    break;
            default: throw new RuntimeException("PrimitiveOrthorhombic.UnitCell.inCell not implemented for given dimension");
        }
        return true;
    }
}//end of UnitCell

}//end of PrimitiveOrthorhombic
    
