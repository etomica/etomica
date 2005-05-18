package etomica.math.geometry;
import java.util.LinkedList;

/**
 * Representation of a mathematical polyhedron, a 3-dimensional polytope.
 */
public abstract class Polyhedron extends Polytope {
    
    protected Polyhedron(Polygon[] faces) {
        super(faces);
        this.faces = faces;
        this.edges = allEdges(faces);
    }
    
    /**
     * Returns all faces defined by the polyhedron.
     */
    public Polygon[] getFaces() {
        updateVertices();
        return faces;
    }
    
    /**
     * Returns all edges defined by the polyhedron.
     */
    public LineSegment[] getEdges() {
        updateVertices();
        return edges;
    }
        
    /**
     * Returns the sum of the length of the edges of the polyhedron
     */
    public double getPerimeter() {
        updateVertices();
        double sum = 0.0;
        for(int i=0; i<edges.length; i++) {
            sum += edges[i].getLength();
        }
        return sum;
    }

    /**
     * Returns the sum the area of the faces of the polyhedron
     */
    public double getSurfaceArea() {
        updateVertices();
        double sum = 0.0;
        for(int i=0; i<faces.length; i++) {
            sum += faces[i].getArea();
        }
        return sum;
    }
    
    private static LineSegment[] allEdges(Polygon[] faces) {
        LinkedList list = new LinkedList();
        for(int i=0; i<faces.length; i++) {
            LineSegment[] edges= faces[i].getEdges();
            for(int j=0; j<edges.length; j++) {
                if(!list.contains(edges[j])) {
                    list.add(edges[j]);
                }
            }
        }
        return (LineSegment[])list.toArray(new LineSegment[0]);
    }

    protected final LineSegment[] edges;
    protected final Polygon[] faces;
    
 }
