package etomica.lattice.crystal;
import etomica.*;
import etomica.lattice.Basis;

/**
 * General basis that makes a crystal on a BravaisLattice
 * having a Cubic primitive.  Subclassed to make cubic forms of
 * fcc, bcc, etc. lattices.
 *
 * @author David Kofke
 */
 
public class BasisCubic implements Basis {
    
    /**
     * @param primitive Primitive of the cubic lattice housing this basis.
     * Needed to ensure that separation of basis atoms is consistent with
     * spacing of atoms on lattice.
     * @param scaledCoordinates basis coordinates for the case in which the
     * primitive lattice constant (getSize) is unity.  Given instance is kept
     * for interal representation of basis, so changes to scaledCoordinates
     * will affect the basis.
     */
    public BasisCubic(PrimitiveCubic primitive, Space.Vector[] scaledCoordinates) {
        this.primitive = primitive;
        size = scaledCoordinates.length;
        this.scaledCoordinates = scaledCoordinates;
        coordinates = new Space.Vector[size];
        for(int i=0; i<size; i++) {
            this.coordinates[i] = (Space.Vector)scaledCoordinates[i].clone();
        }
        oldLatticeConstant = 1.0;
    }
    
    public int size() {
        return size;
    }
    
    /**
     * Calculates coordinates by multiplying scaled coordinates by scalar
     * size (lattice constant) of the cubic primitive.
     */
    public Space.Vector[] positions() {
        double latticeConstant = primitive.getSize();
        if(latticeConstant != oldLatticeConstant) {
            for(int i=0; i<size; i++) {
                coordinates[i].Ea1Tv1(latticeConstant, scaledCoordinates[i]);
            }
            oldLatticeConstant = latticeConstant;
        }
        return coordinates;
    }
    
    private int size;
    private Space.Vector[] scaledCoordinates, coordinates;
    private PrimitiveCubic primitive;
    private double oldLatticeConstant;
    
}//end of BasisCubic