package etomica.lattice.crystal;
import etomica.*;
import etomica.lattice.Basis;

/**
 * A 2-atom basis that makes a hcp crystal on a BravaisLattice
 * having a Hexagonal primitive.
 *
 * @author David Kofke
 */
 
 /* History
  * 09/22/02 (DAK) new
  * 01/19/04 (DAK) revised to extend Basis instead of AtomFactory
  */
 
public class BasisHcp implements Basis {
    
    /**
     * Makes a basis using a default that uses AtomFactoryMono
     * for making atom occupying each site.
     * @param space instance of governing space class
     * @param primitive Primitive of the cubic lattice housing this basis.
     * Needed to ensure that separation of basis atoms is consistent with
     * spacing of atoms on lattice.
     */
    public BasisHcp(PrimitiveHexagonal primitive) {
        this.primitive = primitive;
        coordinates = new Space.Vector[2];
        coordinates[0] = Space.makeVector(3);
        coordinates[1] = Space.makeVector(3);
    }
    
    public int size() {
        return 2;
    }
    
    /**
     * Calculates coordinates by multiplying scaled coordinates by scalar
     * size (lattice constant) of the cubic primitive.
     */
    public Space.Vector[] coordinates() {
        Space.Vector[] a = primitive.vectors();
        coordinates[1].Ea1Tv1(factors[0], a[0]);
        coordinates[1].PEa1Tv1(factors[1], a[1]);
        coordinates[1].PEa1Tv1(factors[2], a[2]);
        return coordinates;
    }
    
    private static final double[] factors = new double[] {1./3., 1./3., 0.5};

    private PrimitiveHexagonal primitive;
    private Space.Vector[] coordinates;

}//end of BasisHcp