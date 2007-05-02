package etomica.lattice;

import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveFcc;
import etomica.potential.P2LennardJones;
import etomica.potential.Potential2SoftSpherical;
import etomica.space.IVector;
import etomica.space3d.Space3D;
import etomica.util.Function;

public class LatticeSum {

    public LatticeSum(SpaceLattice lattice) {
        this.lattice = lattice;
        IndexIteratorTriangularPermutations iteratorT = new IndexIteratorTriangularPermutations(3);
        setIterator(new IndexIteratorReflecting(iteratorT));
        coreIterator = iteratorT.getCoreIterator();
        coreIterator.setMaxElement(50);
    }

    public double calculateSum(Function function) {
        double sum = 0.0;
        double diff = 0.0;
        for(int k=1; k<150; k++) {
            coreIterator.setMaxElement(k);
            coreIterator.setMaxElementMin(k);
            iterator.reset();
            while(iterator.hasNext()) {
                IVector site = (IVector)lattice.site(iterator.next());
                double r2 = site.squared();
                if(r2 > 0.0) {
                    sum += function.f(r2);
                }
            }
            diff = sum - diff;
            System.out.println(k + ". \t " + sum + "\t"+diff);
            diff = sum;
        }
        return sum;
    }
    
    public IndexIterator getIterator() {
        return iterator;
    }
    private void setIterator(IndexIterator iterator) {
        if(iterator.getD() != lattice.D()) {
            throw new IllegalArgumentException("Given iterator does not produce index arrays of a dimension consistent with the lattice: iterator.getD() = "+iterator.getD()+"; lattice.D() = "+lattice.D());
        }
        this.iterator = iterator;
    }
        
    public SpaceLattice getLattice() {
        return lattice;
    }

    public void setLattice(SpaceLattice lattice) {
        this.lattice = lattice;
    }
    
    public static void main(String[] args) {
        Primitive primitive = new PrimitiveFcc(Space3D.getInstance());
        BravaisLattice lattice = new BravaisLattice(primitive);
        LatticeSum summer = new LatticeSum(lattice);
        final Potential2SoftSpherical potential = new P2LennardJones(Space3D.getInstance(), 1.0, 1.0);
        Function function = new Function() {
            public double f(double x) {
                return potential.u(x);
            }
        };
        double sum = summer.calculateSum(function);
        System.out.println(sum);
    }
    
    private SpaceLattice lattice;
    private IndexIterator iterator;
    private IndexIteratorTriangular coreIterator;

}
