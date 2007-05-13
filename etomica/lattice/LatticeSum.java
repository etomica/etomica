package etomica.lattice;

import etomica.data.Data;
import etomica.data.IDataInfo;
import etomica.data.types.DataGroup;
import etomica.space.IVector;
import etomica.util.FunctionGeneral;

public class LatticeSum {

    public LatticeSum(SpaceLattice lattice) {
        this.lattice = lattice;
        IndexIteratorTriangularPermutations iteratorT = new IndexIteratorTriangularPermutations(lattice.getSpace().D());
        setIterator(new IndexIteratorReflecting(iteratorT));
        coreIterator = iteratorT.getCoreIterator();
        coreIterator.setMaxElement(50);
        kVector = lattice.getSpace().makeVector();
    }

    public DataGroup calculateSum(FunctionGeneral function) {
        IDataInfo dataInfo = function.getDataInfo();
        Data sumR = dataInfo.makeData();
        Data sumI = dataInfo.makeData();
        Data work = dataInfo.makeData();
        double diff = 0.0;
        for(int m=1; m<50; m++) {
            coreIterator.setMaxElement(m);
            coreIterator.setMaxElementMin(m);
            iterator.reset();
            while(iterator.hasNext()) {
                IVector site = (IVector)lattice.site(iterator.next());
                Data value = function.f(site);
                double kDotr = kVector.dot(site);
                double ckr = Math.cos(kDotr);
                double skr = Math.sin(kDotr);
                work.E(value);
                work.TE(ckr);
                sumR.PE(work);
                work.E(value);
                work.TE(skr);
                sumI.PE(work);
            }
//            diff = sumR - diff;
//            System.out.println(m + ".\t" + sumR + "\t" + sumI + "\t" + diff/sumR);
//            diff = sumR;
        }
        return new DataGroup(new Data[] {sumR, sumI});
    }
    
    public void setK(IVector k) {
        kVector.E(k);
    }
    
    public IVector getK() {
        return kVector;
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
    
//    public static void main(String[] args) {
//        PrimitiveFcc primitive = new PrimitiveFcc(Space3D.getInstance());
//        primitive.setCubicSize(Math.pow(2.0, 1.0/6.0));//unit density
//        BravaisLattice lattice = new BravaisLattice(primitive);
//        WaveVectorFactory kFactory = new WaveVectorFactoryFcc(primitive);
//        Simulation sim = new Simulation(Space3D.getInstance());
//        Phase phase = new Phase(sim);
//        sim.addPhase(phase);
//        SpeciesSpheresMono species = new SpeciesSpheresMono(sim);
//        sim.getSpeciesManager().addSpecies(species);
//        species.getAgent(phase).setNMolecules(32);
//        int n = 2;
//        phase.getBoundary().setDimensions(new Vector3D(n*primitive.getSize()[0]*Math.sqrt(2),
//                n*primitive.getSize()[1]*Math.sqrt(2),n*primitive.getSize()[2]*Math.sqrt(2)));
//        System.out.println("Density: "+phase.getDensity());
//        kFactory.makeWaveVectors(phase);
//        System.out.println("Number of wave vectors: "+kFactory.getWaveVectors().length);
//        LatticeSum summer = new LatticeSum(lattice);
//        final Potential2SoftSpherical potential = new P2LennardJones(Space3D.getInstance(), 1.0, 1.0);
//        FunctionGeneral function = new FunctionGeneral() {
//            public Data f(Object obj) {
//                Vector3D r = (Vector3D)obj;
//                tensor.x.Ev1v2(r, r);
//                double r2 = r.squared();
//                double dW = potential.du(r2);
//                double d2W = potential.d2u(r2);
//                tensor.TE(1.0/(r2*r2)*(dW - d2W));
//                tensor.x.PEa1Tt1(-dW/r2,identity);
//                return tensor;
//            }
//            public IDataInfo getDataInfo() {
//                return dataInfo;
//            }
//            final DataTensor tensor = new DataTensor(Space3D.getInstance());
//            final DataInfo dataInfo = new DataTensor.DataInfoTensor("",Energy.DIMENSION,Space3D.getInstance());
//            final Tensor3D identity = new Tensor3D(new double[] {1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0});
//        };
//
//        IVector kVector = new Vector3D();
//        kVector.E(0.0);
//        summer.setK(kVector);
//        System.out.println("\n k:"+kVector.toString());
//        Object[] sum = summer.calculateSum(function);
//        DataTensor sum0 = new DataTensor(Space3D.getInstance());
//        sum0.E((Data)sum[0]);
//        Function chopper = new Function.Chop(1e-9);
//        sum0.map(chopper);
//        System.out.println(sum0.toString());
//        for(int i=kFactory.getWaveVectors().length-1; i>=0; i--) {
//            kVector.E(kFactory.getWaveVectors()[i]);
//            summer.setK(kVector);
//            System.out.println("\n k:"+kVector.toString());
//            sum = summer.calculateSum(function);
//            ((Data)sum[0]).map(chopper);
//            ((Data)sum[1]).map(chopper);
//            ((Data)sum[0]).ME(sum0);
//            System.out.println(sum[0].toString());
// //           System.out.println();
////            System.out.println(sum[1].toString());
//        }
//    }
    
    private SpaceLattice lattice;
    private IndexIterator iterator;
    private IndexIteratorTriangular coreIterator;
    private final IVector kVector;
}
