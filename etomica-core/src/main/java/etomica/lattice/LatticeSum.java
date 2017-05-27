/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice;

import etomica.space.Vector;
import etomica.data.FunctionData;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataGroup;

public class LatticeSum {

    public LatticeSum(SpaceLattice lattice) {
        this.lattice = lattice;
        IndexIteratorTriangularPermutations iteratorT = new IndexIteratorTriangularPermutations(lattice.getSpace().D());
        setIterator(new IndexIteratorReflecting(iteratorT));
        coreIterator = iteratorT.getCoreIterator();
        setMaxElement(50);
        kVector = lattice.getSpace().makeVector();
    }

    public DataGroup calculateSum(FunctionData<Object> function) {
        IDataInfo dataInfo = function.getDataInfo();
        IData sumR = dataInfo.makeData();
        IData sumI = dataInfo.makeData();
        IData work = dataInfo.makeData();
//        double diff = 0.0;
        for(int m=1; m<=maxElement; m++) {
            coreIterator.setMaxElement(m);
            coreIterator.setMaxElementMin(m);
            iterator.reset();
            while(iterator.hasNext()) {
                Vector site = (Vector)lattice.site(iterator.next());
                IData value = function.f(site);
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
        return new DataGroup(new IData[] {sumR, sumI});
    }
    
    public void setK(Vector k) {
        kVector.E(k);
    }
    
    public Vector getK() {
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

    public int getMaxElement() {
        return maxElement;
    }

    /**
     * Specifies the largest index element that will be included in the lattice sum.
     * Default is 50. 
     */
    public void setMaxElement(int maxElement) {
        this.maxElement = maxElement;
    }
    
//    public static void main(String[] args) {
//        PrimitiveFcc primitive = new PrimitiveFcc(Space3D.getInstance());
//        primitive.setCubicSize(Math.pow(2.0, 1.0/6.0));//unit density
//        BravaisLattice lattice = new BravaisLattice(primitive);
//        WaveVectorFactory kFactory = new WaveVectorFactoryFcc(primitive);
//        Simulation sim = new Simulation(Space3D.getInstance());
//        Box box = new Box(sim);
//        sim.addBox(box);
//        SpeciesSpheresMono species = new SpeciesSpheresMono(sim);
//        sim.getSpeciesManager().addSpecies(species);
//        box.setNMolecules(species, 32);
//        int n = 2;
//        box.getBoundary().setDimensions(new Vector3D(n*primitive.getSize()[0]*Math.sqrt(2),
//                n*primitive.getSize()[1]*Math.sqrt(2),n*primitive.getSize()[2]*Math.sqrt(2)));
//        System.out.println("Density: "+box.getDensity());
//        kFactory.makeWaveVectors(box);
//        System.out.println("Number of wave vectors: "+kFactory.getWaveVectors().length);
//        LatticeSum summer = new LatticeSum(lattice);
//        final Potential2SoftSpherical potential = new P2LennardJones(Space3D.getInstance(), 1.0, 1.0);
//        FunctionData function = new FunctionData() {
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
//        Vector kVector = new Vector3D();
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
    private final Vector kVector;
    private int maxElement;
}
