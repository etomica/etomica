package etomica.models.clathrates.molecularhma;
import Jama.Matrix;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.chem.elements.Oxygen;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeAgentManager;
import etomica.potential.*;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.OrientationFull3D;
import etomica.space3d.Space3D;
import etomica.space3d.Tensor3D;
import etomica.units.dimensions.Null;

//Meter that calculates HMA energy and heat capacity for type 1 clathrate hydrates

public class MeterDADBWaterTIP4P implements IDataSource {

    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final Space space;
    protected final MoleculeAgentManager latticeCoordinates;
    protected final DataSourceScalar meterPE;
    protected final PotentialCalculationForceSum pcForceSum;
    protected final PotentialMaster potentialMaster;
    protected final AtomLeafAgentManager<Vector> forceManager;
    protected final IteratorDirective id;
    protected final Vector dr;
    protected final Vector dri;
    protected final Vector drj;
    protected double latticeEnergy;
    protected final double temperature;
    public static boolean justDADB = true;
    protected final Vector ph1h2;
    protected final Vector q;
    protected final Vector totalForce;
    private final Vector centerMass;
    private final Vector centerMassi;
    private final Vector centerMassj;
    public boolean doTranslation;
    public boolean doRotation;
    protected   Tensor tmpDrr1;
    protected final Tensor[] tmpAtomicTensor3;//46X4=184 atoms dimentional array of 3dim Tensor
public int basisDim;
    protected final Vector torque;
    protected Potential2SoftSpherical potentialLJ;
    protected Potential2SoftSphericalLS potentialLJLS;
    protected EwaldSummation potentialES;
    protected Box box;
    int[] nC; double[] a0; double rCutLJ;

    public MeterDADBWaterTIP4P(int[] nC,double[] a0, double rCutLJ, Potential2SoftSpherical potentialLJ,Potential2SoftSphericalLS potentialLJLS,EwaldSummation potentialES, Space space, DataSourceScalar meterPE, Box box, int basisDim, PotentialMaster potentialMaster, double temperature, MoleculeAgentManager latticeCoordinates) {

            this.nC=nC;
        this.a0=a0;
        this.rCutLJ=rCutLJ;
        this.potentialLJ = potentialLJ;
        this.potentialLJLS = potentialLJLS;
        this.potentialES = potentialES;
        int nData = 15;
        data = new DataDoubleArray(nData);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.space = space;
        this.basisDim = basisDim;
        int atomsPerMol = box.getLeafList().size() / box.getMoleculeList().size();
        this.latticeCoordinates = latticeCoordinates;
        this.meterPE = meterPE;
        this.potentialMaster = potentialMaster;
        id = new IteratorDirective();
        pcForceSum = new PotentialCalculationForceSum();
        forceManager = new AtomLeafAgentManager<>(a -> space.makeVector(), latticeCoordinates.getBox());
        pcForceSum.setAgentManager(forceManager);
        dr = space.makeVector();
        dri = space.makeVector();
        drj = space.makeVector();
        MeterPotentialEnergy meterPE2 = new MeterPotentialEnergy(potentialMaster);
        meterPE2.setBox(latticeCoordinates.getBox());
        latticeEnergy = meterPE2.getDataAsScalar();
        this.temperature = temperature;
        ph1h2 = space.makeVector();
        q = space.makeVector();
        totalForce = space.makeVector();
        centerMass = space.makeVector();
        centerMassi = space.makeVector();
        centerMassj = space.makeVector();
        torque = space.makeVector();
         this.box = box;

         this.tmpAtomicTensor3 = new Tensor[basisDim*atomsPerMol];//46X4=184 atoms dimensional array of 3D Tensor
        for (int i=0; i<basisDim*atomsPerMol ;i++){
            tmpAtomicTensor3[i] = space.makeTensor();
        }
    }

    public void setLatticeEnergy(double newLatticeEnergy) {
        latticeEnergy = newLatticeEnergy;
    }


    //aTensor is related to atomicTensorAtomicPair

    public IData getData() {
        Box box = latticeCoordinates.getBox();
        double[] x = data.getData();
        double pe=meterPE.getDataAsScalar();
        double pesquare=meterPE.getDataAsScalar()*meterPE.getDataAsScalar();
        double x0 = meterPE.getDataAsScalar() - latticeEnergy;
        pcForceSum.reset();
        potentialMaster.calculate(box, id, pcForceSum);
        IMoleculeList molecules = box.getMoleculeList();
        double fdotdeltar = 0;
        double orientationSum = 0;
        int N = molecules.size();
        double boltzmann=1.0; //ANSWER WOULD BE IN PER TEMPERATURE UNITS OR SO NOT IN SI UNITS
 //       double dNminus1by2beta=(1.5 * (N - 1) + (1.5* N))*temperature*boltzmann;
 //       AtomLeafAgentManager<Vector> atomLatticeCoordinates = new AtomLeafAgentManager<>(new AtomSiteSource(space), box, Vector.class);

        for (int i = 0; i < molecules.size(); i++) {
            IMolecule molecule = molecules.get(i);
            IAtomList leafList = molecule.getChildList();

            Vector h1Force = forceManager.getAgent(leafList.get(0));
            Vector h2Force = forceManager.getAgent(leafList.get(1));
            Vector oForce = forceManager.getAgent(leafList.get(2));
            Vector mForce = forceManager.getAgent(leafList.get(3));

            totalForce.E(h1Force);
            totalForce.PE(h2Force);
            totalForce.PE(mForce);
            totalForce.PE(oForce);
            Vector h1 = leafList.get(0).getPosition();
            Vector h2 = leafList.get(1).getPosition();
            Vector o = leafList.get(2).getPosition();
            Vector m = leafList.get(3).getPosition();
            double hMass = leafList.get(0).getType().getMass();
            double oMass = leafList.get(2).getType().getMass();

            centerMass.Ea1Tv1(hMass, h1);
            centerMass.PEa1Tv1(hMass, h2);
            centerMass.PEa1Tv1(oMass, o);
            centerMass.TE(1 / (2 * hMass + oMass));

            dr.Ev1Mv2(m, centerMass);
            torque.E(mForce);
            torque.XE(dr);   //this is -ve of torque - so in final formula use a minus in front of torque
            q.E(torque); //this is -ve of torque - so in final formula use a minus in front of torque
            dr.Ev1Mv2(h1, centerMass);
            torque.E(h1Force);
            torque.XE(dr);
            q.PE(torque);
            dr.Ev1Mv2(h2, centerMass);
            torque.E(h2Force);
            torque.XE(dr);
            q.PE(torque);
            dr.Ev1Mv2(o, centerMass);
            torque.E(oForce);
            torque.XE(dr);
            q.PE(torque); //this is -ve of torque - so in final formula use a minus in front of torque
            Vector lPos = ((MoleculeSiteSource.LatticeCoordinate) latticeCoordinates.getAgent(molecule)).position;
            dr.Ev1Mv2(centerMass, lPos);
            fdotdeltar += totalForce.dot(dr);

            OrientationFull3D or = ((MoleculeSiteSource.LatticeCoordinate) latticeCoordinates.getAgent(molecule)).orientation;
            Vector a0 = or.getDirection();//om
            Vector a1 = or.getSecondaryDirection();//h1h2
            dr.Ev1Mv2(h2, h1);
            dr.normalize();
            Vector axis = space.makeVector();
            Vector om = space.makeVector();
            Vector hh = space.makeVector();
            Vector a2 = space.makeVector();

            a2.E(a0);
            a2.XE(a1);
            a2.normalize();
            double[][] array = new double[3][3];
            a0.assignTo(array[0]);
            a1.assignTo(array[1]);
            a2.assignTo(array[2]);
            Matrix a = new Matrix(array).transpose();

            om.Ev1Mv2(m, o);
            om.normalize();
            hh.Ev1Mv2(h2, h1);
            hh.normalize();
            a2.E(om);
            a2.XE(hh);
            a2.normalize();
            double[][] array1 = new double[3][3];
            om.assignTo(array1[0]);
            hh.assignTo(array1[1]);
            a2.assignTo(array1[2]);
            Matrix newa = new Matrix(array1).transpose();

            // newa = rotationMatrix * a so rotationMatrix = newa * a^-1
            a = a.inverse();
            Matrix matrix = newa.times(a);

            //https://en.wikipedia.org/wiki/Rotation_matrix#Determining_the_axis
            double d, b, c, g, h, f;
            b = matrix.get(0, 1);
            c = matrix.get(0, 2);
            d = matrix.get(1, 0);
            f = matrix.get(1, 2);
            g = matrix.get(2, 0);
            h = matrix.get(2, 1);
            axis.setX(0, h - f);
            axis.setX(1, c - g);
            axis.setX(2, d - b);
            double beta = Math.asin(0.5 * Math.sqrt(axis.squared()));
            if (beta == 0) continue;
            axis.normalize();


            double DUDT = -q.dot(axis);
            double denominator = 1 - Math.cos(beta);
            if (denominator == 0) continue;
             orientationSum += 1.5 * (beta - Math.sin(beta)) / denominator * DUDT;  //TRY ANOTHER 3/2 THAT U WERE GETTING FOR TORQUE ..DERIVATION
        }

        final double[] a0_sc = new double[]{a0[0] * nC[0], a0[1] * nC[1], a0[2] * nC[2]};
        final int[] nLJShells = new int[]{(int) Math.ceil(rCutLJ / a0_sc[0] - 0.49999), (int) Math.ceil(rCutLJ / a0_sc[1] - 0.49999), (int) Math.ceil(rCutLJ / a0_sc[2] - 0.49999)};
        final double rCutLJ2 = rCutLJ * rCutLJ;
        final Space space = Space3D.getInstance();

        LatticeSumMolecularCrystal.AtomicTensorAtomicPair atomicTensorAtomicPair = new LatticeSumMolecularCrystal.AtomicTensorAtomicPair() {
            final Tensor identity = new Tensor3D(new double[][]{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}});
            Vector Lxyz = space.makeVector();
            Vector dr = space.makeVector();
            Vector drTmp = space.makeVector();
            Tensor D3ESLJ = space.makeTensor();
            Tensor tmpD3LJ = space.makeTensor();
            //second derivative
            //dr == r0 - r1
            public Tensor atomicTensor(IAtom atom0, IAtom atom1) {
                D3ESLJ.E(0);
                //LJ
                if (atom0.getType().getElement() == Oxygen.INSTANCE && atom1.getType().getElement() == Oxygen.INSTANCE) {
                    dr.Ev1Mv2(atom0.getPosition(), atom1.getPosition());
                    box.getBoundary().nearestImage(dr);

                    //LS of LJ
                    for (int nx = -nLJShells[0]; nx <= nLJShells[0]; nx++) {
                        Lxyz.setX(0, nx * a0_sc[0]);
                        for (int ny = -nLJShells[1]; ny <= nLJShells[1]; ny++) {
                            Lxyz.setX(1, ny * a0_sc[1]);
                            for (int nz = -nLJShells[2]; nz <= nLJShells[2]; nz++) {
                                Lxyz.setX(2, nz * a0_sc[2]);
                                drTmp.Ev1Pv2(dr, Lxyz);
                                double dr2Tmp = drTmp.squared();
                                if (dr2Tmp > rCutLJ2 || dr2Tmp == 0) continue;
                                tmpD3LJ.Ev1v2(drTmp, drTmp);
                                double dW = potentialLJ.du(dr2Tmp);
                                double d2W = potentialLJ.d2u(dr2Tmp);//2nd der
                                tmpD3LJ.TE(1.0 / (dr2Tmp * dr2Tmp) * (dW - d2W));
                                tmpD3LJ.PEa1Tt1(-dW / dr2Tmp, identity);
                                D3ESLJ.PE(tmpD3LJ);
                            }
                        }
                    }
                    //ES
                } else if (atom0.getType().getElement() != Oxygen.INSTANCE && atom1.getType().getElement() != Oxygen.INSTANCE) {
                    D3ESLJ.PE(potentialES.secondDerivative(atom0, atom1));
                }
                return D3ESLJ;
            }
        }; // End AtomicTensorAtomicPair

        LatticeSumMolecularCrystal latticeSumMolecularCrystal = new LatticeSumMolecularCrystal(potentialMaster, box, space, basisDim);
        latticeSumMolecularCrystal.reset();

//second derivative related to atomic tensor
                double fac = (doTranslation ? 1.5 : 0) * (N - 1) + (doRotation ? 1.5 : 0) * N;
 //               x[0] = (x0 + latticeEnergy) + (fac * temperature) + 0.5 * fdotdeltar + orientationSum;
//x[0]=pe
                 x[0]=(x0)  + 0.5 * fdotdeltar + orientationSum;    //anh energy mapped
  //               x[1]=(x0)  + 0.5 * fdotdeltar + 1.5*orientationSum; THIS ONE DIDN'T WORK
                 x[1]=x[0]*x[0];   ////anh energysquared mapped

        //x[1]=anharmonicenergy;

        double orientationSumi=0.0;
        double orientationSumj=0.0;
        double sum=0.0;

        //        double delrdotphidotdelrwithuasthetaby2=0.0;
//        double delrdotphidotdelrwithuas2theta=0.0;

        for (int i = 0; i < molecules.size(); i++) {

            for (int j = 0; j < molecules.size(); j++) {
                if(i==j) continue;

                Tensor D6mol = latticeSumMolecularCrystal.atomicToMolecularD(atomicTensorAtomicPair, molecules.get(i), molecules.get(j));
//above method will call         atomicTensorAtomicPair.atomicTensor(); //

                        Tensor3D D3tt  = new Tensor3D();	Tensor3D D3tr  = new Tensor3D();
                        Tensor3D D3rt  = new Tensor3D();	Tensor3D D3rr  = new Tensor3D();

                        for (int ii=0; ii<3; ii++) {
                            for (int jj = 0; jj < 3; jj++) {
                                D3tt.setComponent(ii, jj, (D6mol.component(ii, jj)));
                                D3tr.setComponent(ii, jj, D6mol.component(ii, jj + 3));
                                D3rt.setComponent(ii, jj, D6mol.component(ii + 3, jj));
                                D3rr.setComponent(ii, jj, D6mol.component(ii + 3, jj + 3));
                     //           System.out.println("Dtt "+i+" "+j+" "+ii+" "+jj+" "+D3tt.component(ii,jj));
                     //           System.out.println("Dtr "+i+" "+j+" "+ii+" "+jj+" "+D3tr.component(ii,jj));
                     //           System.out.println("Drt "+i+" "+j+" "+ii+" "+jj+" "+D3rt.component(ii,jj));
                     //           System.out.println("Drr "+i+" "+j+" "+ii+" "+jj+" "+D3rr.component(ii,jj));

                            }
                        }
                                Vector h1i = molecules.get(i).getChildList().get(0).getPosition();
                                Vector h2i = molecules.get(i).getChildList().get(1).getPosition();
                                Vector oi = molecules.get(i).getChildList().get(2).getPosition();
                                Vector mi = molecules.get(i).getChildList().get(3).getPosition();
                                double hMass = molecules.get(i).getChildList().get(0).getType().getMass();
                                double oMass = molecules.get(i).getChildList().get(2).getType().getMass();
                                centerMassi.Ea1Tv1(hMass, h1i);
                                centerMassi.PEa1Tv1(hMass, h2i);
                                centerMassi.PEa1Tv1(oMass, oi);
                                centerMassi.TE(1 / (2 * hMass + oMass));
                                Vector lPosi = ((MoleculeSiteSource.LatticeCoordinate) latticeCoordinates.getAgent(molecules.get(i))).position;
                                dri.Ev1Mv2(centerMassi, lPosi);

                                Vector h1j = molecules.get(j).getChildList().get(0).getPosition();
                                Vector h2j = molecules.get(j).getChildList().get(1).getPosition();
                                Vector oj = molecules.get(j).getChildList().get(2).getPosition();
                                Vector mj = molecules.get(j).getChildList().get(3).getPosition();
                                centerMassj.Ea1Tv1(hMass, h1j);
                                centerMassj.PEa1Tv1(hMass, h2j);
                                centerMassj.PEa1Tv1(oMass, oj);
                                centerMassj.TE(1 / (2 * hMass + oMass));
                                Vector lPosj = ((MoleculeSiteSource.LatticeCoordinate) latticeCoordinates.getAgent(molecules.get(j))).position;
                                drj.Ev1Mv2(centerMassj, lPosj);


                                OrientationFull3D ori = ((MoleculeSiteSource.LatticeCoordinate) latticeCoordinates.getAgent(molecules.get(i))).orientation;
                                Vector a0i = ori.getDirection();//om
                                Vector a1i = ori.getSecondaryDirection();//h1h2
                                Vector axisi = space.makeVector();
                                Vector omi = space.makeVector();
                                Vector hhi = space.makeVector();
                                Vector a2i = space.makeVector();
                                a2i.E(a0i);
                                a2i.XE(a1i);
                                a2i.normalize();
                                double[][] arrayi = new double[3][3];
                                a0i.assignTo(arrayi[0]);
                                a1i.assignTo(arrayi[1]);
                                a2i.assignTo(arrayi[2]);
                                Matrix ai = new Matrix(arrayi).transpose();
                                omi.Ev1Mv2(mi, oi);
                                omi.normalize();
                                hhi.Ev1Mv2(h2i, h1i);
                                hhi.normalize();
                                a2i.E(omi);
                                a2i.XE(hhi);
                                a2i.normalize();
                                double[][] array1i = new double[3][3];
                                omi.assignTo(array1i[0]);
                                hhi.assignTo(array1i[1]);
                                a2i.assignTo(array1i[2]);
                                Matrix newai = new Matrix(array1i).transpose();
                                // newa = rotationMatrix * a so rotationMatrix = newa * a^-1
                                ai = ai.inverse();
                                Matrix matrixi = newai.times(ai);
                                double di, bi, ci, gi, hi, fi;
                                bi = matrixi.get(0, 1);
                                ci = matrixi.get(0, 2);
                                di = matrixi.get(1, 0);
                                fi = matrixi.get(1, 2);
                                gi = matrixi.get(2, 0);
                                hi = matrixi.get(2, 1);

                axisi.setX(0, hi - fi);
                axisi.setX(1, ci - gi);
                axisi.setX(2, di - bi);
                double betai = Math.asin(0.5 * Math.sqrt(axisi.squared()));
                if (betai == 0) continue;
                double denominatori = 1 - Math.cos(betai);
                if (denominatori == 0) continue;
                orientationSumi = 1.5 * (betai - Math.sin(betai)) / denominatori ;
                axisi.normalize();

          OrientationFull3D orj = ((MoleculeSiteSource.LatticeCoordinate) latticeCoordinates.getAgent(molecules.get(j))).orientation;
                                Vector a0j = orj.getDirection();//om
                                Vector a1j = orj.getSecondaryDirection();//h1h2
                                Vector axisj = space.makeVector();
                                Vector omj = space.makeVector();
                                Vector hhj = space.makeVector();
                                Vector a2j = space.makeVector();
                                a2j.E(a0j);
                                a2j.XE(a1j);
                                a2j.normalize();
                                double[][] arrayj = new double[3][3];
                                a0j.assignTo(arrayj[0]);
                                a1j.assignTo(arrayj[1]);
                                a2j.assignTo(arrayj[2]);
                                Matrix aj = new Matrix(arrayj).transpose();
                                omj.Ev1Mv2(mj, oj);
                                omj.normalize();
                                hhj.Ev1Mv2(h2j, h1j);
                                hhj.normalize();
                                a2j.E(omj);
                                a2j.XE(hhj);
                                a2j.normalize();
                                double[][] array1j = new double[3][3];
                                omj.assignTo(array1j[0]);
                                hhj.assignTo(array1j[1]);
                                a2j.assignTo(array1j[2]);
                                Matrix newaj = new Matrix(array1j).transpose();
                                aj = aj.inverse();
                                Matrix matrixj = newaj.times(aj);
                                double dj, bj, cj, gj, hj, fj;
                                bj = matrixj.get(0, 1);
                                cj = matrixj.get(0, 2);
                                dj = matrixj.get(1, 0);
                                fj = matrixj.get(1, 2);
                                gj = matrixj.get(2, 0);
                                hj = matrixj.get(2, 1);

                axisj.setX(0, hj - fj);
                axisj.setX(1, cj - gj);
                axisj.setX(2, dj - bj);
                double betaj = Math.asin(0.5 * Math.sqrt(axisj.squared()));
                if (betaj == 0) continue;
                double denominatorj = 1 - Math.cos(betaj);
                if (denominatorj == 0) continue;
                orientationSumj = 1.5 * (betaj - Math.sin(betaj)) / denominatorj ;
                axisj.normalize();

                Vector drjj=space.makeVector();
                Vector drjj2=space.makeVector();

                Vector axisjj=space.makeVector();
                Vector axisjj2=space.makeVector();

                drjj.E(drj);
                drjj2.E(drj);

                axisjj.E(axisj);
                axisjj2.E(axisj);

                D3tt.transform(drjj);
                sum += dri.dot(drjj);

                D3tr.transform(axisjj);
                sum += orientationSumj*dri.dot(axisjj)*Math.sin(betaj);

                D3rt.transform(drjj2);
                sum += orientationSumi*axisi.dot(drjj2)*Math.sin(betai);

                D3rr.transform(axisjj2);
                sum += Math.sin(betai)*orientationSumj*orientationSumi*axisi.dot(axisjj2)*Math.sin(betaj);    //sum is deltar.phi.delr/theta sum

                    }


            for (int j = 0; j < molecules.size(); j++) {
                if(i!=j) continue;
                //              System.out.println("fine till here");
                Tensor D6mol = latticeSumMolecularCrystal.atomicToMolecularD(atomicTensorAtomicPair, molecules.get(i), molecules.get(j));
//above method will call         atomicTensorAtomicPair.atomicTensor(); //

                Tensor3D D3tt  = new Tensor3D();	Tensor3D D3tr  = new Tensor3D();
                Tensor3D D3rt  = new Tensor3D();	Tensor3D D3rr  = new Tensor3D();

                for (int ii=0; ii<3; ii++) {
                    for (int jj = 0; jj < 3; jj++) {
                        D3tt.setComponent(ii, jj, (D6mol.component(ii, jj)));
                        D3tr.setComponent(ii, jj, D6mol.component(ii, jj + 3));
                        D3rt.setComponent(ii, jj, D6mol.component(ii + 3, jj));
                        D3rr.setComponent(ii, jj, D6mol.component(ii + 3, jj + 3));
                 //       System.out.println("Dtt "+i+" "+j+" "+ii+" "+jj+" "+D3tt.component(ii,jj));
                 //       System.out.println("Dtr "+i+" "+j+" "+ii+" "+jj+" "+D3tr.component(ii,jj));
                 //       System.out.println("Drt "+i+" "+j+" "+ii+" "+jj+" "+D3rt.component(ii,jj));
                 //       System.out.println("Drr "+i+" "+j+" "+ii+" "+jj+" "+D3rr.component(ii,jj));

                    }
                }


                Vector h1i = molecules.get(i).getChildList().get(0).getPosition();
                Vector h2i = molecules.get(i).getChildList().get(1).getPosition();
                Vector oi = molecules.get(i).getChildList().get(2).getPosition();
                Vector mi = molecules.get(i).getChildList().get(3).getPosition();
                double hMass = molecules.get(i).getChildList().get(0).getType().getMass();
                double oMass = molecules.get(i).getChildList().get(2).getType().getMass();
                centerMassi.Ea1Tv1(hMass, h1i);
                centerMassi.PEa1Tv1(hMass, h2i);
                centerMassi.PEa1Tv1(oMass, oi);
                centerMassi.TE(1 / (2 * hMass + oMass));
                Vector lPosi = ((MoleculeSiteSource.LatticeCoordinate) latticeCoordinates.getAgent(molecules.get(i))).position;
                dri.Ev1Mv2(centerMassi, lPosi);
                drj.E(dri);

                OrientationFull3D ori = ((MoleculeSiteSource.LatticeCoordinate) latticeCoordinates.getAgent(molecules.get(i))).orientation;
                Vector a0i = ori.getDirection();//om
                Vector a1i = ori.getSecondaryDirection();//h1h2
                Vector axisi = space.makeVector();
                Vector omi = space.makeVector();
                Vector hhi = space.makeVector();
                Vector a2i = space.makeVector();
                a2i.E(a0i);
                a2i.XE(a1i);
                a2i.normalize();
                double[][] arrayi = new double[3][3];
                a0i.assignTo(arrayi[0]);
                a1i.assignTo(arrayi[1]);
                a2i.assignTo(arrayi[2]);
                Matrix ai = new Matrix(arrayi).transpose();
                omi.Ev1Mv2(mi, oi);
                omi.normalize();
                hhi.Ev1Mv2(h2i, h1i);
                hhi.normalize();
                a2i.E(omi);
                a2i.XE(hhi);
                a2i.normalize();
                double[][] array1i = new double[3][3];
                omi.assignTo(array1i[0]);
                hhi.assignTo(array1i[1]);
                a2i.assignTo(array1i[2]);
                Matrix newai = new Matrix(array1i).transpose();
                // newa = rotationMatrix * a so rotationMatrix = newa * a^-1
                ai = ai.inverse();
                Matrix matrixi = newai.times(ai);
                double di, bi, ci, gi, hi, fi;
                bi = matrixi.get(0, 1);
                ci = matrixi.get(0, 2);
                di = matrixi.get(1, 0);
                fi = matrixi.get(1, 2);
                gi = matrixi.get(2, 0);
                hi = matrixi.get(2, 1);

                axisi.setX(0, hi - fi);
                axisi.setX(1, ci - gi);
                axisi.setX(2, di - bi);
                double betai = Math.asin(0.5 * Math.sqrt(axisi.squared()));
                if (betai == 0) continue;
                double denominatori = 1 - Math.cos(betai);
                if (denominatori == 0) continue;
                orientationSumi = 1.5 * (betai - Math.sin(betai)) / denominatori ;

                Vector axisj=space.makeVector();
                orientationSumj = orientationSumi ;
                axisi.normalize();
                axisj.E(axisi);


                Vector drjj=space.makeVector();
                Vector drjj2=space.makeVector();

                Vector axisjj=space.makeVector();
                Vector axisjj2=space.makeVector();

                drjj.E(drj);
                drjj2.E(drj);

                axisjj.E(axisj);
                axisjj2.E(axisj);

                D3tt.transform(drjj);
                sum += dri.dot(drjj);

                D3tr.transform(axisjj);
                sum += orientationSumj*dri.dot(axisjj)*Math.sin(betai);

                D3rt.transform(drjj2);
                sum += orientationSumi*axisi.dot(drjj2)*Math.sin(betai);

                D3rr.transform(axisjj2);
                sum += orientationSumj*orientationSumi*axisi.dot(axisjj2)*Math.sin(betai)*Math.sin(betai);    //sum is deltar.phi.delr/theta sum

            }


                }

              x[2]=(fac * temperature) - ( (sum/4)+0.5*( 0.5 * fdotdeltar + orientationSum) ); //THIS ONE WORKS

        x[3]=meterPE.getDataAsScalar(); //potential energy conv
        x[4]=meterPE.getDataAsScalar()*meterPE.getDataAsScalar(); //potential energysquared conv

        x[5]=x[3]*x[4];  //u*u2
        x[6]=x[4]*x[4];

        x[7]=x[0]*x[1];  //uanh*
        x[8]=x[1]*x[1];  //uanh2^2

        x[9]=x[0]*x[2];  //
        x[10]=x[1]*x[2];
        x[11]=x[2]*x[2];

//        System.out.println("x[0]= "+x[0]);
//        System.out.println("x[1]= "+x[1]);
//        System.out.println("x[2]= "+x[2]);
//        System.out.println("x[3]= "+x[3]);
//        System.out.println("x[4]= "+x[4]);

        return data;
    }


    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }
}
