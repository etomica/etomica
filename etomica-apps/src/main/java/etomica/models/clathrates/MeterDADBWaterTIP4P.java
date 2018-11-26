package etomica.models.clathrates;
import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeAgentManager;
import etomica.molecule.MoleculePositionCOM;
import etomica.normalmode.AtomSiteSource;
import etomica.normalmode.LatticeSumMolecularCrystal;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Orientation3D;
import etomica.space3d.OrientationFull3D;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Tensor3D;
import etomica.spaceNd.TensorND;
import etomica.units.dimensions.Null;

import java.util.Arrays;
public class MeterDADBWaterTIP4P  implements IDataSource {


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
    protected double latticeEnergy;
    protected final double temperature;
    public static boolean justDADB = true;
    public static boolean justU = false;
    protected final Vector ph1h2;
    protected final Vector q;
    protected final Vector totalForce;
    private final Vector centerMass;
    public boolean doTranslation;
    public boolean doRotation;
    protected   Tensor tmpDrr1;
    protected final Tensor[] tmpAtomicTensor3;//46X4=184 atoms dimentional array of 3dim Tensor

    protected final Vector torque;

    public MeterDADBWaterTIP4P(Space space, DataSourceScalar meterPE, PotentialMaster potentialMaster, double temperature, MoleculeAgentManager latticeCoordinates) {
        int nData = justDADB ? 1 : 9;
        data = new DataDoubleArray(nData);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.space = space;
        this.latticeCoordinates = latticeCoordinates;
        this.meterPE = meterPE;
        this.potentialMaster = potentialMaster;
        id = new IteratorDirective();
        pcForceSum = new PotentialCalculationForceSum();
        forceManager = new AtomLeafAgentManager<>(a -> space.makeVector(), latticeCoordinates.getBox());
        pcForceSum.setAgentManager(forceManager);
        dr = space.makeVector();
        MeterPotentialEnergy meterPE2 = new MeterPotentialEnergy(potentialMaster);
        meterPE2.setBox(latticeCoordinates.getBox());
        latticeEnergy = meterPE2.getDataAsScalar();
        this.temperature = temperature;
        ph1h2 = space.makeVector();
        q = space.makeVector();
        totalForce = space.makeVector();
        centerMass = space.makeVector();
        torque = space.makeVector();
        this.tmpAtomicTensor3 = new Tensor[basisDim*atomsPerMol];//46X4=184 atoms dimensional array of 3D Tensor
        for (int i=0; i<basisDim*atomsPerMol ;i++){
            tmpAtomicTensor3[i] = space.makeTensor();
        }
    }

    public void setLatticeEnergy(double newLatticeEnergy) {
        latticeEnergy = newLatticeEnergy;
    }


    //aTensor is atomicTensorAtomicPair
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
            if (atom0.getType() == sim.species.getOxygenType() && atom1.getType() == sim.species.getOxygenType()) {
                dr.Ev1Mv2(atom0.getPosition(), atom1.getPosition());
                sim.box.getBoundary().nearestImage(dr);

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
                            double dW = sim.potentialLJ.du(dr2Tmp);
                            double d2W = sim.potentialLJ.d2u(dr2Tmp);
                            tmpD3LJ.TE(1.0 / (dr2Tmp * dr2Tmp) * (dW - d2W));
                            tmpD3LJ.PEa1Tt1(-dW / dr2Tmp, identity);
                            D3ESLJ.PE(tmpD3LJ);
                        }
                    }
                }
                //ES
            } else if (atom0.getType() != sim.species.getOxygenType() && atom1.getType() != sim.species.getOxygenType()) {
                D3ESLJ.PE(sim.potentialES.secondDerivative(atom0, atom1));
            }
            return D3ESLJ;
        }
    }; // End AtomicTensorAtomicPair

    public IData getData() {
        Box box = latticeCoordinates.getBox();
        double[] x = data.getData();
        double x0 = meterPE.getDataAsScalar() - latticeEnergy;
        pcForceSum.reset();
        potentialMaster.calculate(box, id, pcForceSum);
        IMoleculeList molecules = box.getMoleculeList();
        double ForceSum = 0;
        double orientationSum = 0;
 //       AtomLeafAgentManager<Vector> atomLatticeCoordinates = new AtomLeafAgentManager<>(new AtomSiteSource(space), box, Vector.class);

        for (int i = 0; i < molecules.size(); i++) {    //#########GET GETATOM/MOLECULE WORKING
            IMolecule molecule = molecules.get(i);
            IAtomList leafList = molecule.getChildList();

            Vector h1Force = forceManager.getAgent(leafList.get(0));
            Vector h2Force = forceManager.getAgent(leafList.get(1));
            Vector oForce = forceManager.getAgent(leafList.get(2));
            Vector mForce = forceManager.getAgent(leafList.get(3));

            totalForce.E(h1Force); //########################WHAT IS E,PE,TE?
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
            centerMass.PEa1Tv1(oMass, o);   //WHAT IS THIS PE AND TE
            centerMass.TE(1 / (2 * hMass + oMass));   //##################WHAT IS TE?

            dr.Ev1Mv2(m, centerMass);
            torque.E(mForce);
            torque.XE(dr);  //#############################WHAT DOES THIS XE DO?
            q.E(torque);
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
            q.PE(torque);
            //for the total torque q
             Vector lPos = ((MoleculeSiteSource.LatticeCoordinate) latticeCoordinates.getAgent(molecule)).position;
            dr.Ev1Mv2(centerMass, lPos);
            ForceSum += totalForce.dot(dr);

//bring class
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
//            orientationSum += 0.5 * betaNew * DUDT;  //Only kappa3


        }

        if (!doTranslation) ForceSum = 0;
        if (!doRotation) orientationSum = 0;
        if (justDADB) {
            if (justU) {
                int N = molecules.size();
                double fac = (doTranslation ? 1.5 : 0) * (N - 1) + (doRotation ? 1.5 : 0) * N;
                x[0] = (x0 + latticeEnergy) + (fac * temperature) + 0.5 * ForceSum + orientationSum;

                for (int i = 0; i < molecules.size(); i++) {    //#########GET GETATOM/MOLECULE WORKING
                    for (int j = 0; j > i; j++) {  //12, 13...


                        tmpDrr1 = space.makeTensor();
                            Tensor3D D3tt  = new Tensor3D();	Tensor3D D3tr  = new Tensor3D();
                            Tensor3D D3rt  = new Tensor3D();	Tensor3D D3rr  = new Tensor3D();
                            Tensor3D D3tt_ = new Tensor3D(); Tensor3D D3tr_ = new Tensor3D();
                            Tensor3D Rk = new Tensor3D();  Tensor3D Rk_ = new Tensor3D();
                            Tensor3D Rkp = new Tensor3D();
                            Vector Xk = space.makeVector();
                            Vector Xkp = space.makeVector();
                            MoleculePositionCOM com_0 = new MoleculePositionCOM(space);
                            Vector com0 = com_0.position(molecules.get(i));
                            MoleculePositionCOM com_1 = new MoleculePositionCOM(space);
                            Vector com1 = com_1.position(molecules.get(j));

//    	com0.E(mol0.getChildList().getAtom(2).getPosition()); // O (-5.970371160466783, 5.978273273935142, -2.996126942837739)
//    	com1.E(mol1.getChildList().getAtom(2).getPosition()); // O (-6.016203213551466, 6.025148464416224, -2.996521341193713)
//    	com0.PE(4.4);
//    	com1.PE(4.4);
                            int numSites0 = molecules.get(i).getChildList().size();
                            int numSites1 = molecules.get(j).getChildList().size();
                            for (int atomk=0; atomk < numSites0; atomk++){
                                Vector posk = molecules.get(i).getChildList().get(atomk).getPosition();
                                Xk.Ev1Mv2(posk, com0);
                                box.getBoundary().nearestImage(Xk);
                                Rk.setComponent(0,0,0.0); Rk.setComponent(1,1,0.0);	Rk.setComponent(2,2,0.0);
                                Rk.setComponent(0,1, Xk.getX(2));  Rk.setComponent(1,0,-Xk.getX(2));
                                Rk.setComponent(0,2,-Xk.getX(1));  Rk.setComponent(2,0, Xk.getX(1));
                                Rk.setComponent(1,2, Xk.getX(0));  Rk.setComponent(2,1,-Xk.getX(0));
                                for (int atomkp=0; atomkp < numSites1; atomkp++){
                                    if(atomk == atomkp && molecules.get(i) == molecules.get(j)) continue;//ADDED:: Non-self
                                    Vector poskp = molecules.get(j).getChildList().get(atomkp).getPosition();
                                    Xkp.Ev1Mv2(poskp, com1);
                                    box.getBoundary().nearestImage(Xkp);
                                    Rkp.setComponent(0,0,0.0); Rkp.setComponent(1,1,0.0);	Rkp.setComponent(2,2,0.0);
                                    Rkp.setComponent(0,1, Xkp.getX(2));  Rkp.setComponent(1,0,-Xkp.getX(2));
                                    Rkp.setComponent(0,2,-Xkp.getX(1));  Rkp.setComponent(2,0, Xkp.getX(1));
                                    Rkp.setComponent(1,2, Xkp.getX(0));  Rkp.setComponent(2,1,-Xkp.getX(0));

                                    D3tt_.E(aTensor.atomicTensor(molecules.get(i).getChildList().get(atomk) , molecules.get(j).getChildList().get(atomkp)));
                                    D3tt.PE(D3tt_);
                                    D3tr_.E(D3tt_); D3tr_.TE(Rkp);  D3tr.PE(D3tr_);
                                    Rk_.E(Rk);  Rk_.TE(D3tt_);  D3rt.ME(Rk_);
                                    Rk_.TE(Rkp); D3rr.ME(Rk_);
                                    tmpAtomicTensor3[molecules.get(i).getChildList().get(atomk).getLeafIndex()].ME(D3tt_);//self summation(sum rule)
                                }//atomkp

                                if(molecules.get(i) == molecules.get(j)){//self
                                    D3tt_.E(tmpAtomicTensor3[molecules.get(i).getChildList().get(atomk).getLeafIndex()]);

                                    D3tt.PE(D3tt_);//Transform to molecular

                                    D3tr_.E(D3tt_);
                                    D3tr_.TE(Rk);
                                    D3tr.PE(D3tr_);

                                    Rk_.E(Rk);
                                    Rk_.TE(D3tt_); //D3tt_=  - SUM(k =/= k')
                                    D3rt.ME(Rk_);
//Drr Original
                                    Rk_.TE(Rk);
                                    D3rr.ME(Rk_);
//Drr Symmetric Part

                                    Tensor3D tmpDrr2 = new Tensor3D();
                                    IAtom atom = molecules.get(i).getChildList().get(atomk);

                                    Vector fk = forceManager.getAgent(atom);//gradient NOT fk
                                    Xk.TE(-1);
                                    tmpDrr1.Ev1v2(Xk, fk);
                                    if(true){ //Symmetrize? SUM_over_atoms{tmpDrr1} should = ZERO; i.e. torque=0
                                        tmpDrr2.E(tmpDrr1);
                                        tmpDrr1.transpose();
                                        tmpDrr1.PE(tmpDrr2);
                                        tmpDrr1.TE(0.5);
                                    }
                                    tmpDrr1.setComponent(0, 0,  -fk.getX(1) * Xk.getX(1) - fk.getX(2) * Xk.getX(2));
                                    tmpDrr1.setComponent(1, 1,  -fk.getX(0) * Xk.getX(0) - fk.getX(2) * Xk.getX(2));
                                    tmpDrr1.setComponent(2, 2,  -fk.getX(0) * Xk.getX(0) - fk.getX(1) * Xk.getX(1));
                                    D3rr.PE(tmpDrr1);
                                }//mol0 == mol1 ?
                            }//atomk




                    }
                }

                for (int i = 0; i < molecules.size(); i++) {    //#########GET GETATOM/MOLECULE WORKING
                    LatticeSumMolecularCrystal.atomicToMolecularD(aTensor, i, i);  //11, 22, etc. (call same molecules)

                    COPY ABOVE STUFF HERE FOR I=J

                }

            } else {
                x[0] = x0 + 0.5 * ForceSum + orientationSum; //translation and rotation
            }

            if (data.isNaN()) {
                throw new RuntimeException();
            }
            return data;
        }


        return data;
    }




    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }










}
