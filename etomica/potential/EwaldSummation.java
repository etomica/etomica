package etomica.potential;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.AtomLeafAgentManager;
import etomica.math.SpecialFunctions;
import etomica.space.ISpace;
import etomica.space.Tensor;
import etomica.space3d.Tensor3D;
import etomica.units.Electron;

/**
 * basic Ewald Sum 
 * U(coulomb) = U(real-space) + U(fourier-space) + U(self-correction)
 * Frenkel&Smit, <Understanding Molecular Simulation from Algorithms to Applications> 
 * p292-300.especially equation(12.1.25)
 * also Allen, computer simulation of liquids Ch5
 * 
 * 
 * @author shu
 * Date:08-09-2012
 */

/**Sabry*/
//Given: rc=L/2  ,  s  ,  
//Calc.: alpha = s/rc    ,  F&S: radius => nc=(s*alpha/PI)*L Sabry: n_max=N  , nx = ny = nz ~ N^1/3   ,
//                            Shu : nx = ny = nz = coefficient_fourier * N^1/3 (coefficient_fourier=4 hERE)

public class EwaldSummation implements PotentialSoft{
    protected final ISpace space;
    protected final AtomLeafAgentManager<MyCharge> atomAgentManager;
    protected final IBox box;
    protected final double alpha, alpha2,alpha3;//sqrt of the Frenkel's alpha, follow other authors' convention
    protected final double boxSize;
    protected final double volume;

    protected final double precision;
    protected final int nKs, nRealShells; 
    protected final IMoleculeList moleculeList;
    protected final int numMolecules;
    protected IVectorMutable[] gradient;
    protected double[] sinkrj, coskrj;
    protected final Tensor secondDerivative;
    protected final Tensor tempTensorkk;
    protected final IVectorMutable rAB;
    protected final IVectorMutable Lxyz;
    protected final IVectorMutable drTmp;
    protected final Tensor identity = new Tensor3D(new double[][] {{1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0}});
    protected final IVectorMutable kVector;
    protected final double rCutSquared;
    protected final double sqrtPI = Math.sqrt(Math.PI);

	// *********************************************** constructor ************************************ // 
    public EwaldSummation(IBox box, AtomLeafAgentManager<MyCharge> atomAgentManager, double precision, ISpace _space, int nKs, double rCut){

        this.box = box;
        this.atomAgentManager = atomAgentManager;
        this.space = _space;
        this.secondDerivative = space.makeTensor();
        this.tempTensorkk = space.makeTensor();
        this.precision = precision; // =???? | SG
        if (precision >= Math.exp(-1)) {
            // you really want better precision than this.
            // also the root finder used below requires this condition
            throw new RuntimeException("precision too large");
        }
        double precision_splus = 0;
        double precision_sminus = Math.sqrt(-Math.log(precision));
        double f = 1;
        double precision_s = 0;
        while (true){
            precision_s = (precision_sminus + precision_splus)/2;
            if(precision_s == precision_sminus || precision_s == precision_splus) break;
            f = Math.exp(-precision_s*precision_s) - precision_s*precision_s*precision;
            if (f==0) break;
            if(f > 0){
                precision_splus = precision_s;
            }
            else {
                precision_sminus = precision_s;
            }
        }
        boxSize = box.getBoundary().getBoxSize().getX(0);
        volume = box.getBoundary().volume();
        Lxyz = space.makeVector();
        rCutSquared = (rCut == 0 ? boxSize*0.49: rCut) ;
//        rCut = boxSize * (0.49 + nRealShells);
        nRealShells = (int) Math.ceil(rCut/boxSize - 0.49);

        alpha = precision_s / rCut;//Separation parameter, obtained from s and L
        alpha2 = alpha*alpha;
        alpha3 = alpha*alpha2;
        this.nKs = nKs;

        moleculeList = box.getMoleculeList();
        numMolecules = moleculeList.getMoleculeCount();
        //nc = (int)Math.ceil(precision_s * boxSize * alpha / Math.PI);//"n = 2s^2/pi" if rC = L/2

        if (false) {
            double precision_sSquared = precision_s * precision_s ;
            double exp_s = Math.exp(-precision_sSquared) / precision_sSquared;
            double e = Electron.UNIT.toSim(1.0); 
            double q_err = numMolecules * e * e ;

            System.out.println("In Ewald sum class");
            System.out.println("exp(-s^2)/s/s:"+exp_s);
            System.out.println("alpha = "+alpha);
            System.out.println("box size is : "+boxSize);
            System.out.println("rCut = "+rCut);
            System.out.println("alpha = s / rCut:"+alpha);
            System.out.println("cf. 5/L:"+5.0/boxSize);
            System.out.println("    3.5/L:"+3.5/boxSize);
            System.out.println("q_err= N*q^2:"+ q_err);

            double basis = 2 * Math.PI / boxSize; // basis(unit) vector magnitude of fourier space vector
            double kCut =  basis * nKs;
            System.out.println("basis vector is :" + basis); 
            System.out.println("kCutoff is :" + kCut); 
        }
        gradient = new IVectorMutable[0];
        sinkrj = new double[0];
        coskrj = new double[0];
        rAB = space.makeVector();
        drTmp = space.makeVector();
        kVector = space.makeVector();
    }

    //////////////////////////////////////////// begin calculating energy //////////////////////////////////////

    // *********************************************************************************************//
    // *************************************  Real-space ******************************************//
    // *********************************************************************************************//
    public double uReal(){
        int nAtoms = box.getLeafList().getAtomCount();

        double uReal = 0.0;
        for (int i=0; i < nAtoms; i++){//H
            IAtom atomA = box.getLeafList().getAtom(i);
            double chargeA = atomAgentManager.getAgent(atomA).charge;
            if (chargeA==0) continue;
            int aIndex = atomA.getParentGroup().getIndex();
            IVectorMutable positionA = atomA.getPosition();
            for (int j=i+1; j < nAtoms; j++){
                IAtom atomB = box.getLeafList().getAtom(j);
                int bIndex = atomB.getParentGroup().getIndex();

                if(aIndex == bIndex && nRealShells==0) continue;//Skip same molecules!

                double chargeB = atomAgentManager.getAgent(atomB).charge;
                if (chargeB==0) continue;

                IVectorMutable positionB = atomB.getPosition();

                rAB.Ev1Mv2(positionA, positionB);// get vector rAB
                box.getBoundary().nearestImage(rAB);// minimum image

//				double rABMagnitudeSquared = rAB.squared() ;// Squared of |rAB|
//				if (rABMagnitudeSquared > rCutSquared) continue ; // check whether the rAB is within the spherical cutoff
//				double rABMagnitude = Math.sqrt(rABMagnitudeSquared);

                for(int nx = -nRealShells; nx <= nRealShells; nx++) {
                    Lxyz.setX(0, nx*box.getBoundary().getBoxSize().getX(0)); 
                    for(int ny = -nRealShells; ny <= nRealShells; ny++) {
                        Lxyz.setX(1, ny*box.getBoundary().getBoxSize().getX(1));
                        for(int nz = -nRealShells; nz <= nRealShells; nz++) {
                            if (aIndex==bIndex && nx*nx+ny*ny+nz*nz == 0) continue;
                            Lxyz.setX(2, nz*box.getBoundary().getBoxSize().getX(2));
                            drTmp.Ev1Pv2(rAB, Lxyz);
                            double r2 = drTmp.squared();
                            if(r2 > rCutSquared) continue;
                            double drTmpM = Math.sqrt(r2);
                            uReal += chargeA * chargeB * SpecialFunctions.erfc( alpha * drTmpM) / drTmpM;//Don't worry about 1/2 factor!
                        }
                    }
                }
            }// close for all sites in j-th molecule
        } // close for the outside loop
        return uReal;
    }

    // *********************************************************************************************//
    // *************************************  Fourier-space ****************************************//
    // *********************************************************************************************//
    public double uFourier(){
//		System.out.println("coefficient_fourier: " + coefficient_fourier); 
//		System.out.println("number of k-vectors along one axis is :"+ nc); 
        // int n =(int)( precision_s * boxSize * alpha / Math.PI);// n(cut): number of vectors along one axis in fourier-space

        double basis = 2 * Math.PI / boxSize; // basis(unit) vector magnitude of fourier space vector
        double kCut =  basis * nKs;
        double kCutSquared = kCut * kCut; // criteria for spherical cutoff in fourier space
        double coefficient = 2.0 / volume * Math.PI ;
        //double uFourier = 0.0;
        double uFourier = 0.0; //>>>>>>>>>>>>>> calculated from cos*cos + sin*sin
        IAtomList atoms = box.getLeafList();
        int nAtoms = atoms.getAtomCount();

        // loop over vectors in k-space. (1)k=(0,0,0) is excluded, (2)within sphere with kCutoff as its radius
        for (int xAxis = -nKs; xAxis < nKs+1; xAxis++){
            kVector.setX(0, (xAxis * basis));// assign value to the x-axis
            for (int yAxis = -nKs; yAxis < nKs+1; yAxis++ ){
                kVector.setX(1, (yAxis * basis));// assign value to the y-axis
                for (int zAxis = -nKs; zAxis < nKs+1; zAxis++ ){
                    if( (xAxis * xAxis + yAxis * yAxis + zAxis * zAxis) == 0) continue;// first check: k is a non-zero vector
                    kVector.setX(2, (zAxis * basis));// assign value to the z-axis, now the vector is specified
                    double kSquared = kVector.squared();

                    if (kSquared > kCutSquared) continue;// k-vector should be within the sphere with kCutoff as the radius

                    double kCoefficientTerm = Math.exp(-0.25 * kSquared / alpha2) / kSquared;//exp(-k*k/4/alpha/alpha) / k/k , a constant for a given kVector
                    double structureFactorReal =  0.0;//>>>>>>>>>>>>> calculated from cos*cos + sin*sin
                    double structureFactorImagine =  0.0;//>>>>>>>>>>>>> calculated from cos*cos + sin*sin
                    for (int i=0; i<nAtoms; i++){
                        IAtom atom = atoms.getAtom(i);
                        IVectorMutable position = atom.getPosition();
                        double charge = atomAgentManager.getAgent(atom).charge;
						if (charge==0) continue;
						double k_r_dot = kVector.dot(position);
						structureFactorReal += charge * Math.cos(k_r_dot);// >>>>>>>>>>>>> calculated from cos*cos + sin*sin
						structureFactorImagine += charge * Math.sin(k_r_dot);// >>>>>>>>>>>>> calculated from cos*cos + sin*sin
                    }// close loop for all sites within 1 molecule
                    double structureFactorSquared = structureFactorReal * structureFactorReal + structureFactorImagine * structureFactorImagine ;//>>>>>>>>>>>>>calculated from cos*cos + sin*sin
                    uFourier +=  kCoefficientTerm * structureFactorSquared; //>>>>>>>>>>>>>calculated from cos*cos + sin*sin
                }// close for z-axis
            }//close for y-axis
        }// close for x-axis(all non-zero kVecors)
        double u = coefficient * uFourier; 

        return u; 
    }

    // *********************************************************************************************//
    // ********************** self-correction Part************************************************* // 
    // *********************************************************************************************//
    public double uSelf(){
        double uSelf = 0.0;

        IAtomList atoms = box.getLeafList();
        int nAtoms = atoms.getAtomCount();
        for (int i=0; i<nAtoms; i++){
            IAtom atom = atoms.getAtom(i);
            double charge = atomAgentManager.getAgent(atom).charge;
            uSelf += charge*charge;
        }

        uSelf *= -alpha/sqrtPI;
        return uSelf;
    }

    public double uBondCorr(){
        double uCorr = 0.0;
        for (int i=0; i< numMolecules; i++){
            IMolecule molecule = moleculeList.getMolecule(i);
            int numSites = molecule.getChildList().getAtomCount();
            for (int siteA=0; siteA<numSites; siteA++){
                IAtom atomA = molecule.getChildList().getAtom(siteA);
                IVectorMutable positionA = atomA.getPosition();
                double chargeA = atomAgentManager.getAgent(atomA).charge;
                if (chargeA==0) continue;
                for (int siteB=siteA+1; siteB<numSites; siteB++){
                    IAtom atomB = molecule.getChildList().getAtom(siteB);
                    IVectorMutable positionB = atomB.getPosition();
                    double chargeB = atomAgentManager.getAgent(atomB).charge;
                    if (chargeB==0) continue;

                    rAB.Ev1Mv2(positionA, positionB);
                    box.getBoundary().nearestImage(rAB);
                    double rABMagnitudeSquared = rAB.squared();
                    double rABMagnitude = Math.sqrt(rABMagnitudeSquared);

                    uCorr += chargeA*chargeB*(1-SpecialFunctions.erfc(alpha*rABMagnitude))/rABMagnitude;
                }
            }
        }		
        return uCorr;
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public int nBody() {
        return 0;
    }

    public void setBox(IBox box) {
    }

    //******************************** inner class ********************************************//
    public static class MyCharge{
        public MyCharge(double charge){
            this.charge = charge;
        }
        public final double charge;
    }

    public double energy(IAtomList atoms) {
        double real = uReal();
        double fourier = uFourier();
        double self = uSelf();
        double bondCorr = uBondCorr();

        double totalEnergy = real + fourier + self + bondCorr;

        if (false) { 
            System.out.println("total:               "+ totalEnergy/numMolecules);
            System.out.println("real   : "+ real/numMolecules);
            System.out.println("fourier: "+ fourier/numMolecules);
            System.out.println("self: "+ self/numMolecules);
            System.out.println("bond correction: "+ bondCorr/numMolecules);
        }

        return totalEnergy;
    }

    public double virial(IAtomList atoms) {
        return 0;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    //////////////////////////////////////////// begin calculating gradient //////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public IVector[] gradient(IAtomList atoms) {
        int nAtoms = box.getLeafList().getAtomCount();
        double coeff = 4.0*Math.PI / box.getBoundary().volume();
        double basis = 2 * Math.PI / boxSize; // basis(unit) vector magnitude of fourier space vector
        double kCut =  basis * nKs;
        double kCutSquared = kCut * kCut; // criteria for spherical cutoff in fourier space

        if(gradient.length < nAtoms){
            gradient = new IVectorMutable[nAtoms];
            sinkrj = new double[nAtoms];
            coskrj = new double[nAtoms];
            for(int i=0; i < nAtoms; i++){
                gradient[i] = space.makeVector();
            }
        }else{
            for(int i=0; i < nAtoms; i++){
                gradient[i].E(0); //for Vecors and scalar
                sinkrj[i] = 0;
                coskrj[i] = 0;
            }
        }

        //Real gradient  //Cross Interaction
        for (int i=0; i < nAtoms; i++){
            IAtom atomA = box.getLeafList().getAtom(i);
            double chargeA = atomAgentManager.getAgent(atomA).charge;
            if (chargeA==0) continue;
            int aIndex = atomA.getParentGroup().getIndex(); // molecule a
            IVectorMutable positionA = atomA.getPosition();
            for (int j=i+1; j < nAtoms; j++){
                IAtom atomB = box.getLeafList().getAtom(j);
                int bIndex = atomB.getParentGroup().getIndex(); // molecule b

                if(nRealShells == 0 && aIndex == bIndex) continue;//Skip same molecules!

                double chargeB = atomAgentManager.getAgent(atomB).charge;
                if (chargeB==0) continue;
                IVectorMutable positionB = atomB.getPosition();
                rAB.Ev1Mv2(positionA, positionB); //rAB == rA - rB
                box.getBoundary().nearestImage(rAB);
                for (int nx = -nRealShells; nx <= nRealShells; nx++) {
                    Lxyz.setX(0, nx*box.getBoundary().getBoxSize().getX(0)); 
                    for (int ny = -nRealShells; ny <= nRealShells; ny++) {
                        Lxyz.setX(1, ny*box.getBoundary().getBoxSize().getX(1));
                        for (int nz = -nRealShells; nz <= nRealShells; nz++) {
                            if (aIndex==bIndex && nx*nx+ny*ny+nz*nz == 0) continue;
                            Lxyz.setX(2, nz*box.getBoundary().getBoxSize().getX(2));
                            drTmp.Ev1Pv2(rAB, Lxyz);

                            double rAB2 = drTmp.squared();
                            if (rAB2 > rCutSquared) continue; 
                            double rABMagnitude = Math.sqrt(rAB2);
                            double rAB3 = rABMagnitude*rAB2;
                            double B = SpecialFunctions.erfc(alpha*rABMagnitude) + 2.0*alpha*rABMagnitude/sqrtPI * Math.exp(-alpha2*rAB2) ;
                            double realCoeff = - chargeA*chargeB * B / rAB3; // gradU = -F
                            gradient[i].PEa1Tv1(realCoeff, drTmp);
                            gradient[j].PEa1Tv1(-realCoeff, drTmp);
                        }
                    }
                }
            }
        }

        //Fourier gradient Part		

        for (int xAxis = -nKs; xAxis < nKs+1; xAxis++){
            kVector.setX(0, (xAxis * basis));// assign value to the x-axis
            for (int yAxis = -nKs; yAxis < nKs+1; yAxis++ ){
                kVector.setX(1, (yAxis * basis));// assign value to the y-axis
                for (int zAxis = -nKs; zAxis < nKs+1; zAxis++ ){
                    if( (xAxis * xAxis + yAxis * yAxis + zAxis * zAxis) == 0) continue;// first check: k is a non-zero vector
                    kVector.setX(2, (zAxis * basis));// assign value to the z-axis, now the vector is specified
                    double kSquared = kVector.squared();
                    if (kSquared > kCutSquared) continue;// k-vector should be within the sphere with kCutoff as the radius
                    double sCoskr = 0.0, sSinkr = 0.0;

                    for (int j=0; j< nAtoms ; j++){ // Loop over atoms (4*nBasis)
                        IAtom atom = box.getLeafList().getAtom(j);
                        IVectorMutable position = atom.getPosition();
                        sinkrj[j] = Math.sin(position.dot(kVector));
                        coskrj[j] = Math.cos(position.dot(kVector));
                        double chargej = atomAgentManager.getAgent(atom).charge;

                        sCoskr += chargej*coskrj[j];
                        sSinkr += chargej*sinkrj[j]; 
                    }//End loop over j

                    double coeffk = coeff / kSquared * Math.exp(-kSquared/4.0/alpha2);
                    for(int i=0; i<nAtoms; i++){
                        IAtom atom = box.getLeafList().getAtom(i);
                        double chargei = atomAgentManager.getAgent(atom).charge;
                        double coeffki = coeffk * chargei * (sinkrj[i] * sCoskr - coskrj[i] * sSinkr); 
                        gradient[i].PEa1Tv1(-coeffki , kVector);  // gradU = -F
                    }
//		            System.out.println("Fourier = " + gradient[0]);	
                }//end of storing Sin and Cos
            }
        }//End loop over ks

        //Intra-Molecular  gradient:
        for (int i=0; i< numMolecules; i++){
            IMolecule molecule = moleculeList.getMolecule(i);	
            int numSites = molecule.getChildList().getAtomCount();
            for (int siteA=0; siteA<numSites; siteA++){
                IAtom atomA = molecule.getChildList().getAtom(siteA); // index = 0, 1, 2, 3|||leafIndex=0...184
                double chargeA = atomAgentManager.getAgent(atomA).charge;
                if (chargeA==0) continue;
                IVectorMutable positionA = atomA.getPosition();
                for (int siteB=siteA+1; siteB<numSites; siteB++){
                    IAtom atomB = molecule.getChildList().getAtom(siteB);
                    double chargeB = atomAgentManager.getAgent(atomB).charge;
                    if (chargeB==0) continue;
                    IVectorMutable positionB = atomB.getPosition();

                    rAB.Ev1Mv2(positionA, positionB);
                    box.getBoundary().nearestImage(rAB);
                    double rAB2 = rAB.squared();
                    double rABMagnitude = Math.sqrt(rAB2);
                    // U = Ur + Uf - Uself - U_intra ====> U_intra = Erf(alpha r)/r
                    // dU = -F = Ur' + Uf' - Uself' - U_intra' ===>   -d[Erf(alpha r)/r]/dx
                    double B = 2*alpha/sqrtPI * Math.exp(-alpha2*rAB2)-(1-SpecialFunctions.erfc(alpha*rABMagnitude))/rABMagnitude; 
                    double coeffAB = - chargeA*chargeB * B / rAB2; // gradU = -F
                    gradient[atomA.getLeafIndex()].PEa1Tv1(coeffAB, rAB);
                    gradient[atomB.getLeafIndex()].PEa1Tv1(-coeffAB, rAB);
                }
            }
        }
        return gradient;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    //////////////////////////          begin calculating secondDerivative               /////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public Tensor secondDerivative(IAtom atom0, IAtom atom1){
        double basis = 2 * Math.PI / boxSize; // basis(unit) vector magnitude of fourier space vector
        double kCut =  basis * nKs;
        double kCutSquared = kCut * kCut; // criteria for spherical cutoff in fourier space
        IVectorMutable pos0 = atom0.getPosition();
        IVectorMutable pos1 = atom1.getPosition();
        rAB.Ev1Mv2(pos0,pos1);
        box.getBoundary().nearestImage(rAB);

        double q0 = atomAgentManager.getAgent(atom0).charge;
        double q1 = atomAgentManager.getAgent(atom1).charge;
        secondDerivative.E(0);
        if (q0*q1 == 0) return secondDerivative;
        double coeff = 4.0*Math.PI/volume;

        //Fourier Space ... Long || ALL pairs of atoms in all unit cells (including atoms in same mol.)		
        for (int xAxis = -nKs; xAxis < nKs+1; xAxis++){
            kVector.setX(0, (xAxis * basis));// assign value to the x-axis
            for (int yAxis = -nKs; yAxis < nKs+1; yAxis++ ){
                kVector.setX(1, (yAxis * basis));// assign value to the y-axis
                for (int zAxis = -nKs; zAxis < nKs+1; zAxis++ ){
                    if( (xAxis * xAxis + yAxis * yAxis + zAxis * zAxis) == 0) continue;// first check: k is a non-zero vector
                    kVector.setX(2, (zAxis * basis));// assign value to the z-axis, now the vector is specified
                    double kSquared = kVector.squared();
                    if (kSquared > kCutSquared) continue;// k-vector should be within the sphere with kCutoff as the radius
                    tempTensorkk.Ev1v2(kVector, kVector);
                    tempTensorkk.TE(Math.exp(-kSquared/4.0/alpha2) * Math.cos(kVector.dot(rAB))/kSquared);
                    secondDerivative.PE(tempTensorkk);
                }
            }
        }
        secondDerivative.TE(coeff);
        //Rreal.Short..diff mol.
        boolean intraMolecular = atom0.getParentGroup() == atom1.getParentGroup();

        for(int nx = -nRealShells; nx <= nRealShells; nx++) {
            Lxyz.setX(0, nx*box.getBoundary().getBoxSize().getX(0)); 
            for(int ny = -nRealShells; ny <= nRealShells; ny++) {
                Lxyz.setX(1, ny*box.getBoundary().getBoxSize().getX(1));
                for(int nz = -nRealShells; nz <= nRealShells; nz++) {
                    Lxyz.setX(2, nz*box.getBoundary().getBoxSize().getX(2));
                    drTmp.Ev1Pv2(rAB, Lxyz);

                    double rAB2 = drTmp.squared();
                    // if we're beyond the cutoff, then skip unless we're 
                    if (rAB2 > rCutSquared && (!intraMolecular || (nx*nx+ny*ny+nz*nz > 0))) continue;
                    double rABMagnitude = Math.sqrt(rAB2);
                    double rAB3 = rABMagnitude*rAB2;
                    double rAB4 = rAB2*rAB2;
                    double rAB5 = rABMagnitude*rAB4;
                    double erfc = SpecialFunctions.erfc(alpha*rABMagnitude);
                    double exp_Alpha2r2 = Math.exp(-alpha2*rAB2);
                    double B0 = (6*alpha/rAB4 + 4*alpha3/rAB2)/sqrtPI * exp_Alpha2r2;

                    tempTensorkk.Ev1v2(drTmp, drTmp);

                    if (intraMolecular && nx*nx+ny*ny+nz*nz==0) {
//                      D = D(real) + D(rec.) - D(Intra)
                        double A = (1-erfc)/rAB3 - 2*alpha/sqrtPI*exp_Alpha2r2/rAB2;
                        double B = B0 - 3 * (1-erfc)/rAB5;
                        tempTensorkk.TE(-B);
                        tempTensorkk.PEa1Tt1(A,identity);
                    }
                    else {
                        double A = erfc/rAB3 + 2*alpha/sqrtPI*exp_Alpha2r2/rAB2;
                        double B = B0 + 3*erfc/rAB5;
                        tempTensorkk.TE(B);
                        tempTensorkk.PEa1Tt1(A,identity);
                    }

                    secondDerivative.PE(tempTensorkk);
                }
            }
        }

        secondDerivative.TE(q0*q1);
        return secondDerivative;
    }

    public IVector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return null;
    }

    public class P2EwaldReal implements PotentialSoft {

        protected final IVectorMutable[] gradient2;

        public P2EwaldReal(ISpace space) {
            gradient2 = new IVectorMutable[2];
            gradient2[0] = space.makeVector();
            gradient2[1] = space.makeVector();
        }

        public double energy(IAtomList atoms) {
            IAtom atomA = atoms.getAtom(0);
            double chargeA = atomAgentManager.getAgent(atomA).charge;
            if (chargeA==0) return 0;

            IAtom atomB = atoms.getAtom(1);
            double chargeB = atomAgentManager.getAgent(atomB).charge;
            if (chargeB==0) return 0;

            IVectorMutable positionA = atomA.getPosition();
            IVectorMutable positionB = atomB.getPosition();

            rAB.Ev1Mv2(positionA, positionB);// get vector rAB
            box.getBoundary().nearestImage(rAB);// minimum image

            double r2 = rAB.squared();
            if(r2 > rCutSquared) return 0;
            double r = Math.sqrt(r2);
            return chargeA * chargeB * SpecialFunctions.erfc( alpha * r) / r;//Don't worry about 1/2 factor!
        }

        public double getRange() {
            return Math.sqrt(rCutSquared);
        }

        public void setBox(IBox box) {}

        public int nBody() {
            return 2;
        }

        public double virial(IAtomList atoms) {
            return 0;
        }

        public IVector[] gradient(IAtomList atoms) {
            //Real gradient  //Cross Interaction
            IAtom atomA = atoms.getAtom(0);
            double chargeA = atomAgentManager.getAgent(atomA).charge;
            if (chargeA==0) {
                gradient2[0].E(0);
                gradient2[1].E(0);
                return gradient2;
            }
            IVectorMutable positionA = atomA.getPosition();
            IAtom atomB = atoms.getAtom(1);

            double chargeB = atomAgentManager.getAgent(atomB).charge;
            if (chargeB==0) {
                gradient2[0].E(0);
                gradient2[1].E(0);
                return gradient2;
            }

            IVectorMutable positionB = atomB.getPosition();
            rAB.Ev1Mv2(positionA, positionB); //rAB == rA - rB
            box.getBoundary().nearestImage(rAB);

            double rAB2 = rAB.squared();
            if (rAB2 > rCutSquared) {
                gradient2[0].E(0);
                gradient2[1].E(0);
                return gradient2;
            }

            double rABMagnitude = Math.sqrt(rAB2);
            double rAB3 = rABMagnitude*rAB2;
            double B = SpecialFunctions.erfc(alpha*rABMagnitude) + 2.0*alpha*rABMagnitude/sqrtPI * Math.exp(-alpha2*rAB2) ;
            double realCoeff = - chargeA*chargeB * B / rAB3; // gradU = -F
            gradient2[0].PEa1Tv1(realCoeff, rAB);
            gradient2[1].PEa1Tv1(-realCoeff, rAB);
            return gradient2;
        }

        public IVector[] gradient(IAtomList atoms, Tensor pressureTensor) {
            return gradient(atoms);
        }
    }
}
