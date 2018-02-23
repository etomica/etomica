/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Tensor3D;
import etomica.units.Joule;
import etomica.util.Constants;
import org.apache.commons.math3.special.Erf;

/**Sabry*/
//Given: rc=L/2  ,  s  ,  
//Calc.: alpha = s/rc    ,  F&S: radius => nc=(s*alpha/PI)*L Sabry: n_max=N  , nx = ny = nz ~ N^1/3   ,
//                            Shu : nx = ny = nz = coefficient_fourier * N^1/3 (coefficient_fourier=4 hERE)

public class EwaldSummation implements PotentialSoft{
    protected final Space space;
    protected final AtomLeafAgentManager<MyCharge> atomAgentManager;
    protected final Box box;
    protected double alpha, alpha2,alpha3;//sqrt of the Frenkel's alpha, follow other authors' convention
    protected final double[] boxSize, basis;
    protected final double volume;
    protected final int[] nKs, nRealShells; 
    protected final IMoleculeList moleculeList;
    protected Vector[] gradient;
    protected double[] sinkrj, coskrj;
    protected final Tensor secondDerivative;
    protected final Tensor tempTensorkk;
    protected final Vector rAB;
    protected final Vector Lxyz;
    protected final Vector drTmp;
    protected final Tensor identity = new Tensor3D(new double[][] {{1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0}});
    protected final Vector kVector;
    protected double rCutRealES, rCutSquared, kCut;
    protected final double sqrtPI = Math.sqrt(Math.PI);
    protected boolean doRealSum = true;

	// *********************************************** constructor ************************************ // 
    public EwaldSummation(Box box, AtomLeafAgentManager<MyCharge> atomAgentManager, Space _space, double kCut, double rCutRealES){

        this.box = box;
        this.atomAgentManager = atomAgentManager;
        this.space = _space;
        this.secondDerivative = space.makeTensor();
        this.tempTensorkk = space.makeTensor();
        this.kCut = kCut;

        moleculeList = box.getMoleculeList();
        boxSize = new double[] {box.getBoundary().getBoxSize().getX(0),  box.getBoundary().getBoxSize().getX(1),  box.getBoundary().getBoxSize().getX(2)};
        volume = box.getBoundary().volume();
        Lxyz = space.makeVector();
        this.rCutRealES = rCutRealES;
        rCutSquared = rCutRealES*rCutRealES;
        nRealShells = new int[] {(int) Math.ceil(rCutRealES/boxSize[0] - 0.49999), (int) Math.ceil(rCutRealES/boxSize[1] - 0.49999), (int) Math.ceil(rCutRealES/boxSize[2] - 0.49999)};
        
		nKs = new int[] {(int) Math.ceil(boxSize[0]/2.0/Math.PI*kCut) , (int) Math.ceil(boxSize[1]/2.0/Math.PI*kCut), (int) Math.ceil(boxSize[2]/2.0/Math.PI*kCut)};

		double s = Math.sqrt(rCutRealES*kCut/2);
		
		alpha = s/rCutRealES; //=0.2406189232882774  //nX=1
		alpha2 = alpha*alpha;
        alpha3 = alpha*alpha2;

        basis = new double[]{2*Math.PI/boxSize[0], 2*Math.PI/boxSize[1], 2*Math.PI/boxSize[2]};
        gradient = new Vector[0];
        sinkrj = new double[0];
        coskrj = new double[0];
        rAB = space.makeVector();
        drTmp = space.makeVector();
        kVector = space.makeVector();
    }

    /**
     * Sets a new value of rCutRealES
     * This will not alter values of alpha or nRealShells!
     */
    public void setRCut(double newRCutRealES) {
        rCutRealES = newRCutRealES;
        rCutSquared = rCutRealES*rCutRealES;
    }

    /**
     * Returns real-space cutoff
     */
    public double getRCut() {
        return rCutRealES;
    }

    //////////////////////////////////////////// begin calculating energy //////////////////////////////////////

    // *********************************************************************************************//
    // *************************************  Real-space ******************************************//
    // *********************************************************************************************//
    public double uReal(){
        int nAtoms = box.getLeafList().size();
        double uReal = 0.0;
        for (int i=0; i < nAtoms; i++){//H
            IAtom atomA = box.getLeafList().get(i);
            double chargeA = atomAgentManager.getAgent(atomA).charge;

            if (chargeA==0) continue;
 
            int aIndex = atomA.getParentGroup().getIndex();
            Vector positionA = atomA.getPosition();
            for (int j=i; j < nAtoms; j++){
                IAtom atomB = box.getLeafList().get(j);
                int bIndex = atomB.getParentGroup().getIndex();
     
                if(aIndex == bIndex && nRealShells[0]==0 && nRealShells[1]==0 && nRealShells[2]==0) continue;//Skip atom-pairs in the same molecule in the orig. cell.
      
                double chargeB = atomAgentManager.getAgent(atomB).charge;
              
                if (chargeB==0) continue;

                Vector positionB = atomB.getPosition();
                rAB.Ev1Mv2(positionA, positionB);// get vector rAB
                box.getBoundary().nearestImage(rAB);// minimum image
            	boolean isSelf = i == j;
                for(int nx = -nRealShells[0]; nx <= nRealShells[0]; nx++) {
                    Lxyz.setX(0, nx*boxSize[0]); 
                    for(int ny = -nRealShells[1]; ny <= nRealShells[1]; ny++) {
                        Lxyz.setX(1, ny*boxSize[1]);
                        for(int nz = -nRealShells[2]; nz <= nRealShells[2]; nz++) {
                        	boolean centerImage = nx*nx+ny*ny+nz*nz == 0;
                        	
                            if (aIndex==bIndex && centerImage) continue;//Skip atom-pairs in the same molecule in the orig. cell & ignores self+centerImage too
                            
                            Lxyz.setX(2, nz*boxSize[2]);
                            drTmp.Ev1Pv2(rAB, Lxyz);
                            double r2 = drTmp.squared();
                            if(r2 > rCutSquared) continue;
                            double drTmpM = Math.sqrt(r2);
                            double tmepReal = chargeA * chargeB * Erf.erfc(alpha * drTmpM) / drTmpM; //Don't worry about 1/2 factor;j>i
                            uReal+= (isSelf ? 0.5 : 1.0)*tmepReal;
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
        double kCutSquared = kCut * kCut; // criteria for spherical cutoff in fourier space
        double coefficient = 2.0*Math.PI/volume;
        double uFourier = 0.0; //>>>>>>>>>>>>>> calculated from cos*cos + sin*sin
        IAtomList atoms = box.getLeafList();
        int nAtoms = atoms.size();

        // loop over vectors in k-space. (1)k=(0,0,0) is excluded, (2)within sphere with kCutoff as its radius
        for (int kX = -nKs[0]; kX < nKs[0]+1; kX++){
            kVector.setX(0, (kX * basis[0]));// assign value to the x-axis
            for (int kY = -nKs[1]; kY < nKs[1]+1; kY++ ){
                kVector.setX(1, (kY * basis[1]));// assign value to the y-axis
                for (int kZ = -nKs[2]; kZ < nKs[2]+1; kZ++ ){
                    if( (kX * kX + kY * kY + kZ * kZ) == 0) continue;// first check: k is a non-zero vector
                    kVector.setX(2, (kZ * basis[2]));// assign value to the z-axis, now the vector is specified
                    double kSquared = kVector.squared();

                    if (kSquared > kCutSquared) continue;// k-vector should be within the sphere with kCutoff as the radius

                    double kCoefficientTerm = Math.exp(-0.25 * kSquared / alpha2) / kSquared;//exp(-k*k/4/alpha/alpha) / k/k , a constant for a given kVector
                    double structureFactorReal =  0.0;//>>>>>>>>>>>>> calculated from cos*cos + sin*sin
                    double structureFactorImagine =  0.0;//>>>>>>>>>>>>> calculated from cos*cos + sin*sin
                    for (int i=0; i<nAtoms; i++){
                        IAtom atom = atoms.get(i);
                        Vector position = atom.getPosition();
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
        int nAtoms = atoms.size();
        for (int i=0; i<nAtoms; i++){
            IAtom atom = atoms.get(i);
            double charge = atomAgentManager.getAgent(atom).charge;
            uSelf += charge*charge;
        }

        uSelf *= -alpha/sqrtPI;
        return uSelf;
    }

    public double uBondCorr(){
        double uCorr = 0.0;
        for (int i = 0; i< moleculeList.size(); i++){
            IMolecule molecule = moleculeList.get(i);
            int numSites = molecule.getChildList().size();
            for (int siteA=0; siteA<numSites; siteA++){
                IAtom atomA = molecule.getChildList().get(siteA);
                Vector positionA = atomA.getPosition();
                double chargeA = atomAgentManager.getAgent(atomA).charge;
                if (chargeA==0) continue;
                for (int siteB=siteA+1; siteB<numSites; siteB++){
                    IAtom atomB = molecule.getChildList().get(siteB);
                    Vector positionB = atomB.getPosition();
                    double chargeB = atomAgentManager.getAgent(atomB).charge;
                    if (chargeB==0) continue;
                    rAB.Ev1Mv2(positionA, positionB);
                    box.getBoundary().nearestImage(rAB);
                    double rABMagnitudeSquared = rAB.squared();
                    double rABMagnitude = Math.sqrt(rABMagnitudeSquared);
                    uCorr -= chargeA*chargeB*Erf.erf(alpha*rABMagnitude)/rABMagnitude;
                }
            }
        }		
        return uCorr;
    }

    public void setAlpha(double alpha){
        this.alpha = alpha;
        alpha2 = alpha*alpha;
        alpha3 = alpha2 * alpha;
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public int nBody() {
        return 0;
    }

    public void setBox(Box box) {
    }

    //******************************** inner class ********************************************//
    public static class MyCharge{
        public MyCharge(double charge){
            this.charge = charge;
        }
        public final double charge;
    }

    public double energy(IAtomList atoms) {
        double real = doRealSum ? uReal() : 0;
        double fourier = uFourier();
        double self = uSelf();
        double bondCorr = uBondCorr();

        double totalEnergy = real + fourier + self + bondCorr;
        if (false) {
            int numMolecules = moleculeList.size();
            System.out.println("total:               "+ totalEnergy/numMolecules);
            System.out.println("real   : "+ real/numMolecules);
            System.out.println("fourier: "+ fourier/numMolecules);
            System.out.println("self: "+ self/numMolecules);
            System.out.println("bond correction: "+ (-bondCorr/numMolecules));
        }

        return totalEnergy;
    }

    public double virial(IAtomList atoms) {
        return 0;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    //////////////////////////////////////////// begin calculating gradient //////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public Vector[] gradient(IAtomList atoms) {
        int nAtoms = box.getLeafList().size();
        double coeff = 4.0*Math.PI/volume;
        double kCutSquared = kCut * kCut; // criteria for spherical cutoff in fourier space

        if(gradient.length < nAtoms){
            gradient = new Vector[nAtoms];
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

        if (doRealSum) {
            //Real gradient  //Cross Interaction
            for (int i=0; i < nAtoms; i++){
                IAtom atomA = box.getLeafList().get(i);
                double chargeA = atomAgentManager.getAgent(atomA).charge;
                if (chargeA==0) continue;
                int aIndex = atomA.getParentGroup().getIndex(); // molecule a
                Vector positionA = atomA.getPosition();
                for (int j=i+1; j < nAtoms; j++){
                    IAtom atomB = box.getLeafList().get(j);
                    int bIndex = atomB.getParentGroup().getIndex(); // molecule b
    
                    if(nRealShells[0] == 0 && nRealShells[1] == 0 && nRealShells[2] == 0 && aIndex == bIndex) continue;//Skip same molecules!
    
                    double chargeB = atomAgentManager.getAgent(atomB).charge;
                    if (chargeB==0) continue;
                    Vector positionB = atomB.getPosition();
                    rAB.Ev1Mv2(positionA, positionB); //rAB == rA - rB
                    box.getBoundary().nearestImage(rAB);
                    for (int nx = -nRealShells[0]; nx <= nRealShells[0]; nx++) {
                        Lxyz.setX(0, nx*boxSize[0]); 
                        for (int ny = -nRealShells[1]; ny <= nRealShells[1]; ny++) {
                            Lxyz.setX(1, ny*boxSize[1]);
                            for (int nz = -nRealShells[2]; nz <= nRealShells[2]; nz++) {
                                if (aIndex==bIndex && nx*nx+ny*ny+nz*nz == 0) continue;
                                Lxyz.setX(2, nz*boxSize[2]);
                                drTmp.Ev1Pv2(rAB, Lxyz);
    
                                double rAB2 = drTmp.squared();
                                if (rAB2 > rCutSquared) continue; 
                                double rABMagnitude = Math.sqrt(rAB2);
                                double rAB3 = rABMagnitude*rAB2;
                                double B = Erf.erfc(alpha*rABMagnitude) + 2.0*alpha*rABMagnitude/sqrtPI * Math.exp(-alpha2*rAB2) ;
                                double realCoeff = - chargeA*chargeB * B / rAB3; // gradU = -F
                                gradient[i].PEa1Tv1(realCoeff, drTmp);
                                gradient[j].PEa1Tv1(-realCoeff, drTmp);
                            }
                        }
                    }
                }
            }
        }

        //Fourier gradient Part		

        for (int xAxis = -nKs[0]; xAxis < nKs[0]+1; xAxis++){
            kVector.setX(0, (xAxis * basis[0]));// assign value to the x-axis
            for (int yAxis = -nKs[1]; yAxis < nKs[1]+1; yAxis++ ){
                kVector.setX(1, (yAxis * basis[1]));// assign value to the y-axis
                for (int zAxis = -nKs[2]; zAxis < nKs[2]+1; zAxis++ ){
                    if( (xAxis * xAxis + yAxis * yAxis + zAxis * zAxis) == 0) continue;// first check: k is a non-zero vector
                    kVector.setX(2, (zAxis * basis[2]));// assign value to the z-axis, now the vector is specified
                    double kSquared = kVector.squared();
                    if (kSquared > kCutSquared) continue;// k-vector should be within the sphere with kCutoff as the radius
                    double sCoskr = 0.0, sSinkr = 0.0;

                    for (int j=0; j< nAtoms ; j++){ // Loop over atoms (4*nBasis)
                        IAtom atom = box.getLeafList().get(j);
                        Vector position = atom.getPosition();
                        sinkrj[j] = Math.sin(position.dot(kVector));
                        coskrj[j] = Math.cos(position.dot(kVector));
                        double chargej = atomAgentManager.getAgent(atom).charge;

                        sCoskr += chargej*coskrj[j];
                        sSinkr += chargej*sinkrj[j]; 
                    }//End loop over j

                    double coeffk = coeff / kSquared * Math.exp(-kSquared/4.0/alpha2);
                    for(int i=0; i<nAtoms; i++){
                        IAtom atom = box.getLeafList().get(i);
                        double chargei = atomAgentManager.getAgent(atom).charge;
                        double coeffki = coeffk * chargei * (sinkrj[i] * sCoskr - coskrj[i] * sSinkr); 
                        gradient[i].PEa1Tv1(-coeffki , kVector);  // gradU = -F
                    }
                }//end of storing Sin and Cos
            }
        }//End loop over ks

        //Intra-Molecular  gradient:
        for (int i = 0; i< moleculeList.size(); i++){
            IMolecule molecule = moleculeList.get(i);
            int numSites = molecule.getChildList().size();
            for (int siteA=0; siteA<numSites; siteA++){
                IAtom atomA = molecule.getChildList().get(siteA); // index = 0, 1, 2, 3|||leafIndex=0...184
                double chargeA = atomAgentManager.getAgent(atomA).charge;
                if (chargeA==0) continue;
                Vector positionA = atomA.getPosition();
                for (int siteB=siteA+1; siteB<numSites; siteB++){
                    IAtom atomB = molecule.getChildList().get(siteB);
                    double chargeB = atomAgentManager.getAgent(atomB).charge;
                    if (chargeB==0) continue;
                    Vector positionB = atomB.getPosition();

                    rAB.Ev1Mv2(positionA, positionB);
                    box.getBoundary().nearestImage(rAB);
                    double rAB2 = rAB.squared();
                    double rABMagnitude = Math.sqrt(rAB2);
                    double B = 2*alpha/sqrtPI * Math.exp(-alpha2*rAB2)-Erf.erf(alpha*rABMagnitude)/rABMagnitude; 
                    double coeffAB = - chargeA*chargeB * B / rAB2; // gradU = -F
                    gradient[atomA.getLeafIndex()].PEa1Tv1(coeffAB, rAB);
                    gradient[atomB.getLeafIndex()].PEa1Tv1(-coeffAB, rAB);
                }
            }
        }
        return gradient;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    //////////////////////////          begin calculating secondDerivatives               /////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public Tensor secondDerivative(IAtom atom0, IAtom atom1){
        double kCutSquared = kCut * kCut; // criteria for spherical cutoff in fourier space
        Vector pos0 = atom0.getPosition();
        Vector pos1 = atom1.getPosition();
        rAB.Ev1Mv2(pos0,pos1);
        box.getBoundary().nearestImage(rAB);

        double q0 = atomAgentManager.getAgent(atom0).charge;
        double q1 = atomAgentManager.getAgent(atom1).charge;
        secondDerivative.E(0);
        if (q0*q1 == 0) return secondDerivative;
        double coeff = 4.0*Math.PI/volume;

        for (int xAxis = -nKs[0]; xAxis < nKs[0]+1; xAxis++){
            kVector.setX(0, (xAxis * basis[0]));// assign value to the x-axis
            for (int yAxis = -nKs[1]; yAxis < nKs[1]+1; yAxis++ ){
                kVector.setX(1, (yAxis * basis[1]));// assign value to the y-axis
                for (int zAxis = -nKs[2]; zAxis < nKs[2]+1; zAxis++ ){
                    if( (xAxis * xAxis + yAxis * yAxis + zAxis * zAxis) == 0) continue;// first check: k is a non-zero vector
                    kVector.setX(2, (zAxis * basis[2]));// assign value to the z-axis, now the vector is specified
                    double kSquared = kVector.squared();
                    if (kSquared > kCutSquared) continue;// k-vector should be within the sphere with kCutoff as the radius
                    tempTensorkk.Ev1v2(kVector, kVector);
                    tempTensorkk.TE(Math.exp(-kSquared/4.0/alpha2) * Math.cos(kVector.dot(rAB))/kSquared);
                    secondDerivative.PE(tempTensorkk);
                }
            }
        }
        secondDerivative.TE(coeff);
        boolean intraMolecular = atom0.getParentGroup() == atom1.getParentGroup();

        for(int nx = -nRealShells[0]; nx <= nRealShells[0]; nx++) {
            Lxyz.setX(0, nx*boxSize[0]); 
            for(int ny = -nRealShells[1]; ny <= nRealShells[1]; ny++) {
                Lxyz.setX(1, ny*boxSize[1]);
                for(int nz = -nRealShells[2]; nz <= nRealShells[2]; nz++) {
                    Lxyz.setX(2, nz*boxSize[2]);
                    drTmp.Ev1Pv2(rAB, Lxyz);

                    double rAB2 = drTmp.squared();
                    if (rAB2 > rCutSquared && (!intraMolecular || (nx*nx+ny*ny+nz*nz > 0))) continue;
                    double rABMagnitude = Math.sqrt(rAB2);
                    double rAB3 = rABMagnitude*rAB2;
                    double rAB4 = rAB2*rAB2;
                    double rAB5 = rABMagnitude*rAB4;
                    double erfc = Erf.erfc(alpha*rABMagnitude);
                    double exp_Alpha2r2 = Math.exp(-alpha2*rAB2);
                    double B0 = (6*alpha/rAB4 + 4*alpha3/rAB2)/sqrtPI * exp_Alpha2r2;

                    tempTensorkk.Ev1v2(drTmp, drTmp);

                    if (intraMolecular && nx*nx+ny*ny+nz*nz==0) {
                        double A = (1.0-erfc)/rAB3 - 2*alpha/sqrtPI*exp_Alpha2r2/rAB2;
                        double B = B0 - 3 * (1-erfc)/rAB5;
                        tempTensorkk.TE(-B);
                        tempTensorkk.PEa1Tt1(-A,identity);//ch
                    }
                    else {
                        double A = erfc/rAB3 + 2*alpha/sqrtPI*exp_Alpha2r2/rAB2;
                        double B = B0 + 3*erfc/rAB5;
                        tempTensorkk.TE(-B);//ch
                        tempTensorkk.PEa1Tt1(A,identity);
                    }

                    secondDerivative.PE(tempTensorkk);
                }//nz
            }//ny
        }//nx

        secondDerivative.TE(q0*q1);
        return secondDerivative;
    }

    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

    public P2EwaldReal makeP2EwaldReal() {
        doRealSum = false;
        return new P2EwaldReal();
    }
    public static void computeScaledSysParams(double sigma_Er_sim,double alpha0,double volume, double Q2){
		double alpha = alpha0/Math.sqrt(7); //nX=2  
        double precision_splus = 0;
        double precision_sminus = Math.pow(-Math.log(sigma_Er_sim),2.0/3.0);
        double f = 1;
        double precision_s = 0;
        while (true){
            precision_s = (precision_sminus + precision_splus)/2;
            if(precision_s == precision_sminus || precision_s == precision_splus) break;
            f = Math.exp(-precision_s*precision_s)*Q2*Math.sqrt(precision_s/2/alpha/volume) - precision_s*precision_s*sigma_Er_sim;
            if (f==0) break;
            if(f > 0){
                precision_splus = precision_s;
            } 
            else {
                precision_sminus = precision_s;
            }
        }
        double scaledAlpha =  alpha;
        System.out.println("Scaled s = "+ precision_sminus); 
        System.out.println("Scaled alpha = "+ scaledAlpha); 
        System.out.println("Scaled rCutES = "+ precision_sminus/scaledAlpha);
        System.out.println("Scaled kCut = "+ 2.0*precision_sminus*scaledAlpha); 
        sigma_Er_sim = Q2*Math.sqrt(precision_sminus/2/scaledAlpha/volume)*Math.exp(-precision_sminus*precision_sminus)/precision_sminus/precision_sminus;
        double sigma_Er = Joule.UNIT.fromSim(sigma_Er_sim)*1.0E-3*Constants.AVOGADRO;
        System.out.println("sigma(Er) = Q/N*(s/2/alpha/L^3)^1/2 * exp(-s^2)/s/s  (kJ/mol) =  " + sigma_Er);
        System.out.println("sigma(Er_sim) = " + sigma_Er_sim);

    }
    public class P2EwaldReal implements PotentialSoft {

        protected final Vector[] gradient2;

        public P2EwaldReal() {
            gradient2 = new Vector[2];
            gradient2[0] = space.makeVector();
            gradient2[1] = space.makeVector();
        }

        public double energy(IAtomList atoms) {
            IAtom atomA = atoms.get(0);
            double chargeA = atomAgentManager.getAgent(atomA).charge;

            IAtom atomB = atoms.get(1);
            double chargeB = atomAgentManager.getAgent(atomB).charge;
            
            Vector positionA = atomA.getPosition();
            Vector positionB = atomB.getPosition();

            rAB.Ev1Mv2(positionA, positionB);// get vector rAB
            box.getBoundary().nearestImage(rAB);// minimum image

            double r2 = rAB.squared();
            if(r2 > rCutSquared) return 0;
            double r = Math.sqrt(r2);
            return chargeA * chargeB * Erf.erfc(alpha * r) / r;//Don't worry about 1/2 factor!
        }

        public double getRange() {
            return rCutRealES;
        }

        public void setBox(Box box) {}

        public int nBody() {
            return 2;
        }

        public double virial(IAtomList atoms) {
            return 0;
        }

        public Vector[] gradient(IAtomList atoms) {
            //Real gradient  //Cross Interaction
            IAtom atomA = atoms.get(0);
            double chargeA = atomAgentManager.getAgent(atomA).charge;
            Vector positionA = atomA.getPosition();
            IAtom atomB = atoms.get(1);

            double chargeB = atomAgentManager.getAgent(atomB).charge;

            Vector positionB = atomB.getPosition();
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
            double B = Erf.erfc(alpha*rABMagnitude) + 2.0*alpha*rABMagnitude/sqrtPI * Math.exp(-alpha2*rAB2) ;
            double realCoeff = - chargeA*chargeB * B / rAB3; // gradU = -F
            gradient2[0].Ea1Tv1(realCoeff, rAB);
            gradient2[1].Ea1Tv1(-realCoeff, rAB);
            return gradient2;
        }

        public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
            return gradient(atoms);
        }
    }
}
