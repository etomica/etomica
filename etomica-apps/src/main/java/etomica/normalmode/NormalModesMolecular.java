/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import java.io.IOException;
import java.io.PrintWriter;

import etomica.box.Box;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;
import etomica.atom.AtomPositionCOM;
import etomica.lattice.crystal.Primitive;
import etomica.models.water.SpeciesWater4P;
import etomica.normalmode.LatticeSumMolecularCrystal.AtomicTensorAtomicPair;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space3d.Tensor3D;
import etomica.spaceNd.TensorND;

/**
 * Uses analysis of 2nd derivatives to compute the normal modes for a Bravais lattice with a basis, 
 * occupied by atoms that interact with a simple, spherically-symmetric soft potential.
 */

public class NormalModesMolecular implements NormalModes {

    public NormalModesMolecular(SpeciesWater4P species, boolean waveVectorMethod , PotentialMaster potentialMaster, Box box, int[] nUniyCellsInSupBox, Primitive primitive, int basisDim, AtomicTensorAtomicPair atomicTensorAtomicPair, Space space) {
    	this.waveVectorMethod = waveVectorMethod;
        this.space = space;
        this.potentialMaster = potentialMaster;
        this.box = box;
        this.basisDim = basisDim;
        this.atomicTensorAtomicPair = atomicTensorAtomicPair;
        this.primitive = primitive;
        this.nC = nUniyCellsInSupBox;
        this.species = species;
    }
    
    public void calculateModes() {
    	int spaceDim = 6;
        int eDim = basisDim * spaceDim;
        WaveVectorFactorySimple kFactory = new WaveVectorFactorySimple(primitive, space);
        kFactory.makeWaveVectors(box);
        double[] kCoefficients = kFactory.getCoefficients(); //kCoefficients=0.5 non-deg.; = 1 degenerate twice!
        
//    	Vector[] kv = kFactory.getWaveVectors();
//    	for (int i=0;i<kv.length;i++){
//        	System.out.println(kv[i]);    		
//    	}
//    	System.exit(0);
    	Tensor[][] Inertia = new Tensor[basisDim][basisDim];
        for(int i=0; i<basisDim; i++) {
            for(int j=0; j<basisDim; j++) {
            	Inertia[i][j] = new TensorND(6);
            }
        }
    	IMoleculeList molList = box.getMoleculeList();
        Tensor tempTensor = space.makeTensor();
        Tensor inertiaTensor = space.makeTensor();
        double massH2O = species.getOxygenType().getMass() + 2.0 * species.getHydrogenType().getMass();
    	AtomPositionCOM comi = new AtomPositionCOM(space);
    	Vector drk = space.makeVector();
        Tensor identity = new Tensor3D(new double[][] {{1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0}});

        for(int i=0; i<basisDim; i++) {
        	IMolecule moleculei = molList.getMolecule(i);
        	for(int j=0;j<space.D();j++){ // 4 NOT 3 but that is fine as the mass of M = 0
            	Inertia[i][i].setComponent(j, j,massH2O);        		
        	}
        	Vector comPos = comi.position(moleculei);
        	inertiaTensor.E(0);
        	for(int j=0; j<moleculei.getChildList().getAtomCount(); j++){
        		drk.Ev1Mv2(moleculei.getChildList().getAtom(j).getPosition(),comPos);
        		//System.out.println(drk.dot(drk));
        		tempTensor.Ev1v2(drk, drk);
        		tempTensor.TE(-1);
        		tempTensor.PEa1Tt1(drk.dot(drk), identity);
        		tempTensor.TE(moleculei.getChildList().getAtom(j).getType().getMass());
        		inertiaTensor.PE(tempTensor);
        	}
    		for(int m=0;m<space.D();m++){
        		for(int n=0;n<space.D();n++){
            		Inertia[i][i].setComponent(m+3,n+3,inertiaTensor.component(m, n));
        		}
    		}
        }
//        double[][] arrayM = new double[eDim][eDim];
//    	outputStream.print("( "+arrayR[n][m]+ " , "+ arrayI[n][m]+" ) ");

        String fileM = "M";
        try{
        	PrintWriter outputStreamM = new PrintWriter(fileM);            	
	        for(int i=0; i<basisDim; i++) {
                for(int alpha=0; alpha<spaceDim; alpha++) {
                	
                	for(int j=0; j<basisDim; j++) {
	                    for(int beta=0; beta<spaceDim; beta++) {
	                        if(waveVectorMethod == true){
	                    		outputStreamM.print("( "+Inertia[i][j].component(alpha, beta)+ " , "+ 0.0 +" ) ");
	                    	}else{
	                    		outputStreamM.print(Inertia[i][j].component(alpha, beta)+ " ");
	                    	}
	                    }
	                }
                	outputStreamM.println("");
	            }
	        }
        	outputStreamM.close();
	    }catch (IOException e) {
	    	throw new RuntimeException(e);
	    }

            


    	
    	
    	
    	
    	
        summer = new LatticeSumMolecularCrystal(potentialMaster,box,space,basisDim,primitive);
        
        Tensor[][][][] sum = summer.calculateSum(atomicTensorAtomicPair);
        double[][] arrayR = new double[eDim][eDim];
        double[][] arrayI = new double[eDim][eDim];
       
        for(int k=0; k<kCoefficients.length; k++) {
            
            for(int j=0; j<basisDim; j++) {
                for(int jp=0; jp<basisDim; jp++) {
                    Tensor tensorj_jpR = sum[0][j][jp][k];
                    Tensor tensorj_jpI = sum[1][j][jp][k];
                    
                    for(int alpha=0; alpha<spaceDim; alpha++) {
                        for(int beta=0; beta<spaceDim; beta++) {
                            double vR = tensorj_jpR.component(alpha, beta);
                            double vI = tensorj_jpI.component(alpha, beta);
                            arrayR[spaceDim*j+alpha][spaceDim*jp+beta] = vR; // spaceDim = 6
                            arrayI[spaceDim*j+alpha][spaceDim*jp+beta] = vI; // spaceDim = 6
                        }
                    }
                }
            }//j
            
            String filek;
            if(waveVectorMethod == true){
            	filek = "D"+nC[0]+nC[1]+nC[2]+k+"_"+kCoefficients[k]+"deg";
            }else{
            	filek = "Dbig"+nC[0]+nC[1]+nC[2];
            }
            try{
            	PrintWriter outputStream = new PrintWriter(filek);            	
            	
            for(int n=0;n<eDim;n++){
                for(int m=0;m<eDim;m++){
//                	outputStream.print(String.format("%20.15e%+20.15ei ", arrayR[n][m], arrayI[n][m]));
                    if(waveVectorMethod == true){
                    	outputStream.print("( "+arrayR[n][m]+ " , "+ arrayI[n][m]+" ) ");
                	}else{
                    	outputStream.print(arrayR[n][m]+ " ");
                	}
                }
            	outputStream.println("");
            }
        	outputStream.close();

            }catch (IOException e) {
            	throw new RuntimeException(e);
            }
        }
    }
        

    public void setHarmonicFudge(double newHarmonicFudge) {
        // we ignore fudge
    }
    
    public void setTemperature(double newTemperature) {
        // we ignore temperature
    }
    
	public String getFileName() {
		return fileName;
	}

	public void setFileName(String filename) {
		this.fileName = filename;
	}

	@Override
	public WaveVectorFactory getWaveVectorFactory() {
		return summer.getWaveVectorFactory();
	}

	@Override
	public double[][] getOmegaSquared() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double[][][] getEigenvectors() {
		// TODO Auto-generated method stub
		return null;
	}

	
	
	
	
    protected final Space space;
    private boolean needToCalculateModes;
    private String fileName;
    protected Box box;
    protected int basisDim;
    protected final AtomicTensorAtomicPair atomicTensorAtomicPair;
    protected final Primitive primitive;
	protected LatticeSumMolecularCrystal summer;
	protected PotentialMaster potentialMaster;
	protected final boolean waveVectorMethod;
	protected SpeciesWater4P species;
	protected int[] nC;
    
}



