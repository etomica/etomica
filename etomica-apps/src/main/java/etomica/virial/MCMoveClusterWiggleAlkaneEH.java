/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.api.IMoleculeList;
import etomica.potential.PotentialMaster;
import etomica.util.random.IRandom;
import etomica.simulation.Simulation;
import etomica.api.ISpecies;
import etomica.space.Vector;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.space.Space;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Vector3D;
import etomica.util.Debug;

/**
 * Wiggle move for alkane, TraPPE-EH
 * hydrogens attached to the 3 carbons are updated, keeping bond lengths and angles fixed
 *   *********** outline ************
 *   // update j if j == 0 or j == n-1
 *   // update 3 H on j 
 *   // if j==0, update 2H on (j+1)
 *   // if j==n-1, update 2H on (j-1)
 *   
 *   && 0<j<(n-1) update j
 *   && update 2 H on j
 *   && if j==1, update 3 H on (j-1)
 *   && if j==(n-2), update 3 H on (j+1) 
 *   && if ( j!= (n-2), update 2H on (j+1) // j==1 or else
 *   && if ( j!= 1, update 2H on (j-1) // j==(n-1) or else
 *   ********************************
 * 
 * @author shu
 *
 */

public class MCMoveClusterWiggleAlkaneEH extends MCMoveMolecule {
    public MCMoveClusterWiggleAlkaneEH(Simulation sim, PotentialMaster potentialMaster, int nAtoms, Space _space) {
        this(potentialMaster,sim.getRandom(), 1.0, nAtoms, _space);

    }
  
    /**
     * Constructor for MCMoveAtomMulti.
     * @param parentIntegrator
     * @param nAtoms number of atoms to move in a trial.  Number of atoms in
     * box should be at least one greater than this value (greater
     * because first atom is never moved)
     */

    public MCMoveClusterWiggleAlkaneEH(PotentialMaster potentialMaster, IRandom random, double stepSize, int nAtoms, Space _space) {
        super(potentialMaster,random,_space, stepSize,Double.POSITIVE_INFINITY);
        this.space = _space;
        setStepSizeMax(Math.PI);
        energyMeter = new MeterPotentialEnergy(potential);
        rotateTensor = (RotationTensor3D)(space.makeRotationTensor());
        work1 = _space.makeVector();
        work2 = _space.makeVector();
        work3 = _space.makeVector();
    }

    public void setBox(Box p) {
    	super.setBox(p);
    	int nMolecules = box.getMoleculeList().getMoleculeCount();
    	selectedAtoms = new IAtom[nMolecules][8];// [0] is carbon, other 5 or 6 or 7 are H
    	positionSelectedAtoms = new Vector[8];   // [0] is carbon, other 5 or 6 or 7 are H
    	translationVectors = new Vector3D[nMolecules][8];// [0] is carbon, other 5 or 6 or 7 are H
    	for (int i=0; i<nMolecules; i++) {
    		for (int m = 0; m < 8 ; m++){
    			translationVectors[i][m] = space.makeVector();// declare all translationVectors here
    		}
    	}
    	energyMeter.setBox(p);
    }

    public void setSpecies(ISpecies newSpecies) {
    	species = newSpecies;
    }
    
    //note that total energy is calculated
    public boolean doTrial() {
    	uOld = energyMeter.getDataAsScalar();
    	wOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
    	IMoleculeList moleculeList = box.getMoleculeList();
        for(int i=0; i<moleculeList.getMoleculeCount(); i++) {//loop over all the molecules in the box
        	if (species != null && moleculeList.getMolecule(i).getType() != species) {
        		continue;
            }

        	IAtomList childList = moleculeList.getMolecule(i).getChildList();
        	int numChildren = childList.getAtomCount(); // number of total atoms in the i-th molecule
            numCarbons = (numChildren - 2)/3; // number of carbons in the i-th molecule
            
            int j = random.nextInt(numCarbons);
            while (    ( numCarbons == 3)  && (j == 1) ) {
            	j = random.nextInt(numCarbons);
            }
//            System.out.println("========================= Wiggle move =====================");
//            System.out.println("Due to the flex correction, the total nbr  of the molecule : " + moleculeList.getMoleculeCount());            
 //           System.out.println("the label of the molecule: " + i);            
//            System.out.println("random number j is : " + j);   
            
            ///// ######################################################################### ///////////////////
            //////////////////////////////// selectedAtoms & positionSelectedAtoms & translationVectors ///////
            ///// ######################################################################### ///////////////////
            selectedAtoms[i][0] = childList.getAtom(j);// j carbon
            selectedAtoms[i][1] = childList.getAtom(j+numCarbons);// H1 on j
            selectedAtoms[i][2] = childList.getAtom(j+numCarbons*2);// H2 on j
            
            for (int t = 3; t<8; t++){
            	selectedAtoms[i][t]=null;
            }
            
            if (  (j==0) || ( j==(numCarbons-1)  )|| ( j==1 )||  (  j==(numCarbons-2)  )    ){
            	numAtomsMove = 6;// 1 carbon and 5 H in total
            	if ( j == 0 ){
            		selectedAtoms[i][3] = childList.getAtom(numCarbons*3);//H3 on j=0
            		selectedAtoms[i][4] = childList.getAtom(j+1+numCarbons);//H1 on j+1=1
            		selectedAtoms[i][5] = childList.getAtom(j+1+numCarbons*2);//H2 on j+1=1
            	}
           
            	if (   j==(numCarbons-1)   ){
            		selectedAtoms[i][3] = childList.getAtom(numCarbons*3+1);//H3 on j=n-1
            		selectedAtoms[i][4] = childList.getAtom(j-1+numCarbons);//H1 on j-1=n-2
            		selectedAtoms[i][5] = childList.getAtom(j-1+numCarbons*2);//H2 on j-1=n-2
            	}
            	
            	if ( j == 1 ){
                	numAtomsMove = 8;// 1 carbon and 7 H in total
            		selectedAtoms[i][3] = childList.getAtom(j-1+numCarbons);//H1 on j-1=0   
            		selectedAtoms[i][4] = childList.getAtom(j-1+numCarbons*2);//H2 on j-1=0   
            		selectedAtoms[i][5] = childList.getAtom(j-1+numCarbons*3);//H3 on j-1=0
            		selectedAtoms[i][6] = childList.getAtom(j+1+numCarbons);//H1 on j+1=2
            		selectedAtoms[i][7] = childList.getAtom(j+1+numCarbons*2);//H2 on j+1=2
 
            	}
            
            	if ( j == (numCarbons-2) ){
                	numAtomsMove = 8;// 1 carbon and 7 H in total
                	selectedAtoms[i][3] = childList.getAtom(numCarbons*2-1);//H1 on j+1   
            		selectedAtoms[i][4] = childList.getAtom(numCarbons*3-1);//H2 on j+1   
            		selectedAtoms[i][5] = childList.getAtom(numCarbons*3+1);//H3 on j+1
            		selectedAtoms[i][6] = childList.getAtom(j-1+numCarbons);//H1 on j-1=n-3
            		selectedAtoms[i][7] = childList.getAtom(j-1+numCarbons*2);//H2 on j-1=n-3
   
            	}
                       	
            }
            else{ 
            	numAtomsMove = 7;// 1 carbon and 5 H in total
            	selectedAtoms[i][3] = childList.getAtom(j - 1 + numCarbons);//H1 on j-1   
            	selectedAtoms[i][4] = childList.getAtom(j - 1 + numCarbons*2);//H2 on j-1   
            	selectedAtoms[i][5] = childList.getAtom(j + 1 + numCarbons);//H1 on j+1   
            	selectedAtoms[i][6] = childList.getAtom(j + 1 + numCarbons*2);//H2 on j+1   
            }
            
            for ( int m = 0;m < numAtomsMove; m++){
                positionSelectedAtoms[m] = selectedAtoms[i][m].getPosition();// assign positionSelectedAtoms with each atom's position
                // let translationVectors hold the information of j carbon and hydrogens need to be moved( with minus sign)
                translationVectors[i][m].Ea1Tv1(-1, positionSelectedAtoms[m]);

            }
         
            double oldBondLength1 = 0, oldBondLength2 = 0;
            double angleHCH = 107.8 * Math.PI/180;
            double alpha = 0.5 * angleHCH;
            double CHBond = 1.10 * 0.5 ; // distance between carbon and middle point of carbon-hydrogen bond
            double theta = (random.nextDouble()-0.5)*stepSize;

            ///// ######################################################################### ///////////////////
            //////////////////////////////  if j is C(H3), work1, work2 and work3 /////////////////////////////
            ///// ######################################################################### ///////////////////
            if ( (j==0) || (j== (numCarbons-1))  ) {
            	
                // this puts atom j in a random orientation without changing the bond length
            	// work1 is the current vector from the bonded atom to atom j
            	work1.E(positionSelectedAtoms[0]);
            	
 //           	System.out.println("------Check---------------------");
 //           	System.out.println("before move, carbon j position is >>>>>>>> : " + positionSelectedAtoms[0]);
 //           	System.out.println("before move, work1 : same with carbon j >>>>>> : "+work1);
 
               
                if (j == 0) {// C(H3) on the left
                	work1.ME(childList.getAtom(j+1).getPosition());
                                   
                }
                
                else { // j==(n-1),C(H3) on the right
                    work1.ME(childList.getAtom(j-1).getPosition());
                }
                
  //              System.out.println("work1 = r(C0)-r(C1)"+work1);
                                
                double bondLength = Math.sqrt(work1.squared());// get the bond length
 //               System.out.println("C-C bond length:is |work1|:"+bondLength);
                if (Debug.ON && Debug.DEBUG_NOW) {
                    oldBondLength1 = bondLength;
                }
                //work2 is a vector perpendicular to work1.  it can be any perpendicular vector, but that just makes it harder!
                if (work1.getX(0)*work1.getX(0) < 0.5*bondLength*bondLength) {
                    // if work1 doesn't point in the X direction (mostly) then
                    // find a vector in the plane containing the X axis and work1
                    double a = -work1.getX(0)/bondLength;// sin(0.5*angle CCC)
                    work2.Ea1Tv1(a,work1);
                    work2.setX(0,work2.getX(0)+bondLength);
//                    System.out.println("work2 is a vector perpendicular to work1");

                }
                else {
                    // work1 does point in the X direction (mostly) so find a vector in the plane containing the Y axis and work1
                    double a = -work1.getX(1)/bondLength;
                    work2.Ea1Tv1(a,work1);
                    work2.setX(1,work2.getX(1)+bondLength);
 //                   System.out.println("work1 does point in the X direction (mostly) so find a vector in the plane containing the Y axis and work1");

                }

                //normalize
                work2.TE(bondLength/Math.sqrt(work2.squared()));
                //work3 is a vector normal to both work1 and work2
                work3.E(work1);
                work3.XE(work2);
                
                work3.TE(bondLength/Math.sqrt(work3.squared()));

                double phi = (random.nextDouble()-0.5)*Math.PI;
                work2.TE(Math.cos(phi));
                work2.PEa1Tv1(Math.sin(phi),work3);
                
                Vector axis = space.makeVector();
                Vector jNeighbor = space.makeVector();
                jNeighbor.Ev1Mv2(positionSelectedAtoms[0], work1); //position of C1, r(C[1]) = r(C[0]) - work1 or  r(C[n-2]) = r(C[n-1]) - work1
                axis.E(work2);// copy of work2
                axis.normalize();
                rotateTensor.setRotationAxis(axis, theta);

 //               IVectorMutable positioncopy = space.makeVector();
//                positioncopy.E(positionSelectedAtoms[0]);
                
                ///// ######################################################################### ///////////////////
                ///// #########  GET the new position of j          ########################### ///////////////////
                ///// ######################################################################### ///////////////////
  //          	System.out.println("c["+j+"] position before move is >>>>>>: "+positionSelectedAtoms[0]);

               
                positionSelectedAtoms[0].ME(jNeighbor);
                rotateTensor.transform(positionSelectedAtoms[0]);
                positionSelectedAtoms[0].PE(jNeighbor);
 //           	System.out.println("c["+j+"] position after move is : "+positionSelectedAtoms[0]);


                if (j == 0) {
                	///// ######################################################################### ///////////////////
                	///// #########  GET the new positions of 3 H on j(C[n-1])   ###################### ///////////////////
                    ///// ######################################################################### ///////////////////
                	IAtom[] jH = new IAtom[3];
         
                	for (int s=0;s<3; s++){
     //           		System.out.println("s is:"+s);
                    
                		jH[s] =  childList.getAtom(numCarbons * (s+1));
                        Vector r = jH[s].getPosition();
                	//	System.out.println("position of H before move is:"+r);
                		
 //               		IVectorMutable vector_chBefore = space.makeVector();
//                		vector_chBefore.Ev1Mv2(positioncopy, r);
//                		System.out.println("distance between C and H before move: " + Math.sqrt(vector_chBefore.squared()));
//                		if ( Math.abs((Math.sqrt(vector_chBefore.squared()) - 0.55 )) > 0.01){
 //                           throw new RuntimeException("|Math.sqrt(vector_chBefore.squared()) - 0.55 | > 0.01");
  //                    	}
                        r.ME(jNeighbor);
                        rotateTensor.transform(r);
                        r.PE(jNeighbor);
                	//	System.out.println("position of H after move is:"+r);
                		
  //              		IVectorMutable vector_ch = space.makeVector();
  //              		vector_ch.Ev1Mv2(positionSelectedAtoms[0], r);
  //              		System.out.println("distance between C and H after move: " + Math.sqrt(vector_ch.squared()));
   //             		if ( ( Math.abs( Math.sqrt(vector_ch.squared()) - 0.55 )) > 0.01){
   //                         throw new RuntimeException("|Math.sqrt(vector_ch.squared()) - 0.55 | > 0.01");
  //                    	}
                		

                	}
                	
                    ///// ######################################################################### ///////////////////
                    ///// #########  GET the new positions of 2 H on (j+1)= C1 ###################### ///////////////////
                    ///// ######################################################################### ///////////////////

                	Vector jPlus1 = childList.getAtom(j+1).getPosition();// (j+1) carbon position
                    Vector jPlus2 = childList.getAtom(j+2).getPosition();// (j+2) carbon position
                    Vector jPlus1_H1 = childList.getAtom(numCarbons+j+1).getPosition();//H1 on (j+1) carbon
                    Vector jPlus1_H2 = childList.getAtom(numCarbons*2+j+1).getPosition();//H2 on (j+1) carbon
                    Vector crossV = space.makeVector();
                    Vector jPlus1H1_q = space.makeVector();// r(H1)-r(q)
                    Vector jPlus1H2_q = space.makeVector();// r(H2)-r(q)

                    // 1.(r(H) - r(q) and [ r(H') - r(q)]
                    Vector jPlus1_j = space.makeVector();
                    Vector jPlus1_jPlus2 = space.makeVector();
                    jPlus1_j.Ev1Mv2(jPlus1, positionSelectedAtoms[0]);// r(j+1)-r(j)(new position)
                    jPlus1_jPlus2.Ev1Mv2(jPlus1, jPlus2);// r(j+1)-r(j+2)
                    crossV.E(jPlus1_j);
                    crossV.XE(jPlus1_jPlus2);
                    crossV.normalize();// normal vector of j-(j+1)-(j+2)
                    crossV.TE(Math.sin(alpha)*CHBond);
                    jPlus1H1_q.E(crossV);
                    jPlus1H2_q.E(crossV);
                    jPlus1H2_q.TE(-1.0);// has opposite sign
                    // 2.r(q) - r(j+1)
                    Vector q_jPlus1 = space.makeVector();
                    q_jPlus1.Ev1Pv2(jPlus1_jPlus2 , jPlus1_j);//[r(j+1)-r(j+2)] - [r(j+1)-r(j)]
                    q_jPlus1.normalize();
                    q_jPlus1.TE(CHBond*Math.cos(alpha));//r(q)-r(j+1)
                    // r(H) = r1 + r2 + r(j+1); r(H') = r1' + r2 + r(j+1)
                    jPlus1_H1.E(q_jPlus1);
                    jPlus1_H2.E(q_jPlus1);
                    jPlus1_H1.PE(jPlus1H1_q);
                    jPlus1_H2.PE(jPlus1H2_q);
                    jPlus1_H1.PE(jPlus1);
                    jPlus1_H2.PE(jPlus1);
                    
                    
                    /*check part
                  	IVectorMutable cNeighbor = childList.getAtom(j+1).getPosition();
                  	IVectorMutable vector_c0_c1 = space.makeVector(); 
                  	vector_c0_c1.Ev1Mv2(cNeighbor, positionSelectedAtoms[0]);
                  	System.out.println("the distance between j and its neighbor after move is: "+Math.sqrt(vector_c0_c1.squared()));
                	if ( ( Math.abs(Math.sqrt(vector_c0_c1.squared()) - 1.535 ) )> 0.01){
                        throw new RuntimeException("|Math.sqrt(vector_c0_c1.squared()) - 1.535 | > 0.01");
                  	}
                	
                	IVectorMutable v_c1_H1 = space.makeVector();
                	IVectorMutable v_c1_H2 = space.makeVector();
                	v_c1_H1.Ev1Mv2(cNeighbor, jPlus1_H1);
                	v_c1_H2.Ev1Mv2(cNeighbor, jPlus1_H2);
                  	System.out.println("the distance between j+1 and its H1 after move is: "+Math.sqrt(v_c1_H1.squared()));
                  	System.out.println("the distance between j+1 and its H2 after move is: "+Math.sqrt(v_c1_H2.squared()));
                  	if ( ( Math.abs(  Math.sqrt(v_c1_H1.squared()) - 0.55 ))  > 0.01){
                        throw new RuntimeException("|Math.sqrt(v_c1_H1.squared()) - 0.55 | > 0.01");
                  	}
                  	if ( (Math.sqrt(v_c1_H2.squared()) - 0.55 ) > 0.01){
                        throw new RuntimeException("|Math.sqrt(v_c1_H2.squared()) - 0.55 | > 0.01");
                  	}
                   */
                }           


                if (j==(numCarbons-1)){
                	///// ######################################################################### ///////////////////
                	///// #########  GET the new positions of 3 H on j = n-1  ###################### ///////////////////
                	///// ######################################################################### ///////////////////
                	IAtom[] hydrogen = new IAtom[3];
                	hydrogen[0]=childList.getAtom(numCarbons*2-1);
                	hydrogen[1]=childList.getAtom(numCarbons*3-1);
                	hydrogen[2]=childList.getAtom(numCarbons*3+1);

                	for (int s=0;s<3; s++){
                		Vector r = hydrogen[s].getPosition();
                		
                	//	r.ME(positionNeighbor);//position C[n-2]
                	//	rotateTensor.transform(r);
                	//	r.PE(positionNeighbor);
                		
                		r.ME(jNeighbor);//position C[n-2]
                		rotateTensor.transform(r);
                		r.PE(jNeighbor);
                	}
                	//////////////////// check the distance between c[n-2]-H[n-1],1,2,3 //////////////////////////////////////
                	//////////////////// check the distance between c[n-2]-H[n-1],1,2,3 //////////////////////////////////////

                	///// ######################################################################### ///////////////////
                	///// #########  GET the new positions of 2 H on (j-1)= (n-2) ###################### ///////////////////
                	///// ######################################################################### ///////////////////
                	Vector jMinus1 = childList.getAtom(j-1).getPosition();// (j-1) position
                	Vector jMinus2 = childList.getAtom(j-2).getPosition();// (j-2)  position
                	Vector jMinus1_H1 = childList.getAtom(j-1 + numCarbons).getPosition();// H1 on (j-1)
                	Vector jMinus1_H2 = childList.getAtom(j-1 + numCarbons * 2 ).getPosition();// H2 on (j-1)

                	Vector H1_q = space.makeVector();
                	Vector H2_q = space.makeVector();
                	Vector jMinus1_jMinus2 = space.makeVector();
                	Vector jMinus1_j = space.makeVector();
                	jMinus1_jMinus2.Ev1Mv2(jMinus1, jMinus2);// r(j-1)-r(j-2)
                	jMinus1_j.Ev1Mv2(jMinus1, positionSelectedAtoms[0]);//r(j-1)-r(j)
                	
                	Vector crossV = space.makeVector();
                	crossV.E(jMinus1_jMinus2);
                	crossV.XE(jMinus1_j);
                	crossV.normalize();// unit normal vector of (j-2)-(j-1)-j plane
                	crossV.TE(Math.sin(alpha)*CHBond);
                	H1_q.E(crossV);
                	H2_q.E(crossV);
                	H2_q.TE(-1.0);// has opposite sign
                	
                	
                	Vector q_jMinus1 = space.makeVector();
                	q_jMinus1.Ev1Pv2(jMinus1_j, jMinus1_jMinus2);//r(q)= [r(j-1)-r(j)] + [r(j-1)-r(j-2)]
                	q_jMinus1.PE(jMinus1_jMinus2);//r(q)= [r(j-1)-r(j)] + [r(j-1)-r(j-2)]
                	q_jMinus1.normalize();
                	q_jMinus1.TE(CHBond*Math.cos(alpha));

                	// r(H) = r1 + r2 + r(j+1); r(H') = r1' + r2' + r(j+1)
                	jMinus1_H1.E(q_jMinus1);
                	jMinus1_H2.E(q_jMinus1);
                	jMinus1_H1.PE(H1_q);
                	jMinus1_H2.PE(H2_q);
                	jMinus1_H1.PE(jMinus1);
                	jMinus1_H2.PE(jMinus1);
                	
                            	
                }        
                
                //  translationVectors[i][0].PE(moveMolecules[0])
                for ( int s= 0; s<numAtomsMove; s++){
                	translationVectors[i][s].PE(positionSelectedAtoms[s]);
                }
                work1.E(translationVectors[i][0]);// based on carbon
                work1.TE(1.0/numCarbons);//numCarbons// keep COM at the same place
                for (int k=0; k<childList.getAtomCount(); k++) {
   	             	childList.getAtom(k).getPosition().ME(work1);
                
                } 
              
            }// end loop of j==0 or j==(n-1)
            
            
///////////////////////// BEGIN j:1~(n-2)///////////////////// ///////////////// ///////////// 
///////////////////////// BEGIN j:1~(n-2)///////////////////// ///////////////// ///////////// 
///////////////////////// BEGIN j:1~(n-2)///////////////////// ///////////////// ///////////// 
            
            else {
  //          	System.out.println("crankshaft move, the position of the atom to be moved is : "+positionSelectedAtoms[0]);
            	// crankshaft move.  atom j is rotated around the j-1 - j+1 bond.
            	// j-1 - j and j - j+1 bond lengths are unaltered.
//                System.out.println("middle move "+j+" "+position);

            	Vector jMinus1 = childList.getAtom(j-1).getPosition();// (j-1) carbon position
            	Vector jPlus1  = childList.getAtom(j+1).getPosition();// (j+1) carbon position
            	
            	work1.Ev1Mv2(jMinus1, positionSelectedAtoms[0]);
                work2.Ev1Mv2(jPlus1, positionSelectedAtoms[0]);
                if (Debug.ON && Debug.DEBUG_NOW) {
                    oldBondLength1 = Math.sqrt(work1.squared());
                    oldBondLength2 = Math.sqrt(work2.squared());
                }
                double cosTheta = work1.dot(work2)/(Math.sqrt(work1.squared()*work2.squared()));
                if (cosTheta < -0.999) {
                    // current bond angle is almost 180degrees, making crankshaft
                    // difficult to do precisely, so skip it.  we'll explore this
                    // degree of freedom some other time when the bond angle is
                    // different
                	for ( int t = 0; t<positionSelectedAtoms.length; t++){
                        translationVectors[i][t].E(0);
                	}
                    continue;
                }
                work1.TE(-1);
                work2.Ev1Mv2(jPlus1, jMinus1);
                work2.TE(work1.dot(work2)/work2.squared());
                // work2 is now the projection of r01 onto r02
                // place atom1... well, here
                work1.ME(work2);
                work2.XE(work1);
                work2.TE(Math.sqrt(work1.squared()/work2.squared()));
                // work1 is perpendicular to r02 and goes from the line
                // connecting 0 and 2 to 1.  work2 is perpendicular to r02 and
                // work1 and has the same length as work1.
                                
                ///// ######################################################################### ///////////////////
                ///// #########  GET the new position of j=positionSelectedAtoms[0]  ########## ///////////////////
                ///// ######################################################################### ///////////////////
                Vector axis = space.makeVector();
                axis.Ev1Mv2(jPlus1, jMinus1); // r(j+1)-r(j-1)
                axis.normalize();
                rotateTensor.setRotationAxis(axis, theta);
                positionSelectedAtoms[0].ME(jMinus1);
                rotateTensor.transform(positionSelectedAtoms[0]);
                positionSelectedAtoms[0].PE(jMinus1); 
             
              	///// ######################################################################### ///////////////////
            	///// #########  GET the new positions of 2H on j       ####################### ///////////////////
            	///// ######################################################################### ///////////////////
                Vector j_H1 = childList.getAtom(j + numCarbons ).getPosition();  // H1 on j
                Vector j_H2 = childList.getAtom(j + numCarbons * 2).getPosition();// H2 on j
                
                Vector j_jMinus1 = space.makeVector();// r(j)-r(j-1)
                j_jMinus1.Ev1Mv2(positionSelectedAtoms[0], jMinus1);
                Vector j_jPlus1 = space.makeVector();// r(j)-(j+1)
                j_jPlus1.Ev1Mv2(positionSelectedAtoms[0], jPlus1);
                
            	Vector crossV = space.makeVector();// helper
            	Vector H1_q = space.makeVector();// helper
            	Vector H2_q = space.makeVector();//helper
            	
            	crossV.E(j_jMinus1);
                crossV.XE(j_jPlus1);
                crossV.normalize(); // unit normal vector of (j-1)--j--(j+1) plane
                crossV.TE(Math.sin(alpha)*CHBond);
                H1_q.E(crossV);
                H2_q.E(crossV);
                H2_q.TE(-1.0);// has opposite sign with H1_q

                Vector q_j = space.makeVector();// r(q) - r(j)
                q_j.Ev1Pv2(j_jMinus1, j_jPlus1);
                q_j.normalize();
                q_j.TE(CHBond*Math.cos(alpha));

                // r(H) = r1 + r2 + r(j); r(H') = r1' + r2' + r(j)
                j_H1.E(q_j);
                j_H2.E(q_j);
                j_H1.PE(H1_q);
                j_H1.PE(positionSelectedAtoms[0]);
                j_H2.PE(H2_q);
                j_H2.PE(positionSelectedAtoms[0]);
                 
                ///// ######################################################################### //////////////
                //////   GET the new positions of 3 H related with j=1 or j=(n-1)              ///////////////
            	///// ######################################################################### //////////////
            	IAtom[] hydrogen = new IAtom[3];
            	if (j==1){
            		for (int s=0;s<3; s++){
                		hydrogen[s] = childList.getAtom(numCarbons * (s+1));
                		Vector r = hydrogen[s].getPosition();
                		r.ME(jMinus1);
                		rotateTensor.transform(r);
                		r.PE(jMinus1);
                	}

            	}
                if (j ==(numCarbons-2)) {
                	hydrogen[0]=childList.getAtom(numCarbons*2-1);
                	hydrogen[1]=childList.getAtom(numCarbons*3-1);
                	hydrogen[2]=childList.getAtom(numCarbons*3+1);
                	for (int s=0;s<3; s++){
                		Vector r = hydrogen[s].getPosition();
                		r.ME(jMinus1);
                   		rotateTensor.transform(r);
                   		r.PE(jMinus1);
                	} 
                	
                }
                
                if (  j != (numCarbons -2 )){
                	///// ######################################################################### ///////////////////
                	///// #########  GET the new positions of 2H on (j+1)              ############ ///////////////////
                	///// ######################################################################### ///////////////////
                	Vector jPlus2 = childList.getAtom(j+2).getPosition();// (j+2) carbon position
                	Vector jPlus1_H1 = childList.getAtom( j + 1 + numCarbons).getPosition();     // H1 on j+1
                	Vector jPlus1_H2 = childList.getAtom( j + 1 + numCarbons * 2).getPosition(); // H2 on j+1
                    	
                   	Vector jPlus1_j = space.makeVector();//r(j+1)-r(j)
                    Vector jPlus1_jPlus2 = space.makeVector();// r(j+1)-r(j+2)
                	jPlus1_j.Ev1Mv2(jPlus1, positionSelectedAtoms[0]);
                	jPlus1_jPlus2.Ev1Mv2(jPlus1, jPlus2);
                	
                    crossV.E(jPlus1_j);
                    crossV.XE(jPlus1_jPlus2);
                    crossV.normalize();// unit normal vector of j-(j+1)-(j+2)
                    crossV.TE(Math.sin(alpha)*CHBond);
                    H1_q.E(crossV);
                    H2_q.E(crossV);
                    H2_q.TE(-1.0);
                    Vector q_jPlus1 = space.makeVector();//r(q)-r(j+1)
                    q_jPlus1.Ev1Pv2(jPlus1_jPlus2,jPlus1_j); 
                    q_jPlus1.normalize();
                    q_jPlus1.TE(CHBond*Math.cos(alpha));
                    // r(H) = r1 + r2 + r(j+1); r(H') = r1' + r2' + r(j+1)
                    jPlus1_H1.E(q_jPlus1);
                    jPlus1_H2.E(q_jPlus1);
                    jPlus1_H1.PE(H1_q);
                    jPlus1_H2.PE(H2_q);
                    jPlus1_H1.PE(jPlus1);
                    jPlus1_H2.PE(jPlus1);
                	
                	
                }
                if (  j != 1 ){                	
                	///// ######################################################################### ///////////////////
                	///// #########  GET the new positions of 2H on (j-1)              ############ ///////////////////
                	///// ######################################################################### ///////////////////
                	Vector jMinus2 = childList.getAtom(j-2).getPosition();// (j-2) carbon position
                	Vector jMinus1_H1 = childList.getAtom(j - 1 + numCarbons).getPosition();     // H1 on (j-1)
                	Vector jMinus1_H2 = childList.getAtom(j - 1 + numCarbons * 2 ).getPosition();// H2 on (j-1)
                
                	Vector jMinus1_jMinus2 = space.makeVector();// r(j-1)-r(j-2)
                   	Vector jMinus1_j = space.makeVector();//r(j-1)-r(j)
                	jMinus1_jMinus2.Ev1Mv2(jMinus1, jMinus2);
                	jMinus1_j.Ev1Mv2(jMinus1, positionSelectedAtoms[0]);//r(j-1)-r(j)
                
                    crossV.E(jMinus1_jMinus2);
                    crossV.XE(jMinus1_j);
                    crossV.normalize();// unit normal vector of (j-2)-(j-1)-j plane
                    crossV.TE(Math.sin(alpha)*CHBond);
                    H1_q.E(crossV);
                    H2_q.E(crossV);
                    H2_q.TE(-1.0);// has opposite sign with H1_q
                   
                    Vector q_jMinus1 = space.makeVector();
                    q_jMinus1.Ev1Pv2(jMinus1_j,jMinus1_jMinus2);// r(q) = [r(j-1)+r(j)]+[r(j-1)+r(j-2)]
                    q_jMinus1.normalize();
                    q_jMinus1.TE(CHBond*Math.cos(alpha));
                    // r(H) = r1 + r2 + r(j+1); r(H') = r1' + r2' + r(j+1)
                    jMinus1_H1.E(q_jMinus1);
                    jMinus1_H2.E(q_jMinus1);
                    jMinus1_H1.PE(H1_q);
                    jMinus1_H2.PE(H2_q);
                    jMinus1_H1.PE(jMinus1);
                    jMinus1_H2.PE(jMinus1);
                	
                	
                }
                
                
                ////// move the 1st carbon back to the origin///////////////////////////
                ////// keep the center of mass unchanged ///////////////////////////////
                for ( int s= 0; s<numAtomsMove; s++){
                    translationVectors[i][s].PE(positionSelectedAtoms[s]);//translationVectors is now position(new)-position(old)
                }
                work1.E(translationVectors[i][0]);//work1: translation of COM based on carbon
                work1.TE(1.0/numCarbons);
                for (int k=0; k<childList.getAtomCount(); k++) {
                    childList.getAtom(k).getPosition().ME(work1);
                }
                
            }

            if (Debug.ON && Debug.DEBUG_NOW) {
                if (j > 0) {
                    work1.Ev1Mv2(positionSelectedAtoms[0], childList.getAtom(j-1).getPosition());
                    double bondLength = Math.sqrt(work1.squared());
                    if (Math.abs(bondLength - oldBondLength1)/oldBondLength1 > 0.000001) {
                        throw new IllegalStateException("wiggle "+i+" "+j+" bond length should be close to "+oldBondLength1+" ("+bondLength+")");
                    }
                }
                if (j < numCarbons-1) {
                    work1.Ev1Mv2(positionSelectedAtoms[0], childList.getAtom(j+1).getPosition());
                    double bondLength = Math.sqrt(work1.squared());
                    double oldBondLength = oldBondLength2 == 0 ? oldBondLength1 : oldBondLength2;
                    if (Math.abs(bondLength - oldBondLength)/oldBondLength > 0.000001) {
                        throw new IllegalStateException("wiggle "+i+" "+j+" bond length should be close to "+oldBondLength+" ("+bondLength+")");
                    }
                }
            }
        }

        ((BoxCluster)box).trialNotify();
        wNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        uNew = energyMeter.getDataAsScalar();
        return true;
    }
    public void acceptNotify() {
        ((BoxCluster)box).acceptNotify();
    }

    public void rejectNotify() {
  //  	System.out.println("reject");
        IMoleculeList moleculeList = box.getMoleculeList();
        for(int i=0; i<moleculeList.getMoleculeCount(); i++) {
            if (species != null && moleculeList.getMolecule(i).getType() != species) continue;
            IAtomList childList = moleculeList.getMolecule(i).getChildList();
            work1.E(translationVectors[i][0]);//based on carbon 
            work1.TE(1.0/numCarbons);// Center of mass is based on all carbons
            for (int k=0; k<childList.getAtomCount(); k++) {
                childList.getAtom(k).getPosition().PE(work1);
            }
            
            for (int t = 0; t<8; t++){// put all the moved atoms (carbon and hydrogens) back
               if (selectedAtoms[i][t]==null){
            	   break;
               }
            	selectedAtoms[i][t].getPosition().ME(translationVectors[i][t]);
                
            }  
 //           System.out.println("position of carbon is : "+ selectedAtoms[i][0]);
 //           System.out.println("position of carbon is : "+ selectedAtoms[i][0].getPosition());
 //           System.out.println("position of carbon is : "+ positionSelectedAtoms[0]);

        }
        ((BoxCluster)box).rejectNotify();
    }

    public double getB() {
        return -(uNew - uOld);
    }

    public double getA() {
        return wNew/wOld;
    }
    private static final long serialVersionUID = 1L;
    protected final MeterPotentialEnergy energyMeter;
    protected IAtom[][] selectedAtoms;
    protected final Vector work1, work2, work3;
    protected Vector[][] translationVectors;
    protected double wOld, wNew;
    protected final Space space;
    protected ISpecies species;
    protected RotationTensor3D rotateTensor;
    protected int numAtomsMove ;
    protected int numCarbons;
    protected Vector[] positionSelectedAtoms;
}
