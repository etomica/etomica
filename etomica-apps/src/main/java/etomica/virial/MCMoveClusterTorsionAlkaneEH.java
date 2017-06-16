/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.molecule.*;
import etomica.potential.P4BondTorsion;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;

/**
 * MC move for alkane-TraPPE-EH torsion
 * pick up carbons in the box only, but move update all the atoms(carbons and hydrogens)
 * 
 * @author shu
 * Feb 2013
 */
public class MCMoveClusterTorsionAlkaneEH extends MCMoveMolecule {
    public MCMoveClusterTorsionAlkaneEH(Simulation sim, PotentialMaster potentialMaster, Space space, P4BondTorsion torsionPotential) {
    this(potentialMaster, space, sim.getRandom(), 1.0,torsionPotential, 20);
        setBondLength(1.0);
    }
    /**
     * Constructor for MCMoveAtomMulti.
     * @param parentIntegrator
     * @param nAtoms number of atoms to move in a trial.  Number of atoms in box should be at least one greater than this value (greater
     * because first atom is never moved)
     */
    public MCMoveClusterTorsionAlkaneEH(PotentialMaster potentialMaster, Space space,
                                        IRandom random, double stepSize, P4BondTorsion torsionPotential, int nBins) {
    	super(potentialMaster,random,space,stepSize,Double.POSITIVE_INFINITY);
        ((MCMoveStepTracker)getTracker()).setTunable(false);
        probabilityBins = new double[nBins+1];
        binSize = new double[nBins];
        probabilityReverseMap = new int[nBins+1];
        this.torsionPotential = torsionPotential;
        setStepSizeMax(Math.PI);
        energyMeter = new MeterPotentialEnergy(potential);
        work1 = space.makeVector();
        work2 = space.makeVector();
        work3 = space.makeVector();
        dr21 = space.makeVector();
        dr23 = space.makeVector();
        dr34 = space.makeVector();
        oldCenter = space.makeVector();
    }
    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
        ((MCMoveStepTracker)getTracker()).setTunable(false);
    }
    public void setTemperature(double temperature) {
        int nBins = probabilityBins.length - 1;
        int nSubBins = 100;
        double beta = 1.0/temperature;
        // numerically integrate P = exp(-beta U) from cosphi=1 to cosphi=-1 (0 to pi radians)
        double totP = 0.5 * Math.exp(-beta*torsionPotential.energyAtAngle(1));
        for (int i=1; i<nSubBins * nBins; i++) {
            double cosphi = Math.cos((Math.PI*i)/(nSubBins*nBins));
            double u = torsionPotential.energyAtAngle(cosphi);
            totP +=  Math.exp(-beta*u);
         
        }
        
        totP += 0.5 * Math.exp(-beta * torsionPotential.energyAtAngle(-1));
        // pPerBin is the amount of probability allocated to each bin.
        double pPerBin = totP / nBins;
        double thisBinP = 0;
        double previousP = 1;
        double newTot = 0;
        int iBin = 1;
        // last bin always ends at phi=PI
        int lastMapBin = 0;
        // probabilityBins is an array of the boundaries of the bins (in terms of cos phi)
        // bin i goes from probabilityBins[i] to probabilityBins[i+1]
        // the first bin is bin 0 and the last bin is bin nBins-1
        probabilityBins[0] = 0;
        probabilityBins[nBins] = Math.PI;
        // probabilityReverseMap maps an angle to a bin
        // probabilityReverseMap[phi / PI * nBins] will be less than or
        //   equal to the appropriate bin.
        // probabilityReverseMap[(phi / PI * nBins) + 1] will be greater than or
        //   equal to the appropriate bin.
        probabilityReverseMap[0] = 0;
        // sum up probability (=exp(-beta U)) until we exceed pPerBin.  Then
        // interpolate to find approximate when the sum is equal to pPerBin.
        // Call that the bin boundary and then begin summing for the next bin.
        for (int i=1; i<nSubBins * nBins + 1; i++) {
            double cosphi = Math.cos((Math.PI*i)/(nSubBins*nBins));
            double u = torsionPotential.energyAtAngle(cosphi);
            double thisP = Math.exp(-beta*u);
            double newP = 0.5 * (previousP + thisP);
            newTot += newP;
            
                  
   //         System.out.println(i+" "+cosphi+" "+u+" "+newTot+" "+thisBinP+" "+pPerBin);
            // if i is too much probability, then ix would be the value
            // of i (based on interpolation) where we hit the limit
            double ix = (i-1) + (pPerBin-thisBinP)/newP;
            while (ix < i) {
                // we have too much probability for this bin.  back up (interpolate)
                // to find where we hit the limit for this bin.  dump the rest into the
                // next bin
                probabilityBins[iBin] = (ix*Math.PI)/(nSubBins*nBins);//????????????????????/
                
                // ::::::::::::::::::::: check probability distribution ::::::::::::::::::::::::::::: //
               // System.out.print((Math.PI*i)/(nSubBins*nBins));// phi
              //  System.out.print("    ");
               // System.out.print(cosphi);
               // System.out.print("    ");
               // System.out.println(thisP);//probability exp(-beta*U)
               // System.out.print(iBin);//
                //System.out.print("    ");
             //   System.out.print(((double)iBin)/nBins);//x-axis
             //   System.out.print("    ");
            //    System.out.println(probabilityBins[iBin]);//angle,y-axis
                // ::::::::::::::::::::: end of check probability distribution ::::::::::::::::::::::::::::: //

                
                
                // we'll use this to correct acceptance (which is not correct due
                // to energy inhomogeniety within the bin)
                binSize[iBin-1] = (probabilityBins[iBin]-probabilityBins[iBin-1])/(Math.PI/nBins);
                // we found the lower bound of iBin.  Now mark 
                int nextMapBin = (int)((probabilityBins[iBin])/Math.PI * nBins);
                for (int j=lastMapBin+1; j<nextMapBin+1; j++) {
                    probabilityReverseMap[j] = iBin-1;
//                    System.out.println(j+" "+Math.acos(-1.0+(2.0*j)/nBins)+" "+(iBin-1));
                }
                lastMapBin = nextMapBin;
//                System.out.println(iBin+" "+probabilityBins[iBin]+" "+Math.acos(probabilityBins[iBin]) +" "+ torsionPotential.energyAtAngle(probabilityBins[iBin]));
//                double phi = Math.acos((probabilityBins[iBin-1]+probabilityBins[iBin])*0.5);
//                System.out.println(phi+" "+ torsionPotential.energyAtAngle(Math.cos(phi)) + " "+ (probabilityBins[iBin-1]-probabilityBins[iBin]));
                iBin++;
                if (iBin == nBins) {
                    // we just found the end of the next-to-the-last bin
                    // we know the last bin ends at PI (and already marked it as such)
//                    phi = Math.acos((probabilityBins[iBin-1]+probabilityBins[iBin])*0.5);
//                    System.out.println(phi+" "+ torsionPotential.energyAtAngle(Math.cos(phi))+" "+(probabilityBins[iBin-1]-probabilityBins[iBin]));
                    // we'll use this to correct acceptance (which is not correct due
                    // to energy inhomogeniety within the bin)
                    binSize[iBin-1] = (probabilityBins[iBin]-probabilityBins[iBin-1])/(Math.PI/nBins);
                    for (int j=lastMapBin+1; j<nBins+1; j++) {
                        probabilityReverseMap[j] = iBin-1;
//                        System.out.println(j+" "+Math.acos(-1.0+(2.0*j)/nBins)+" "+(iBin-1));
                    }
                    return;
                }
                thisBinP -= pPerBin;
                // if we still have too much probability (for the next bin)
                // ix will be the the new value where we hit the limit
                ix += pPerBin/newP;
            }
            thisBinP += newP;
            previousP = thisP;
        }
        throw new RuntimeException("oops");
    } 
    
    
    public void setSpecies(ISpecies newSpecies) {
        species = newSpecies;
        positionDefinition = new MoleculePositionGeometricCenterAlkaneEH(space,species);

    }
    
    //note that total energy is calculated
    public boolean doTrial() {
   // 	System.out.println("--------------- torsion move -------------------------------");
        if (selectedMolecules == null) selectMolecules();
    //    System.out.println("<<<<<<<<<<<<<<<<<<<<<<<<<<<<before move");

        uOld = energyMeter.getDataAsScalar();

        wOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
  //      System.out.println("wOld:" + wOld);
        bias = 1;
        for(int i=0; i<selectedMolecules.getMoleculeCount(); i++) {
            oldCenter.E(positionDefinition.position(selectedMolecules.getMolecule(i)));
            IAtomList childList = selectedMolecules.getMolecule(i).getChildList();
            int numChildren = childList.getAtomCount();
            int numCarbons = (numChildren-2)/3;
            int j = random.nextInt(numCarbons-3); // pick up an atom from all the carbons

  //          System.out.println("j is : "+j);
            // j=0 ==> first torsion bond (atoms 0,1,2,3)
            // rotate about atoms j+1 and j+2 (stay fixed).  atoms 0 to j and (j+3) to (N-1): move
            IAtom atom0 = childList.getAtom(j+0);
            IAtom atom1 = childList.getAtom(j+1);
            IAtom atom2 = childList.getAtom(j+2);
            IAtom atom3 = childList.getAtom(j+3);
            Vector atom1Position = atom1.getPosition();// use for rotation tensor, the "fixed point" on the rotating axis

            // axis: (j+2)-(j+1)
            RotationTensor3D rotationTensor = (RotationTensor3D)(space.makeRotationTensor());
            Vector axis = space.makeVector();
            axis.Ev1Mv2(atom2.getPosition(), atom1Position);//r(j+2)-r(j+1)
            axis.normalize();
            dr21.Ev1Mv2(atom0.getPosition(), atom1.getPosition());
            dr23.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
            dr34.Ev1Mv2(atom3.getPosition(), atom2.getPosition());

            double dr23Sq = dr23.squared();
            dr23.TE(1.0/Math.sqrt(dr23Sq));
            dr21.PEa1Tv1(-dr21.dot(dr23), dr23);
            dr34.PEa1Tv1(-dr34.dot(dr23), dr23);
            // current angle
            double oldcosphi = dr21.dot(dr34)/Math.sqrt(dr21.squared()*dr34.squared());
  //          System.out.println("j:"+atom0.getPosition());
  //          System.out.println("j+1:"+atom1.getPosition());
  //          System.out.println("j+2:"+atom2.getPosition());
  //          System.out.println("j+3:"+atom3.getPosition());

  //          System.out.println("oldcosphi:"+oldcosphi);
  //          System.out.println("dr21:"+dr21);
   //         System.out.println("dr34:"+dr34);
         

            if (oldcosphi > 0.99999 || oldcosphi < -0.99999) {
                // torsion angle is almost exactly 0 or PI.  skip this molecule
                for (int k=0; k<numChildren; k++) {
                    // save positions anyway since rejectNotify won't know to skip this
                    IAtom atomk = childList.getAtom(k);
//                    System.out.println("skipping "+i+" "+k+" "+atomk+" "+atomk.getPosition());
                    oldPositions[i][k].E(atomk.getPosition());
                }
                continue;
            }
            double oldphi = Math.acos(oldcosphi);
            int oldTorsionBin = probabilityReverseMap[(int)((oldphi/Math.PI)*(probabilityBins.length-1))];
//            System.out.println("oldTorsionBin0 = "+oldTorsionBin);
            if (probabilityBins[oldTorsionBin] > oldphi) {
                throw new RuntimeException("reversemap overshot "+oldphi+" "+oldTorsionBin+" "+probabilityBins[oldTorsionBin]);
            }
            while (oldTorsionBin < probabilityReverseMap.length && probabilityBins[oldTorsionBin+1] < oldphi) {
                oldTorsionBin++;
            }
            work1.E(dr21);
            work1.XE(dr34);
            // do we need to rotate the first part of the molecule clockwise or counter-clockwise?

            boolean flippit = work1.dot(dr23) > 0;
            // pick a new bin
            int torsionBin = random.nextInt(probabilityBins.length-1);
            // pick some angle within that bin
            double newphi = probabilityBins[torsionBin] + random.nextDouble()*(probabilityBins[torsionBin+1]-probabilityBins[torsionBin]);
            double deltaphi = 0.5*(newphi - oldphi);
            // if we move from a high-U (wide) bin to a low-U (narrow) bin, the nominal
            // acceptance probability will be too high.  correct for that here. 
            bias *= binSize[torsionBin]/binSize[oldTorsionBin];
            if (flippit) deltaphi = -deltaphi;
   //         System.out.println("new phi: "+ newphi);
            
   //         System.out.println("cos new phi: "+ Math.cos(newphi));

            /////////////////////////////////////////// update all atoms /////////////////////////////////////////////////////////
            for (int k=0; k<numChildren; k++) {// loop over all the atoms in the picked molecule, regardless it is carbon or hydrogen
  //          	System.out.println("k value is : " + k);
  //          	System.out.println("deltaphi value is : " + deltaphi);
            	
            	 IAtom atom = childList.getAtom(k);
                 Vector atomPosition = atom.getPosition();
                 oldPositions[i][k].E(atomPosition);
                 
                if ( (k==(j+1)) ||  (k==(j+2)) ){//1
                    continue;
                }
                if ( k==(j+3) ) {//  2
                    deltaphi = -deltaphi;
                }
                if ( k==numCarbons ){ //3
                    deltaphi = -deltaphi;
                }
                if ( k ==(j+2+numCarbons) ){//4
                    deltaphi = -deltaphi;
                }
                if ( k==(2*numCarbons) ){// 5
                    deltaphi = -deltaphi;
                }
                if ( k ==(j + 2 + 2*numCarbons) ){ //6
                    deltaphi = -deltaphi;
                }
                if ( k ==(3 * numCarbons)){//7
                    deltaphi = -deltaphi;
                }
                if ( k ==(3*numCarbons+1) ){//8
                    deltaphi = -deltaphi;
                }
  //          	System.out.println("deltaphi value is : " + deltaphi);
  //          	System.out.println("newphi"+newphi+ "    oldphi:"+oldphi);
                rotationTensor.setRotationAxis(axis, deltaphi);
                atomPosition.ME(atom1Position);//atom1Position: r(j+1)
                rotationTensor.transform(atomPosition);
                atomPosition.PE(atom1Position);
    
            }
       
            oldCenter.ME(positionDefinition.position(selectedMolecules.getMolecule(i)));// COM of i-th molecule
            for (int k=0; k<numChildren; k++) {
                // shift the whole molecule so that the center of mass (or whatever the position definition uses) doesn't change
            	IAtom atomk = childList.getAtom(k);
                atomk.getPosition().PE(oldCenter);
            }
        }
        ((BoxCluster)box).trialNotify();
        wNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);

        uNew = energyMeter.getDataAsScalar();
   //     System.out.println(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>after move");

//        System.out.println(uOld+" => "+uNew+"   "+wOld+" => "+wNew+" "+bias+" "+stepSize);
        return true;
    }
    
     public void setBondLength(double b) {
        bondLength = b;
    }
    protected void selectMolecules() {
        IMoleculeList molecules = box.getMoleculeList();
        selectedMolecules = new MoleculeArrayList();
        oldPositions = new Vector[molecules.getMoleculeCount()][0];
        int i=0;
        for (int k=0; k < molecules.getMoleculeCount();k++) {
            IMolecule a = molecules.getMolecule(k);
            int numChildren = a.getChildList().getAtomCount();
            if (numChildren<14) { //at least C4H10, total atoms:14
                continue;
            }
            oldPositions[i] = new Vector[numChildren];
            for (int j=0; j<numChildren; j++) {
                oldPositions[i][j] = space.makeVector();
            }
            i++;
            selectedMolecules.add(a);
        } 
    }
    
    public void acceptNotify() {
    	((BoxCluster)box).acceptNotify();
    }

    public void rejectNotify() {
        for(int i=0; i<selectedMolecules.getMoleculeCount(); i++) {
            IAtomList childList = selectedMolecules.getMolecule(i).getChildList();
            for (int j=0; j<childList.getAtomCount(); j++) {
                IAtom atomj = childList.getAtom(j);
                atomj.getPosition().E(oldPositions[i][j]);// move every atom back
            }
//            System.out.println(selectedAtoms[i]+" rejected => "+selectedAtoms[i].coord.position());
        }
        ((BoxCluster)box).rejectNotify();
    }

    public double getB() {
        return -(uNew - uOld);
    }
 
    public double getA() {
    	return bias*wNew/wOld;
    }
    private static final long serialVersionUID = 1L;
    protected final MeterPotentialEnergy energyMeter;
    protected final P4BondTorsion torsionPotential;
    protected IMoleculePositionDefinition positionDefinition;
    protected final double[] probabilityBins;
    protected final double[] binSize;
    protected final int[] probabilityReverseMap;
    protected MoleculeArrayList selectedMolecules;
    protected double bondLength;
    protected final Vector work1, work2, work3;
    protected final Vector dr21, dr23, dr34;
    protected Vector[][] oldPositions;
    protected final Vector oldCenter;
    protected double wOld, wNew, bias;
    protected ISpecies species;

}
