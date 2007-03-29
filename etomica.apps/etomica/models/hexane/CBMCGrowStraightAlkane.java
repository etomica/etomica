package etomica.models.hexane;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.integrator.IntegratorMC;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.space.IVector;
import etomica.space.IVectorRandom;
import etomica.space3d.Vector3D;
import etomica.species.Species;
import etomica.util.IRandom;

/**
 * Grows a straight-chain alkane of specified length.
 * Follows Frenkel & Smit algorithm 25, 26, 27, 28 in Chapter 13.
 * 
 * Programmer will need to know about the potentials used by the model in order
 * to write the necessary subclass.
 * 
 * The following methods are functional in this code, but may overridden 
 * as needed to implement a specific model:
 * <ul>
 * <li>calcRandomBond()
 * <li>calcRandomBondWithAngle()
 * <li>calcRandomBondWithAngleAndTorsion()
 * <li>calcStartIndex()
 * </ul>
 * @author cribbin
 *
 */

public abstract class CBMCGrowStraightAlkane extends MCMoveCBMC {

    public CBMCGrowStraightAlkane(PotentialMaster potentialMaster, IRandom random, IntegratorMC integrator, Phase p, Species species, int n, int NTrials){
        super(potentialMaster, random, integrator, p, n, NTrials);
        
        setChainlength(n);
        sumW = 0.0;

        Species[] sp = new Species[1];
        sp[0] = species;

        vex = (IVectorRandom)potentialMaster.getSpace().makeVector();
        temp = potentialMaster.getSpace().makeVector();
        
        a = new double[numTrial];
        b = new double[chainlength];   //used to store old rosenbluth factors
    

        storePos = new IVector[numTrial];
        for(int k = 0; k < numTrial; k++){
            storePos[k] = potentialMaster.getSpace().makeVector();
        }
        
        
        tempCloser = potentialMaster.getSpace().makeVector();
        tempFarther = potentialMaster.getSpace().makeVector();
    }
    
    /**
     * Calculates wOld and wNew; does not have a return value.
     */
    protected void calcRosenbluthFactors(){
        //Pick a direction and an atom to start with
        forward = random.nextInt(2) == 0;
        int dir, endIndex, beginIndex;
        if(forward){
            dir = 1;
            endIndex = chainlength;
            beginIndex = 0;
        } else {
            dir = -1;
            endIndex = -1;
            beginIndex = chainlength-1;
        }
        int startIndex = calcStartIndex();
        //nan
        
        double uExt;
        wOld = 1.0;
        wNew = 1.0;
        double sumA = 0.0;;
        
    
//OLD OLD OLD OLD OLD OLD OLD
        //Calculate the OLD Rosenbluth factor
        //the for loops in these if statements do not calculate using the
        // startIndex.  The startIndex value is calculated in the code that follows
        // the loops.  No loop is needed for trials, as trials are not occurring.
        for(int i = startIndex; i != endIndex; i += dir){//This loops through the atoms
            sumA = 0.0;
            for(int k = 0; k < numTrial-1; k++){  //This loops through the trials
                
                if(i == beginIndex){    //If we're placing the first atom of a molecule
                    (((AtomLeaf)atomList.get(i)).getCoord().getPosition()).E(phase.getBoundary().randomPosition());
                } else if(i == beginIndex + dir){  //If we're placing the second atom of a molecule
                    ((AtomLeaf)atomList.get(i)).getCoord().getPosition().E(calcRandomBond());
                    ((AtomLeaf)atomList.get(i)).getCoord().getPosition().PE(((AtomLeaf)atomList.get(i-dir)).getCoord().getPosition());
                } else if(i == beginIndex + dir * 2){//If we're placing the third atom of a molecule
                    ((AtomLeaf)atomList.get(i)).getCoord().getPosition().E(calcRandomBondWithAngle(
                            (AtomLeaf)atomList.get(i-dir),
                            (AtomLeaf)atomList.get(i-2*dir)));
                    ((AtomLeaf)atomList.get(i)).getCoord().getPosition().PE(((AtomLeaf)atomList.get(i-dir)).getCoord().getPosition());                 
                } else {//For the rest of the atoms in a molecule
                    ((AtomLeaf)atomList.get(i)).getCoord().getPosition().E(calcRandomBondWithAngleAndTorsion(
                            (AtomLeaf)atomList.get(i-dir), 
                            (AtomLeaf)atomList.get(i-2*dir), 
                            (AtomLeaf)atomList.get(i-3*dir)));
                    ((AtomLeaf)atomList.get(i)).getCoord().getPosition().PE(((AtomLeaf)atomList.get(i-dir)).getCoord().getPosition());
                }
                
                //evaluate the Boltzmann factor of this configuration
                // (configuration of this molecule, for this trial)
                // and store it.
                uExt = calcExternalEnergy(((AtomLeaf)atomList.get(i)));
                if(i == endIndex || i == 0){
                    a[k] = numTrial * Math.exp(-beta*uExt);
                } else {
                    a[k] = Math.exp(-beta*uExt);
                }
                sumA += a[k];
            }//end of k loop
            
            //do the k-loop stuff for the actual position of the molecule, since we are in the old section
            ((AtomLeaf)atomList.get(i)).getCoord().getPosition().E(positionOld[i]);
            uExt = calcExternalEnergy(((AtomLeaf)atomList.get(i)));
            if(i == endIndex || i == 0){
                a[numTrial-1] = numTrial * Math.exp(-beta*uExt);
            } else {
                a[numTrial-1] = Math.exp(-beta*uExt);
            }
            sumA += a[numTrial-1];
            
//          Calculate the probablilities
            for(int k = 0; k < numTrial; k++){
                a[k] /= sumA;
            }

            //Increment the Rosenbluth factor for the system.
            wOld *= sumA;

        }//end of i loop
        
        
        
//NEW NEW NEW NEW NEW NEW NEW     
        //Calculate the NEW Rosenbluth factor

        for(int i = startIndex; i != endIndex; i += dir){//This loops through the atoms
            sumA = 0.0;
            for(int k = 0; k < numTrial; k++){  //This loops through the trials
                
                if(i == beginIndex){    //If we're placing the first atom of a molecule
                    (((AtomLeaf)atomList.get(i)).getCoord().getPosition()).E(phase.getBoundary().randomPosition());
                } else if(i == beginIndex + dir){  //If we're placing the second atom of a molecule
                    ((AtomLeaf)atomList.get(i)).getCoord().getPosition().E(calcRandomBond());
                    ((AtomLeaf)atomList.get(i)).getCoord().getPosition().PE(((AtomLeaf)atomList.get(i-dir)).getCoord().getPosition());
                } else if(i == beginIndex + dir * 2){//If we're placing the third atom of a molecule
                    ((AtomLeaf)atomList.get(i)).getCoord().getPosition().E(calcRandomBondWithAngle(
                            (AtomLeaf)atomList.get(i-dir),
                            (AtomLeaf)atomList.get(i-2*dir)));
                    ((AtomLeaf)atomList.get(i)).getCoord().getPosition().PE(((AtomLeaf)atomList.get(i-dir)).getCoord().getPosition());                 
                } else {//For the rest of the atoms in a molecule
                    ((AtomLeaf)atomList.get(i)).getCoord().getPosition().E(calcRandomBondWithAngleAndTorsion(
                            (AtomLeaf)atomList.get(i-dir), 
                            (AtomLeaf)atomList.get(i-2*dir), 
                            (AtomLeaf)atomList.get(i-3*dir)));
                    ((AtomLeaf)atomList.get(i)).getCoord().getPosition().PE(((AtomLeaf)atomList.get(i-dir)).getCoord().getPosition());
                }
                
//              store new position
                storePos[k].E(((AtomLeaf)atomList.get(i)).getCoord().getPosition());
                
                //evaluate the Boltzmann factor of this configuration
                // (configuration of this molecule, for this trial)
                // and store it.
                uExt = calcExternalEnergy(((AtomLeaf)atomList.get(i)));  
                if(i == endIndex || i == 0){
                    a[k] = numTrial * Math.exp(-beta*uExt);
                } else {
                    a[k] = Math.exp(-beta*uExt);
                }
//                if(i == 0)  {a[k] *= getPrefactor();}
                sumA += a[k];
                
//                System.out.println("Trial# " + k);
//                System.out.println("a[k] " + a[k]);
            }//end of k loop
            
//            System.out.println("SumA = " + sumA);
            //What to do if none of the trials was accepted.
            if(sumA == 0.0){
                wNew = 0.0;
//                System.out.println("Bailing from CBMCGrowStraightAlkane.calcRosenbluthFactor() now!");
                return;
            }
           
//          Calculate the probablilities
            for(int k = 0; k < numTrial; k++){
                a[k] /= sumA;
            }
            
            //Per discussion with Andrew on 3/26/07 & Algorithm 41 in F&S p.577)
            double rand = random.nextDouble();
            double sum = 0.0;
            int pickThisOne = numTrial-1;
            
            for(int j = 0; j < a.length; j++){
                sum += a[j];
                if(rand < sum){ 
                    pickThisOne = j;
                    break;
                }
            }
            
            System.out.println("I picked trial " + pickThisOne);
            //Move the atom to the selected position
            ((AtomLeaf)atomList.get(i)).getCoord().getPosition().E(storePos[pickThisOne]);
            
            //Increment the Rosenbluth factor for the system.
            wNew *= sumA;
        }//end of i loop
    }

    protected int calcStartIndex(){
        return random.nextInt(chainlength);
    }

    /**
     * Returns the bond length for the new atom in a test
     * Should take bonded potential into account.
     * @return a new bond length
     */
    protected abstract double calcBondL();
    
    /**
     * Generates a random bond.
     * @return a new bond vector
     */
    protected IVector calcRandomBond(){
        vex.setRandomSphere(random);
        vex.TE(calcBondL());
        return vex;
    }
    
    /**
     * XXX METHOD NOT TESTED!  USE AT OWN PERIL!
     * Generates a random bond that accounts for bond angle in the calculation
     * This method will change its parameter.
     * @param b the atom farther from the new atom
     * @param a the atom nearer to the new atom
     * @return a new bond vector
     */
    //Based on algorithm 45 in Frenkel & Smit
    protected IVector calcRandomBondWithAngle(AtomLeaf a, AtomLeaf b){
        double phi;
        double ubb;

        tempCloser.E(a.getCoord().getPosition());
        tempCloser.ME(b.getCoord().getPosition());
        tempCloser.normalize();
        
        do{
           vex.setRandomSphere(random);
           phi = Math.acos(vex.dot(tempCloser));
           ubb = calcBondAngleEnergy(phi);
        } while(random.nextDouble() < Math.exp(-beta*ubb));
        
        vex.TE(calcBondL());
        return vex;
    }
    
    /**
     * XXX METHOD NOT TESTED!  USE AT OWN PERIL!l dy
     * Generates a random bond that accounts for bond angle and torsion in 
     * the calculation.
     * This method will change its parameters.
     * @param c the atom farthest from the new atom
     * @param b the atom nearer to the new atom
     * @param a the atom nearest to the new atom
     * @return a new bond vector
     */
    //Based on algorithm 46 in Frenkel & Smit
    protected IVector calcRandomBondWithAngleAndTorsion(AtomLeaf a, AtomLeaf b,
            AtomLeaf c){
        if(phase.getSpace().D() != 3){
            throw new IllegalArgumentException("Torsional bond is only used in 3D simulations");
        }

        tempFarther.E(b.getCoord().getPosition());
        tempFarther.ME(c.getCoord().getPosition()); 
        tempCloser.E(a.getCoord().getPosition());
        tempCloser.ME(b.getCoord().getPosition());
        double phi, theta;
        double ubb, utors, usum;
        
        tempFarther.normalize();
        tempCloser.normalize();
        
        do{
            vex.setRandomSphere(random);
            temp.E(vex);

            phi = Math.acos(temp.dot(tempCloser));
            ubb = calcBondAngleEnergy(phi);
            //nan put if ubb is larger than an exceptionally large ## in here, bail! (continue)
            
            temp.E(vex);
            ((Vector3D)temp).XE((Vector3D)tempCloser);
            ((Vector3D)tempCloser).XE((Vector3D)tempFarther);
            theta = temp.dot(tempCloser);
            
            //nan pass all info to this method or calcBTE depends on theta only & pass only theta in.
            //nan this dude needs an actual length
            //nan calc BTE needs actual vector from 1-4, not vector from 3-4
            utors = calcBondTorsionalEnergy(vex);
            
            usum = utors + ubb;
            
        } while (random.nextDouble() < Math.exp(-beta*usum));
        
        vex.TE(calcBondL());
        return vex;
    }
    
    protected abstract double calcExternalEnergy(Atom a);
    
    //uub in algorithm 45, 46
    protected abstract double calcBondAngleEnergy(double dub);
    
    //utors in algorithm 46
    protected abstract double calcBondTorsionalEnergy(IVector v);
    
    public abstract double energyChange();

    
    boolean forward;    //indicates direction of growth
    protected double bondlength;
    double sumW;
    IVector temp;
    IVectorRandom vex;
    double[] a;
    double[] b;
    IVector[] storePos;
    IVector tempCloser;
    IVector tempFarther;
}
