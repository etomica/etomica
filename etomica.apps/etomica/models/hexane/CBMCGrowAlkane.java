package etomica.models.hexane;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorMolecule;
import etomica.atom.iterator.IteratorDirective.Direction;
import etomica.integrator.IntegratorMC;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.species.Species;

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

public abstract class CBMCGrowAlkane extends MCMoveCBMC {

    public CBMCGrowAlkane(PotentialMaster potentialMaster, IntegratorMC integrator, Species species){
        super(potentialMaster, integrator);
        
        sumW = 0.0;
        setPrefactor(1.0);

        Species[] sp = new Species[1];
        sp[0] = species;
        affectedAtomIterator = new AtomIteratorMolecule(sp);
        affectedAtomIterator.setPhase(phase);
        
    }
    
    /**
     * Calculates wOld and wNew; does not have a return value.
     */
    protected void calcRosenbluthFactors(){
        
        //Pick a direction and an atom to start with
        forward = Simulation.random.nextBoolean();
        int startIndex = calcStartIndex();
        int numTrial = calcNumberOfTrials() + 1;
        affectedAtomIterator.setTarget(atomList.get(startIndex));
        if(forward){
            affectedAtomIterator.setDirection(Direction.UP);
        } else {
            affectedAtomIterator.setDirection(Direction.DOWN);
        }
        
        Vector vex = phase.space().makeVector();
        double uExt;
        double[] a = new double[numTrial];
        wOld = 1.0;
        wNew = 1.0;
        double sumA = 0.0;;
        double[] b = new double[chainlength];   //used to store old rosenbluth factors
    
//OLD OLD OLD OLD OLD OLD OLD
        //Calculate the OLD Rosenbluth factor
        //the for loops in these if statements do not calculate using the
        // startIndex.  The startIndex value is calculated in the code that follows
        // the loops.  No loop is needed for trials, as trials are not occurring.
        sumA = 0.0;
        if(forward){
            for(int i = chainlength-1; i >= startIndex; i--){
                uExt = externalMeter.getDataAsScalar();               
                b[i] = Math.exp(-beta*uExt);
                if(i == 0)  {b[i] *= getPrefactor();}
                sumA += b[i];
                wOld *= sumA;
            }// end of i loop
            
        } else {
            for(int i = 0; i >= startIndex; i++){
                uExt = externalMeter.getDataAsScalar();
                b[i] = Math.exp(-beta*uExt);
                if(i == chainlength - 1)  {b[i] *= getPrefactor();}
                sumA += b[i];
                wOld *= sumA;
            }//end of i loop
        }//end of forward/reverse if-else statement
        
        
        
//NEW NEW NEW NEW NEW NEW NEW     
        //Calculate the NEW Rosenbluth factor
        Vector[] storePos = new Vector[numTrial];
        for(int k = 0; k < numTrial; k++){
            storePos[k] = phase.space().makeVector();
        }
        
        sumA = 0.0;
        if(forward){
            for(int i = startIndex; i < chainlength; i++){//This loops through the atoms
                for(int k = 0; k < numTrial; k++){  //This loops through the trials
                    bondlength = calcBondL(); //possibly generated
                    
                    if (i == 0) { //if we're starting with methane
                        vex.setRandomSphere();
                        (((AtomLeaf)atomList.get(0)).getCoord().getPosition()).E(vex);
                    } else if(i == 1) { //ethane
                        vex.E(calcRandomBond());
                        vex.PE(((AtomLeaf)atomList.get(0)).getCoord().getPosition());
                        ((AtomLeaf)atomList.get(1)).getCoord().getPosition().E(vex);
                    } else if(i == 2) { //propane
                        Vector temp = phase.space().makeVector();
                        temp.E(((AtomLeaf)atomList.get(1)).getCoord().getPosition());
                        temp.ME(((AtomLeaf)atomList.get(0)).getCoord().getPosition());        
                        vex.E(calcRandomBondWithAngle(temp));
                        vex.PE(((AtomLeaf)atomList.get(1)).getCoord().getPosition());
                        ((AtomLeaf)atomList.get(2)).getCoord().getPosition().E(vex);
                    } else {  //butane and up; has a torsional thing going.
                        vex.E(calcRandomBondWithAngleAndTorsion(
                                (AtomLeaf)atomList.get(i-1), 
                                (AtomLeaf)atomList.get(i-2), 
                                (AtomLeaf)atomList.get(i-3)));
                        vex.PE(((AtomLeaf)atomList.get(i-1)).getCoord().getPosition());
                        ((AtomLeaf)atomList.get(i)).getCoord().getPosition().E(vex);
                    }
                    
                    //store new position
                    storePos[k].E(((AtomLeaf)atomList.get(i)).getCoord().getPosition());
                    
                    //evaluate the Boltzmann factor of this configuration
                    // (configuration of this molecule, for this trial)
                    // and store it.
                    uExt = externalMeter.getDataAsScalar();
                    a[k] = Math.exp(-beta*uExt);
                    if(i == 0)  {a[k] *= getPrefactor();}
                    sumA += a[k];
                }//end of k loop
                
                //Calculate the probablilities
                for(int k = 0; k < numTrial; k++){
                    a[k] /= sumA;
                }
                
                //Select the best (most probable move) 
                int best = 0;
                for(int k = 1; k < numTrial; k++){
                    if(a[k] > a[best]){
                        best = k;
                    }
                }
                
                //Move the atom to the selected position
                ((AtomLeaf)atomList.get(i)).getCoord().getPosition().E(storePos[best]);
                
                //Increment the Rosenbluth factor for the system.
                wNew *= sumA;
            }//end of i loop

        } else {
            for(int i = startIndex; i >= 0; i--){
                for(int k = 0; k < numTrial; k++){
                    bondlength = calcBondL();
                    
                    if (i == chainlength-1) { //if we're starting with methane
                        vex.setRandomSphere();
                        (((AtomLeaf)atomList.get(chainlength-1)).getCoord().getPosition()).E(vex);
                    } else if(i == chainlength - 2) {  //ethane
                        vex.E(calcRandomBond());
                        vex.PE(((AtomLeaf)atomList.get(chainlength-1)).getCoord().getPosition());
                        ((AtomLeaf)atomList.get(chainlength-2)).getCoord().getPosition().E(vex);
                    } else if(i == chainlength - 3) {  //propane
                        Vector temp = phase.space().makeVector();
                        temp.E(((AtomLeaf)atomList.get(chainlength-1)).getCoord().getPosition());
                        temp.ME(((AtomLeaf)atomList.get(chainlength-2)).getCoord().getPosition());        
                        vex.E(calcRandomBondWithAngle(temp));
                        vex.PE(((AtomLeaf)atomList.get(chainlength-2)).getCoord().getPosition());
                        ((AtomLeaf)atomList.get(chainlength-3)).getCoord().getPosition().E(vex);
                    } else {    //butane and up
                        vex.E(calcRandomBondWithAngleAndTorsion(
                                (AtomLeaf)atomList.get(i+1),
                                (AtomLeaf)atomList.get(i+2),
                                (AtomLeaf)atomList.get(i+3)));
                        vex.PE(((AtomLeaf)atomList.get(i+1)).getCoord().getPosition());
                        ((AtomLeaf)atomList.get(i)).getCoord().getPosition().E(vex);
                    }
                    
//                  store new position
                    storePos[k].E(((AtomLeaf)atomList.get(i)).getCoord().getPosition());
                    
                    //evaluate the Boltzmann factor of this configuration
                    // (configuration of this molecule, for this trial)
                    // and store it.
                    uExt = externalMeter.getDataAsScalar();
                    a[k] = Math.exp(-beta*uExt);
                    if(i == chainlength-1) {a[k] *= getPrefactor();}
                    sumA += a[k];
                }//end of k loop
                
                //Calculate the probablilities
                for(int k = 0; k < numTrial; k++){
                    a[k] /= sumA;
                }
                
                //Select the best (most probable move) 
                int best = 0;
                for(int k = 1; k < numTrial; k++){
                    if(a[k] > a[best]){
                        best = k;
                    }
                }
                
                //Move the atom to the selected position
                ((AtomLeaf)atomList.get(i)).getCoord().getPosition().E(storePos[best]);
                
                //Increment the Rosenbluth factor for the system.
                wNew *= sumA;
                
            }//end of i loop
        }//end of forward/reverse if-else statement
    }

    public double getA(){
        return 1.0;
    }
    
    public double getB(){
        uNew = energyMeter.getDataAsScalar();
        return -(uNew - uOld);
    }
    
    public void setPrefactor(double p){
        pref = p;
    }
    public double getPrefactor(){
        return pref;
    }
    protected int calcStartIndex(){
        return Simulation.random.nextInt(chainlength);
    }

    protected abstract int calcNumberOfTrials();
    
    /**
     * Returns the bond length for the new atom in a test
     * Should take bonded potential into account.
     * @return a new bond length
     */
    protected abstract double calcBondL();
    
    /**
     * Generates a random bond on a unit sphere
     * @return a new bond vector
     */
    protected Vector calcRandomBond(){
        Vector vax = phase.space().makeVector();
        vax.setRandomSphere();
        vax.TE(calcBondL());
        return vax;
    }
    
    /**
     * Generates a random bond that accounts for bond angle in the calculation
     * This method will change its parameter.
     * @param v a vector that represents the prior bond
     * @return a new bond vector
     */
    //Based on algorithm 45 in Frenkel & Smit
    protected Vector calcRandomBondWithAngle(Vector v){
        Vector vax = phase.space().makeVector();
        double phi;
        double ubb;
        v.normalize();
        
        boolean ready = false;
        do {
           vax.setRandomSphere();
           Vector vux = phase.space().makeVector();
           vux.E(vax);
           
           phi = Math.acos(vux.dot(v));
           ubb = calcBondBendingEnergy(phi);
           if(Simulation.random.nextDouble() < Math.exp(-beta*ubb)) {ready = true;}
        } while (ready == false);
        
        vax.TE(calcBondL());
        return vax;
    }
    
    /**
     * Generates a random bond that accounts for bond angle and torsion in 
     * the calculation.
     * This method will change its parameters.
     * @param c the atom farthest to the new atom
     * @param b the atom nearer to the new atom
     * @param a the atom nearest to the new atom
     * @return a new bond vector
     */
    //Based on algorithm 46 in Frenkel & Smit
    protected Vector calcRandomBondWithAngleAndTorsion(AtomLeaf a, AtomLeaf b,
            AtomLeaf c){
        if(phase.space().D() != 3){
            throw new IllegalArgumentException("Torsional bond is only used in 3D simulations");
        }
        
        Vector vax = phase.space().makeVector();
        Vector vux = phase.space().makeVector();
        
        Vector tempCloser = phase.space().makeVector();
        Vector tempFarther = phase.space().makeVector();
        tempFarther.E(b.getCoord().getPosition());
        tempFarther.ME(c.getCoord().getPosition()); 
        tempCloser.E(a.getCoord().getPosition());
        tempCloser.ME(b.getCoord().getPosition());
        double phi, theta;
        double ubb, utors, usum;
        
        tempFarther.normalize();
        tempCloser.normalize();
        
        boolean ready = false;
        do{
            vax.setRandomSphere();
            vux.E(vax);

            phi = Math.acos(vux.dot(tempCloser));
            ubb = calcBondBendingEnergy(phi);
            
            vux.E(vax);
            ((Vector3D)vux).XE((Vector3D)tempCloser);
            ((Vector3D)tempCloser).XE((Vector3D)tempFarther);
            theta = vux.dot(tempCloser);
            
            utors = calcBondTorsionalEnergy(vax);
            
            usum = utors + ubb;
            if(Simulation.random.nextDouble() < Math.exp(-beta*usum)) {
                ready = true;
            }
            
        } while (ready == false);
        
        vax.TE(calcBondL());
        return vax;
    }
    
    protected abstract double calcExternalEnergy(Atom a);
    
    //uub in algorithm 45, 46
    protected abstract double calcBondBendingEnergy(double dub);
    
    //utors in algorithm 46
    protected abstract double calcBondTorsionalEnergy(Vector v);
    
    public AtomIterator affectedAtoms(){
        return affectedAtomIterator;
    }
    
    boolean forward;    //indicates direction of growth
    protected double bondlength;
    protected double pref;
    
    double sumW;
    private AtomIteratorMolecule affectedAtomIterator;
}
