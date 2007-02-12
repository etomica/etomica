package etomica.models.hexane;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.integrator.IntegratorMC;
import etomica.potential.Potential;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.species.Species;

/**
 * Contains specifics for the solid hexane model:
 * <ul>
 * <li>bond length = 0.4
 * <li>rejects moving the whole molecule as a CBMC move
 * </ul>
 * @author cribbin
 *
 */
public class CBMCGrowSolidHexane extends CBMCGrowAlkane {

    public CBMCGrowSolidHexane(PotentialMaster p, IntegratorMC integrator, Species species){
        super(p, integrator, species);
        setBondLength(0.4);
        setPrefactor(1.0);
        setNumberOfTrials(20);  //number picked randomly/ off F&S's reference.
        setChainlength(6);
        
        phi = 109.47 / 360.0 * 2 * Math.PI;
        
        pots = p.getPotentials();
        potMast = p;
        }
    
    protected int calcStartAtom(){
        int si;
        do{
            si = Simulation.random.nextInt(chainlength);//function will not pick chainlength
        } while (si != 0);
        return si;
    }
    
    //Different because we know the bond angle
    //All moves are accepted,
    protected Vector calcRandomBondWithAngle(Vector v){
        Vector vax = phase.space().makeVector();
        double ubb;
        v.normalize();
        
        boolean ready = false;
        do {
           double ang = Simulation.random.nextDouble();
           ang *= 2*Math.PI;
           ubb = calcBondBendingEnergy(phi);
           if(Simulation.random.nextDouble() < Math.exp(-beta*ubb)) {ready = true;}
        } while (ready == false);
        
        vax.TE(getBondlength());
        return vax;
    }
    
    protected Vector calcRandomBondWithAngleAndTorsion(AtomLeaf a, AtomLeaf b, 
            AtomLeaf c){
        if(phase.space().D() != 3){
            throw new IllegalArgumentException("Torsional bond is only used in 3D simulations");
        }
        
        Vector vax = phase.space().makeVector();
        Vector vux = phase.space().makeVector();
        double theta;
        double ubb, utors, usum;
        Vector tempCloser = phase.space().makeVector();
        Vector tempFarther = phase.space().makeVector();
        
        tempFarther.E(b.getCoord().position());
        tempFarther.ME(c.getCoord().position()); 
        tempCloser.E(a.getCoord().position());
        tempCloser.ME(b.getCoord().position());
     
        boolean ready = false;
        /*
         * This bit here takes care of all the pairs that the internal energy
         * needs to have called.  We aren't really calling pairs, we're just 
         * constructing the system to comply with what we need as far as the 
         * placement of the new atom in the molecule.  We do have a hard
         * sphere system with rigid bond lengths and set bond angles, so we can
         * "get away" with doing this.
         * 
         */
        do{
            vax = calcRandomBondWithAngle(tempCloser);
            vux.E(vax);
            vux.PE(tempCloser);
            vux.PE(tempFarther);
            utors = calcBondTorsionalEnergy(vux);
            if(utors == 0.0){ready = true;}
        } while (ready == false);
        
        ubb = calcBondBendingEnergy(phi);
        usum = ubb + utors;
        vax.PE(a.getCoord().position());
        
        return vax;
    }
    
    protected double calcBondBendingEnergy(double d){
        /*
         * We are using a set bond length, and a rigid bond angle,
         * we don't need to calculate an energy because we wrote this program
         * to only use the proper bond length and angle
         */
        return 0.0;
    }
    
    protected double calcBondTorsionalEnergy(Vector v){
        double denom = Math.sqrt(v.x(0) * v.x(0) + v.x(1) * v.x(1) + v.x(2) * v.x(2));
        
        if (denom > 1.0){
            return 0.0;
        }
        
        return 1.0;
    }
    
    protected double calcExternalEnergy(Atom a){
        //make sure that the proper potential is enabled.  Really,
        //we only have the one potential, so this line is unnecessary
        //but I want it in here for reference when I am extending
        //this code.
        potMast.setEnabled(pots[0], true);
        externalMeter.setPhase(phase);
        externalMeter.setTarget(a);
        
        return externalMeter.getDataAsScalar();
    }
    //Since the bondlength is constant, this returns an ordinary getter
    protected double calcBondL(){
        return getBondlength();
    }
    public void setBondLength(double d){
        bondlength = d;
    }
    public double getBondlength(){
        return bondlength;
    }
    //I didn't do any work to calculate the number of trials.
    public int calcNumberOfTrials(){
        return getNumberOfTrials();
    }
    public void setNumberOfTrials(int in){
        numTrial = in;
    }
    public int getNumberOfTrials(){
        return numTrial;
    }
    
    private static final long serialVersionUID = 1L;
    private Potential[] pots;
    private PotentialMaster potMast;
    private double phi;
}