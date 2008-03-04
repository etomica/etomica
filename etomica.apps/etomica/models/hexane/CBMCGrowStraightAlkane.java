package etomica.models.hexane;

import etomica.api.IBox;

import etomica.api.IAtom;
import etomica.api.IAtomPositioned;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.integrator.IntegratorMC;
import etomica.space.IVectorRandom;
import etomica.space3d.Vector3D;


/**
 * Grows a straight-chain alkane of specified length. Follows Frenkel & Smit
 * algorithm 25, 26, 27, 28 in Chapter 13.
 * 
 * Programmer will need to know about the potentials used by the model in order
 * to write the necessary subclass.
 * 
 * The following methods are functional in this code, but may overridden as
 * needed to implement a specific model:
 * <ul>
 * <li>calcRandomBond()
 * <li>calcRandomBondWithAngle()
 * <li>calcRandomBondWithAngleAndTorsion()
 * <li>calcStartIndex()
 * </ul>
 * 
 * @author cribbin
 * 
 */

public abstract class CBMCGrowStraightAlkane extends MCMoveCBMC {

    public CBMCGrowStraightAlkane(IPotentialMaster potentialMaster,
            IRandom random, IntegratorMC integrator, IBox p, ISpecies species,
            int n, int NTrials) {
        super(potentialMaster, random, integrator, p, n, NTrials);

        setChainlength(n);
        sumW = 0.0;

        ISpecies[] sp = new ISpecies[1];
        sp[0] = species;

        vex = (IVectorRandom) potentialMaster.getSpace().makeVector();
        temp = potentialMaster.getSpace().makeVector();
        a = new double[numTrial];
        b = new double[chainlength]; // used to store old rosenbluth factors
        storePos = new IVector[numTrial];
        // angleSet = new double[numTrial];
        for (int k = 0; k < numTrial; k++) {
            storePos[k] = potentialMaster.getSpace().makeVector();
        }
        tempCloser = potentialMaster.getSpace().makeVector();
        tempFarther = potentialMaster.getSpace().makeVector();

    }

    /**
     * Calculates wOld and wNew
     */
    protected boolean calcRosenbluthFactors() {
        // Pick a direction and an atom to start with
        forward = random.nextInt(2) == 0;
        int dir, endIndex, beginIndex;
        if (forward) {
            dir = 1;
            endIndex = chainlength;
            beginIndex = 0;
        } else {
            dir = -1;
            endIndex = -1;
            beginIndex = chainlength - 1;
        }
        int startIndex = calcStartIndex();

        double uExt;
        wOld = 1.0;
        wNew = 1.0;
        double sumA = 0.0;

        // System.out.println("Direction " + dir + " startIndex "+ startIndex);

        // NEW NEW NEW NEW NEW NEW NEW
        // Calculate the NEW Rosenbluth factor

        for (int i = startIndex; i != endIndex; i += dir) {// This loops
            // through the atoms
            sumA = 0.0;
            for (int k = 0; k < numTrial; k++) { // This loops through the
                // trials
                if (i == beginIndex) { // If we're placing the first atom of a
                    // molecule
                    (((IAtomPositioned) atomList.getAtom(i)).getPosition()).E(box
                            .getBoundary().randomPosition());
                } else if (i == beginIndex + dir) { // If we're placing the
                    // second atom of a molecule
                    ((IAtomPositioned) atomList.getAtom(i)).getPosition().E(
                            calcRandomBond());
                    ((IAtomPositioned) atomList.getAtom(i)).getPosition().PE(
                            ((IAtomPositioned) atomList.getAtom(i - dir)).getPosition());
                } else if (i == beginIndex + dir * 2) {// If we're placing the
                    // third atom of a
                    // molecule
                    ((IAtomPositioned) atomList.getAtom(i)).getPosition().E(
                            calcRandomBondWithAngle((IAtomPositioned) atomList.getAtom(i
                                    - dir), (IAtomPositioned) atomList
                                    .getAtom(i - 2 * dir)));
                    ((IAtomPositioned) atomList.getAtom(i)).getPosition().PE(
                            ((IAtomPositioned) atomList.getAtom(i - dir)).getPosition());
                } else {// For the rest of the atoms in a molecule
                    ((IAtomPositioned) atomList.getAtom(i)).getPosition().E(
                            calcRandomBondWithAngleAndTorsion(
                                    (IAtomPositioned) atomList.getAtom(i - dir),
                                    (IAtomPositioned) atomList.getAtom(i - 2 * dir),
                                    (IAtomPositioned) atomList.getAtom(i - 3 * dir)));
                    ((IAtomPositioned) atomList.getAtom(i)).getPosition().PE(
                            ((IAtomPositioned) atomList.getAtom(i - dir)).getPosition());
                }

                // store new position & angle
                storePos[k].E(((IAtomPositioned) atomList.getAtom(i)).getPosition());

                // evaluate the Boltzmann factor of this configuration
                // (configuration of this molecule, for this trial)
                // and store it.
                uExt = calcExternalEnergy(atomList.getAtom(i));
                if (i == endIndex || i == 0) {
                    a[k] = numTrial * Math.exp(-beta * uExt);
                } else {
                    a[k] = Math.exp(-beta * uExt);
                }
                sumA += a[k];
            }// end of k loop

            // What to do if none of the trials was accepted.
            if (sumA == 0.0) {
                wNew = 0.0;
                // System.out.println("Bailing from
                // CBMCGrowStraightAlkane.calcRosenbluthFactor() now!");
                rejectNotify();
                return false;
            }

            // Calculate the probablilities
            for (int k = 0; k < numTrial; k++) {
                a[k] /= sumA;
            }

            // Per discussion with Andrew on 3/26/07 & Algorithm 41 in F&S
            // p.577)
            double rand = random.nextDouble();
            double sum = 0.0;
            int pickThisOne = numTrial - 1;

            for (int j = 0; j < a.length; j++) {
                sum += a[j];
                if (rand < sum) {
                    pickThisOne = j;
                    break;
                }
            }

            // Move the atom to the selected position
            ((IAtomPositioned) atomList.getAtom(i)).getPosition().E(storePos[pickThisOne]);

            // System.out.println("New pos: "+
            // ((AtomLeaf)atomList.get(i)).getPosition());

            // Increment the Rosenbluth factor for the system.
            wNew *= sumA;

        }// end of i loop

        // OLD OLD OLD OLD OLD OLD OLD
        // Calculate the OLD Rosenbluth factor
        // the for loops in these if statements do not calculate using the
        // startIndex. The startIndex value is calculated in the code that
        // follows the loops. No loop is needed for trials, as trials are
        // not occurring.
        for (int i = startIndex; i != endIndex; i += dir) {// This loops
            // through the atoms
            sumA = 0.0;
            for (int k = 0; k < numTrial - 1; k++) { // This loops through
                // the trials

                if (i == beginIndex) { // If we're placing the first atom of a
                    // molecule
                    (((IAtomPositioned) atomList.getAtom(i)).getPosition()).E(box
                            .getBoundary().randomPosition());
                } else if (i == beginIndex + dir) { // If we're placing the
                    // second atom of a molecule
                    ((IAtomPositioned) atomList.getAtom(i)).getPosition().E(
                            calcRandomBond());
                    ((IAtomPositioned) atomList.getAtom(i)).getPosition().PE(
                            ((IAtomPositioned) atomList.getAtom(i - dir)).getPosition());
                } else if (i == beginIndex + dir * 2) {// If we're placing the
                    // third atom of a
                    // molecule
                    ((IAtomPositioned) atomList.getAtom(i)).getPosition().E(
                            calcRandomBondWithAngle((IAtomPositioned) atomList.getAtom(i
                                    - dir), (IAtomPositioned) atomList
                                    .getAtom(i - 2 * dir)));
                    ((IAtomPositioned) atomList.getAtom(i)).getPosition().PE(
                            ((IAtomPositioned) atomList.getAtom(i - dir)).getPosition());
                } else {// For the rest of the atoms in a molecule
                    ((IAtomPositioned) atomList.getAtom(i)).getPosition().E(
                            calcRandomBondWithAngleAndTorsion(
                                    (IAtomPositioned) atomList.getAtom(i - dir),
                                    (IAtomPositioned) atomList.getAtom(i - 2 * dir),
                                    (IAtomPositioned) atomList.getAtom(i - 3 * dir)));
                    ((IAtomPositioned) atomList.getAtom(i)).getPosition().PE(
                            ((IAtomPositioned) atomList.getAtom(i - dir)).getPosition());
                }

                // evaluate the Boltzmann factor of this configuration
                // (configuration of this molecule, for this trial)
                // and store it.
                uExt = calcExternalEnergy(atomList.getAtom(i));
                if (i == endIndex || i == 0) {
                    a[k] = numTrial * Math.exp(-beta * uExt);
                } else {
                    a[k] = Math.exp(-beta * uExt);
                }
                sumA += a[k];
            }// end of k loop

            // do the k-loop stuff for the actual position of the molecule,
            // since we are in the old section
            ((IAtomPositioned) atomList.getAtom(i)).getPosition().E(positionOld[i]);
            uExt = calcExternalEnergy(atomList.getAtom(i));
            if (i == endIndex || i == 0) {
                a[numTrial - 1] = numTrial * Math.exp(-beta * uExt);
            } else {
                a[numTrial - 1] = Math.exp(-beta * uExt);
            }
            sumA += a[numTrial - 1];

            // Calculate the probablilities
            for (int k = 0; k < numTrial; k++) {
                a[k] /= sumA;
            }

            // Increment the Rosenbluth factor for the system.
            wOld *= sumA;

        }// end of i loop
        
        return true;
    }

    protected int calcStartIndex() {
        return random.nextInt(chainlength);
    }

    /**
     * Returns the bond length for the new atom in a test Should take bonded
     * potential into account.
     * 
     * @return a new bond length
     */
    protected abstract double calcBondL();

    /**
     * Generates a random bond.
     * 
     * @return a new bond vector
     */
    protected IVector calcRandomBond() {
        vex.setRandomSphere(random);
        vex.TE(calcBondL());
        return vex;
    }

    /**
     * XXX METHOD NOT TESTED! USE AT OWN PERIL! Generates a random bond that
     * accounts for bond angle in the calculation This method will change its
     * parameter.
     * 
     * @param b
     *            the atom farther from the new atom
     * @param a
     *            the atom nearer to the new atom
     * @return a new bond vector
     */
    // Based on algorithm 45 in Frenkel & Smit
    protected IVector calcRandomBondWithAngle(IAtomPositioned a, IAtomPositioned b) {
        double phi;
        double ubb;

        tempCloser.E(a.getPosition());
        tempCloser.ME(b.getPosition());
        tempCloser.normalize();

        do {
            vex.setRandomSphere(random);
            phi = Math.acos(vex.dot(tempCloser));
            ubb = calcBondAngleEnergy(phi);
        } while (random.nextDouble() < Math.exp(-beta * ubb));

        vex.TE(calcBondL());
        return vex;
    }

    /**
     * XXX METHOD NOT TESTED! USE AT OWN PERIL!l dy Generates a random bond that
     * accounts for bond angle and torsion in the calculation. This method will
     * change its parameters.
     * 
     * @param c
     *            the atom farthest from the new atom
     * @param b
     *            the atom nearer to the new atom
     * @param a
     *            the atom nearest to the new atom
     * @return new bond vector
     */
    // Based on algorithm 46 in Frenkel & Smit
    protected IVector calcRandomBondWithAngleAndTorsion(IAtomPositioned a, IAtomPositioned b,
            IAtomPositioned c) {
/*
        if (dim != 3) {
            throw new IllegalArgumentException("Torsional bond is only used "
                    + "in 3D simulations");
        }
*/
        tempFarther.E(b.getPosition());
        tempFarther.ME(c.getPosition());
        tempCloser.E(a.getPosition());
        tempCloser.ME(b.getPosition());
        double phi, theta;
        double ubb, utors, usum;

        tempFarther.normalize();
        tempCloser.normalize();

        do {
            vex.setRandomSphere(random);
            temp.E(vex);

            phi = Math.acos(temp.dot(tempCloser));
            ubb = calcBondAngleEnergy(phi);
            // nan put if ubb is larger than an exceptionally large ## in here,
            // bail! (continue)

            temp.E(vex);
            ((Vector3D) temp).XE((Vector3D) tempCloser);
            ((Vector3D) tempCloser).XE((Vector3D) tempFarther);
            theta = temp.dot(tempCloser);

            // nan pass all info to this method or calcBTE depends on theta only
            // & pass only theta in.
            // nan this dude needs an actual length
            // nan calc BTE needs actual vector from 1-4, not vector from 3-4
            utors = calcBondTorsionalEnergy(vex);

            usum = utors + ubb;

        } while (random.nextDouble() < Math.exp(-beta * usum));

        vex.TE(calcBondL());
        return vex;
    }

    protected abstract double calcExternalEnergy(IAtom a);

    // uub in algorithm 45, 46
    protected abstract double calcBondAngleEnergy(double dub);

    // utors in algorithm 46
    protected abstract double calcBondTorsionalEnergy(IVector v);

    public abstract double energyChange();

    boolean forward; // indicates direction of growth

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
