/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.hexane;

import etomica.api.IAtom;
import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISpecies;
import etomica.space.Vector;
import etomica.box.RandomPositionSource;
import etomica.box.RandomPositionSourceRectangular;
import etomica.integrator.IntegratorMC;
import etomica.space.Space;


/**
 * Grows a straight-chain alkane of specified length. Follows Frenkel {@literal &} Smit
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

    public CBMCGrowStraightAlkane(PotentialMaster potentialMaster,
                                  IRandom random, IntegratorMC integrator, Box p, ISpecies species,
                                  Space _space, int n, int NTrials) {
        super(potentialMaster, random, _space, integrator, p, n, NTrials);

        setChainlength(n);
        sumW = 0.0;

        ISpecies[] sp = new ISpecies[1];
        sp[0] = species;

        vex = _space.makeVector();
        temp = _space.makeVector();
        a = new double[numTrial];
        b = new double[chainlength]; // used to store old rosenbluth factors
        storePos = new Vector[numTrial];
        // angleSet = new double[numTrial];
        for (int k = 0; k < numTrial; k++) {
            storePos[k] = _space.makeVector();
        }
        tempCloser = _space.makeVector();
        tempFarther = _space.makeVector();

        positionSource = new RandomPositionSourceRectangular(_space, random);
    }

    public void setBox(Box newBox) {
        super.setBox(newBox);
        positionSource.setBox(box);
    }

    public void setPositionSource(RandomPositionSource newPositionSource) {
        positionSource = newPositionSource;
        if (box != null) {
            positionSource.setBox(box);
        }
    }
    
    public RandomPositionSource getPositionSource() {
        return positionSource;
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
                    (atomList.getAtom(i).getPosition()).E(
                            positionSource.randomPosition());
                } else if (i == beginIndex + dir) { // If we're placing the
                    // second atom of a molecule
                    atomList.getAtom(i).getPosition().E(
                            calcRandomBond());
                    atomList.getAtom(i).getPosition().PE(
                             atomList.getAtom(i - dir).getPosition());
                } else if (i == beginIndex + dir * 2) {// If we're placing the
                    // third atom of a
                    // molecule
                    atomList.getAtom(i).getPosition().E(
                            calcRandomBondWithAngle(atomList.getAtom(i
                                    - dir), atomList
                                    .getAtom(i - 2 * dir)));
                    atomList.getAtom(i).getPosition().PE(
                             atomList.getAtom(i - dir).getPosition());
                } else {// For the rest of the atoms in a molecule
                    atomList.getAtom(i).getPosition().E(
                            calcRandomBondWithAngleAndTorsion(
                                    atomList.getAtom(i - dir),
                                    atomList.getAtom(i - 2 * dir),
                                    atomList.getAtom(i - 3 * dir)));
                    atomList.getAtom(i).getPosition().PE(
                             atomList.getAtom(i - dir).getPosition());
                }

                // store new position & angle
                storePos[k].E( atomList.getAtom(i).getPosition());

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
             atomList.getAtom(i).getPosition().E(storePos[pickThisOne]);

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
                    (atomList.getAtom(i).getPosition()).E(
                            positionSource.randomPosition());
                } else if (i == beginIndex + dir) { // If we're placing the
                    // second atom of a molecule
                    atomList.getAtom(i).getPosition().E(
                            calcRandomBond());
                    atomList.getAtom(i).getPosition().PE(
                             atomList.getAtom(i - dir).getPosition());
                } else if (i == beginIndex + dir * 2) {// If we're placing the
                    // third atom of a
                    // molecule
                    atomList.getAtom(i).getPosition().E(
                            calcRandomBondWithAngle(atomList.getAtom(i
                                    - dir), atomList
                                    .getAtom(i - 2 * dir)));
                    atomList.getAtom(i).getPosition().PE(
                             atomList.getAtom(i - dir).getPosition());
                } else {// For the rest of the atoms in a molecule
                    atomList.getAtom(i).getPosition().E(
                            calcRandomBondWithAngleAndTorsion(
                                    atomList.getAtom(i - dir),
                                    atomList.getAtom(i - 2 * dir),
                                    atomList.getAtom(i - 3 * dir)));
                    atomList.getAtom(i).getPosition().PE(
                             atomList.getAtom(i - dir).getPosition());
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
             atomList.getAtom(i).getPosition().E(positionOld[i]);
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
    protected Vector calcRandomBond() {
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
    protected Vector calcRandomBondWithAngle(IAtom a, IAtom b) {
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
    protected Vector calcRandomBondWithAngleAndTorsion(IAtom a, IAtom b,
                                                       IAtom c) {
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
            temp.XE(tempCloser);
            tempCloser.XE(tempFarther);
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
    protected abstract double calcBondTorsionalEnergy(Vector v);

    public abstract double energyChange();

    boolean forward; // indicates direction of growth

    protected double bondlength;

    double sumW;

    Vector temp;

    Vector vex;

    double[] a;

    double[] b;

    Vector[] storePos;

    Vector tempCloser;

    Vector tempFarther;
    
    protected RandomPositionSource positionSource;

}
