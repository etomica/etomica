package etomica.normalmode.nptdemo.fluid;

import java.awt.Color;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.graphics.ColorSchemeCollectiveAgent;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.space.ISpace;

/**
 * Color atoms based on being neighbors of the reference atom
 *
 * @author Andrew Schultz
 */
public class ColorSchemeFluidOverlap extends ColorSchemeCollectiveAgent {
    
    public ColorSchemeFluidOverlap(ISpace space, PotentialMasterList potentialMaster, IBox box, IStuff stuff) {
        super(box);
        this.box = box;
        nOverlaps = new int[box.getLeafList().getAtomCount()];
        neighborManager = potentialMaster.getNeighborManager(box);
        pi = space.makeVector();
        pj = space.makeVector();
        dr = space.makeVector();
        this.stuff = stuff;
    }

    public void setDisplayDensity(double newDisplayDensity) {
        displayDensity = newDisplayDensity;
    }

    public double getDisplayDensity() {
        return displayDensity;
    }
    
    public synchronized void colorAllAtoms() {
        
        IAtomList leafList = box.getLeafList();
        double vOld = box.getBoundary().volume();
        int nAtoms = box.getLeafList().getAtomCount();
        double vNew = nAtoms/displayDensity;
//        System.out.println(displayDensity+" "+vOld+" "+vNew);
        
        IVector[] xScaling = stuff.stuff();
        double rScale = Math.sqrt(vNew/vOld);
        // T=1, D=2
        
        
        double sigma = 1.0;
        double scaledSig = sigma/rScale;
        double sig2 = scaledSig*scaledSig;
        
        //color all atoms according to their type
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            nOverlaps[iLeaf] = 0;
        }
        IBoundary boundary = box.getBoundary();
        for (int i=0; i<leafList.getAtomCount(); i++) {
            //color blue the neighbor atoms in same group
            IAtom atom = leafList.getAtom(i);
            pi.Ea1Tv1((vNew-vOld)/rScale, xScaling[atom.getLeafIndex()]);
            pi.PE(atom.getPosition());
//            System.out.println("dpi "+i+" "+pi);

            IAtomList list = neighborManager.getDownList(atom)[0];
            for (int j=0; j<list.getAtomCount(); j++) {
                IAtom jAtom = list.getAtom(j);

                pj.Ea1Tv1((vNew-vOld)/rScale, xScaling[jAtom.getLeafIndex()]);
//                System.out.println("dpj "+j+" "+pj);
                pj.PE(jAtom.getPosition());

                dr.Ev1Mv2(pi, pj);
                boundary.nearestImage(dr);
                double r2 = dr.squared();
                if (r2 < sig2) {
                    nOverlaps[i]++;
                    nOverlaps[jAtom.getLeafIndex()]++;
                }
            }
        }
        for (int i=0; i<leafList.getAtomCount(); i++) {
            //color green the target atom 
            agentManager.setAgent(leafList.getAtom(i), colors[nOverlaps[i]]);
        }
    }

    private final NeighborListManager neighborManager;
    protected final IVectorMutable dr;
    protected final IVectorMutable pi, pj;
    protected final int[] nOverlaps;
    protected double displayDensity;
    protected final IBox box;
    protected final IStuff stuff;
    protected Color[] colors = new Color[]{Color.RED, Color.BLUE, Color.GREEN, Color.BLACK, Color.CYAN, Color.PINK, new Color(0.5f, 0.0f, 0.5f)};
}