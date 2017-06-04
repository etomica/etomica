package etomica.normalmode.nptdemo;

import java.awt.Color;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.api.IBoundary;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.graphics.ColorSchemeCollectiveAgent;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.CoordinateDefinition;
import etomica.space.Space;

/**
 * Color atoms based on being neighbors of the reference atom
 *
 * @author Andrew Schultz
 */
public class ColorSchemeScaledOverlap extends ColorSchemeCollectiveAgent {
    
    public ColorSchemeScaledOverlap(Space space, PotentialMasterList potentialMaster, CoordinateDefinition coordinateDefinition) {
        super(coordinateDefinition.getBox());
        this.coordinateDefinition = coordinateDefinition;
        Box box = coordinateDefinition.getBox();
        nOverlaps = new int[box.getLeafList().getAtomCount()];
        neighborManager = potentialMaster.getNeighborManager(box);
        pi = space.makeVector();
        pj = space.makeVector();
        dr = space.makeVector();
    }

    public void setPressure(double newPressure) {
        pressure = newPressure;
    }
    
    public double getPressure() {
        return pressure;
    }

    public void setDisplayDensity(double newDisplayDensity) {
        displayDensity = newDisplayDensity;
    }

    public double getDisplayDensity() {
        return displayDensity;
    }
    
    public synchronized void colorAllAtoms() {
        
        Box box = coordinateDefinition.getBox();
        IAtomList leafList = box.getLeafList();
        double vOld = box.getBoundary().volume();
        int nAtoms = box.getLeafList().getAtomCount();
        double vNew = nAtoms/displayDensity;
        double rScale = Math.sqrt(vNew/vOld);
        double latticeScale = Math.exp((pressure*(vNew-vOld))/((nAtoms-1)*1*2))/rScale;
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
            pi.E(atom.getPosition());
            Vector l = coordinateDefinition.getLatticePosition(atom);
            pi.ME(l);
            pi.TE(latticeScale);
            pi.PE(l);

            IAtomList list = neighborManager.getDownList(atom)[0];
            for (int j=0; j<list.getAtomCount(); j++) {
                IAtom jAtom = list.getAtom(j);

                pj.E(jAtom.getPosition());
                Vector lj = coordinateDefinition.getLatticePosition(jAtom);
                pj.ME(lj);
                pj.TE(latticeScale);
                pj.PE(lj);

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

    private static final long serialVersionUID = 1L;
    private final NeighborListManager neighborManager;
    protected final Vector dr;
    protected final Vector pi, pj;
    protected final int[] nOverlaps;
    protected double pressure, displayDensity;
    protected final CoordinateDefinition coordinateDefinition;
    protected Color[] colors = new Color[]{Color.RED, Color.BLUE, Color.GREEN, Color.BLACK, Color.CYAN, Color.PINK, new Color(0.5f, 0.0f, 0.5f)};
}
