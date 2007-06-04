package etomica.graphics;

import java.awt.Color;
import java.awt.TextArea;

import javax.swing.JFrame;

import etomica.action.Action;
import etomica.atom.AtomSet;
import etomica.atom.IAtomKinetic;
import etomica.phase.Phase;
import etomica.space.IVector;

/**
 * Action that opens a new window and dumps the velocities into the window.
 * @author andrew
 */
public class ActionVelocityWindow implements Action {
    private final AtomSet leafList;
    
    public ActionVelocityWindow(Phase phase) {
        leafList = phase.getSpeciesMaster().getLeafList();
    }
    
    public void actionPerformed() {
        JFrame f = new JFrame();
        TextArea textArea = new TextArea();
        textArea.setEditable(false);
        textArea.setBackground(Color.white);
        textArea.setForeground(Color.black);
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            IVector vel = a.getVelocity();
            String str = Double.toString(vel.x(0));
            for (int i=1; i<vel.getD(); i++) {
                str += " "+Double.toString(vel.x(i));
            }
            textArea.append(str+"\n");
        }
        f.add(textArea);
        f.pack();
        f.setSize(400,600);
        f.setVisible(true);
    }
}