/**
 * 
 */
package etomica.graphics;

import java.awt.Color;
import java.awt.TextArea;

import javax.swing.JFrame;

import etomica.action.Action;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.phase.Phase;
import etomica.space.ICoordinateKinetic;
import etomica.space.IVector;

/**
 * Action that opens a new window and dumps the velocities into the window.
 * @author andrew
 */
public class ActionVelocityWindow implements Action {
    private final AtomIteratorLeafAtoms iterator;
    
    public ActionVelocityWindow(Phase phase) {
        iterator = new AtomIteratorLeafAtoms(phase);
    }
    
    public String getLabel() {
        return "a label";
    }
    
    public void actionPerformed() {
        JFrame f = new JFrame();
        TextArea textArea = new TextArea();
        textArea.setEditable(false);
        textArea.setBackground(Color.white);
        textArea.setForeground(Color.black);
        iterator.reset();
        while (iterator.hasNext()) {
            AtomLeaf atom = (AtomLeaf)iterator.nextAtom();
            IVector vel = ((ICoordinateKinetic)atom.getCoord()).getVelocity();
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