
/**
 * Pop up a (modal) error dialog and wait for a user to press "continue".
 */

package simulate.gui;

import java.awt.*;
import java.awt.event.*;
import java.beans.*;

public class ErrorDialog extends MessageDialog {

    public ErrorDialog(Frame frame, String message) {
	super(frame, "Error", message);
    }

}
