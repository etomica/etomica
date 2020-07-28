/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;


import net.miginfocom.swing.MigLayout;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.Arrays;

import javax.swing.*;
import javax.swing.border.MatteBorder;
import javax.swing.event.HyperlinkEvent;


public class AboutBoxWindow extends JDialog {

    private static final JPanel etomicaGraphic = makeEtomicaGraphic();

    private Component owner = null;

	private String[] creationCredits = null;
	private String[] softwareCredits = null;

	private static final String ETOMICA_IMAGE = "etomica.jpg";

    private static final String SUNY_TEXT = "University at Buffalo, The State University of New York";
    private static final String CE_TEXT = "Department of Chemical & Biological Engineering";
    private static final String NSF_TEXT = "Funding support provided by the National Science Foundation.";
    private static final String WWW_TEXT = "http://www.etomica.org";

    private static final int ABOUT_X          = 10;
    private static final int ABOUT_Y          = 10;

	public AboutBoxWindow(Component owner, String appName) {

	    this.setTitle(appName);
		this.owner = owner;
		createGUI();
	}

	public AboutBoxWindow(Component owner) {

		this.owner = owner;
		createGUI();
	}

	public AboutBoxWindow(Component owner, String appName,
			              String[] created, String[] software) {

		this.setTitle(appName);
		this.owner = owner;
		creationCredits = created;
		softwareCredits = software;
		createGUI();

	}

	private void createGUI() {

		this.getContentPane().setLayout(new MigLayout("flowy, insets dialog", "[align center]"));
		this.setResizable(false);
		this.setLocation(ABOUT_X, ABOUT_Y);
		this.setModal(true);

    	getContentPane().setBackground(Color.WHITE);

		getContentPane().add(etomicaGraphic);

		// Add link to homepage
		JEditorPane wwwLabel = new JEditorPane();
		wwwLabel.setEditable(false);
		wwwLabel.setBackground(null);
		wwwLabel.setContentType("text/html");
		wwwLabel.setText("<a href=\"" + WWW_TEXT + "\">" + WWW_TEXT + "</a>");
		wwwLabel.addHyperlinkListener(e -> {
			if (e.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
				if (e.getURL() != null) {
					try {
						Desktop.getDesktop().browse(e.getURL().toURI());
					} catch (IOException | URISyntaxException ioException) {
						ioException.printStackTrace();
					}
				}
			}
		});
		getContentPane().add(wwwLabel);

		// Add personal acknowledgements
		if(creationCredits != null || softwareCredits != null) {
		    JPanel acknowledgementPanel = new JPanel(new MigLayout("gapx 30px"));
		    acknowledgementPanel.setBackground(Color.WHITE);

			// Add Created By Credits
			if (creationCredits != null) {
                JTextArea textArea = makeTextArea();
                textArea.append("Created By:");

				for (String creationCredit : creationCredits) {
					textArea.append("\n" + creationCredit);
				}
		    	acknowledgementPanel.add(textArea);
		    }

		    // Add Software Support Credits
		    if (softwareCredits != null) {
		    	JTextArea textArea = makeTextArea();
		    	textArea.append("Software Support:");

				for (String softwareCredit : softwareCredits) {
                    textArea.append("\n" + softwareCredit);
				}
		    	acknowledgementPanel.add(textArea);
		    }

		    getContentPane().add(acknowledgementPanel);
		}

		// Add SUNY and CE
        JTextArea sunyArea = makeTextArea();
		sunyArea.append(SUNY_TEXT + "\n");
		sunyArea.append(CE_TEXT);
        getContentPane().add(sunyArea);

		// Add NSF
        JTextArea nsfArea = makeTextArea();
        nsfArea.append(NSF_TEXT);
    	getContentPane().add(nsfArea);

	   	// Add Close button
		JButton closeBtn = new JButton("OK");
		closeBtn.addActionListener(event -> this.setVisible(false));
		closeBtn.setBackground(Color.LIGHT_GRAY);
    	getContentPane().add(closeBtn, "gaptop 15px");
	}

	public void setVisible(boolean b) {
	    this.pack();
	    // prevent jumping to new location on close
	    if (b) {
			this.setLocation(owner.getLocationOnScreen().x + 10,
					owner.getLocationOnScreen().y + 10);
		}
		super.setVisible(b);
	}

	private static JTextArea makeTextArea() {
		JTextArea f = new JTextArea();
		f.setEditable(false);
		f.setBorder(null);
		f.setBackground(null);
		return f;
	}

	private static JPanel makeEtomicaGraphic() {
		JPanel panel = new JPanel();
		panel.setBackground(Color.WHITE);

		final Icon etomicaImage = new ImageIcon(AboutBoxWindow.class.getResource(ETOMICA_IMAGE));
		JLabel pic = new JLabel();
		pic.setIcon(etomicaImage);
		panel.add(pic);
		return panel;
	}
}
