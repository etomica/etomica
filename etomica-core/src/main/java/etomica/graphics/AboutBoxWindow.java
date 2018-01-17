/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;


import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.MatteBorder;


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

	private static final int GRAPHIC_WIDTH    = 200;
	private static final int GRAPHIC_HEIGHT   = 40;
	private static final int GRAPHIC_X        = 30;
    private static final int GRAPHIC_Y        = 5;

    private static final int ABOUT_WIDTH      = 400;
    private static final int ABOUT_HEIGHT     = 340;
    private static final int ABOUT_X          = 10;
    private static final int ABOUT_Y          = 10;

    private static final int CREDIT_FONT_SIZE       = 12;

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

		this.getContentPane().setLayout(new GridBagLayout());
		this.setResizable(false);
		this.setLocation(ABOUT_X, ABOUT_Y);
		this.setModal(true);

		GridBagConstraints gbc = new GridBagConstraints();

		Font creditFont = new Font("", Font.PLAIN, CREDIT_FONT_SIZE);

    	setSize(ABOUT_WIDTH, ABOUT_HEIGHT);
    	getContentPane().setBackground(Color.WHITE);

		JLabel[] creditLabels = null;

		gbc.gridx = 0; gbc.gridy = 0;
		getContentPane().add(etomicaGraphic, gbc);
		gbc.gridy++;

	    gbc.insets = new Insets(10, 0, 10, 0);

		// Add personal acknowledgements
		if(creationCredits != null || softwareCredits != null) {
		    JPanel acknowledgementPanel = new JPanel();			
		    acknowledgementPanel.setBackground(Color.WHITE);

			// Add Created By Credits
			if (creationCredits != null) {

		    	JPanel creditPanel = new JPanel();
		    	 
		    	creditPanel.setBorder(new MatteBorder(0, 0, 0, 15, Color.WHITE));
		    	creditPanel.setBackground(Color.WHITE);
		    	creditPanel.setLayout(new GridLayout(creationCredits.length+1, 1));

		    	JLabel cb = new JLabel("Created By : ");
		    	cb.setFont(creditFont);
		    	creditPanel.add(cb);

		    	creditLabels = new JLabel[creationCredits.length];

		    	for(int i = 0; i < creationCredits.length; i++) {
		    		creditLabels[i] = new JLabel(creationCredits[i]);
		    	}

		    	for(int i = 0; i < creditLabels.length; i++) {
		    		creditLabels[i].setFont(creditFont);
		    		creditPanel.add(creditLabels[i]);
		    	}
		    	acknowledgementPanel.add(creditPanel);
		    }

		    creditLabels = null;

		    // Add Software Support Credits
		    if (softwareCredits != null) {

		    	JPanel softwareCreditPanel = new JPanel();
		    	softwareCreditPanel.setBackground(Color.WHITE);
		    	softwareCreditPanel.setBorder(new MatteBorder(0, 15, 0, 0, Color.WHITE));
		    	softwareCreditPanel.setLayout(new GridLayout(softwareCredits.length+1, 1));

		    	creditLabels = new JLabel[softwareCredits.length];

		    	for(int i = 0; i < softwareCredits.length; i++) {
		    		creditLabels[i] = new JLabel(softwareCredits[i]);
		    	}

		    	JLabel ss = new JLabel("Software Support :");
		    	ss.setFont(creditFont);
		    	softwareCreditPanel.add(ss);

		    	for(int i = 0; i < creditLabels.length; i++) {
		    		creditLabels[i].setFont(creditFont);
		    		softwareCreditPanel.add(creditLabels[i]);
		    	}
		    	acknowledgementPanel.add(softwareCreditPanel);
		    }

		    getContentPane().add(acknowledgementPanel, gbc);
			gbc.gridy++;
		}

		// Add SUNY and DE
        JPanel sunyPanel = new JPanel();
        sunyPanel.setBackground(Color.WHITE);
    	sunyPanel.setLayout(new GridLayout(2, 0));
        
        JLabel sunyLabel = new JLabel(SUNY_TEXT);
    	sunyLabel.setFont(creditFont);
    	sunyPanel.add(sunyLabel);
        JLabel ceLabel = new JLabel(CE_TEXT);
    	ceLabel.setFont(creditFont);
    	sunyPanel.add(ceLabel);
        
        getContentPane().add(sunyPanel, gbc);
		gbc.gridy++;

		// Add NSF
        JLabel nsfLabel = new JLabel(NSF_TEXT);
        nsfLabel.setFont(creditFont);
    	getContentPane().add(nsfLabel, gbc);
		gbc.gridy++;

		// Add link to homepage
		// For now, this is just text.  Not sure how
		// to add a link in any java version except for 6.
		// In java 1.6, can use the java.awt.Desktop, but that
		// is not available before java 6.
        JLabel wwwLabel = new JLabel(WWW_TEXT);
        wwwLabel.setFont(creditFont);
    	getContentPane().add(wwwLabel, gbc);
		gbc.gridy++;

	   	// Add Close button
		JButton closeBtn = new JButton("OK");
		closeBtn.addActionListener(event -> this.setVisible(false));
		closeBtn.setBackground(Color.LIGHT_GRAY);
    	getContentPane().add(closeBtn, gbc);
		gbc.gridy++;
	}

	public void setVisible(boolean b) {
        this.setLocation(owner.getLocationOnScreen().x + 10,
        		         owner.getLocationOnScreen().y + 10);
		super.setVisible(b);
	}

	private static JPanel makeEtomicaGraphic() {
		JPanel panel = new JPanel();
		panel.setSize(new Dimension(GRAPHIC_WIDTH, GRAPHIC_HEIGHT));
		panel.setLocation(GRAPHIC_X, GRAPHIC_Y);
		panel.setBackground(Color.WHITE);

		final Icon etomicaImage = new ImageIcon(AboutBoxWindow.class.getResource(ETOMICA_IMAGE));
		JLabel pic = new JLabel();
		pic.setIcon(etomicaImage);
		panel.add(pic);
		return panel;
	}
}
