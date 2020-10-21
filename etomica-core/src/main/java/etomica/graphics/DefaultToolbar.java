/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;


public class DefaultToolbar {

	private JMenuBar mBar;
	private String[] DEFAULT_AUTHORS = {"Dr. David A. Kofke", "Dr. Andrew Schultz" };
    private String[] DEFAULT_CONTRIB = { "Alex Kofke", "Robert Rassler" };
    private String[] authors = DEFAULT_AUTHORS;
    private String[] contributors = DEFAULT_CONTRIB;

    public DefaultToolbar() {
    	this(null, "");
    }

    public DefaultToolbar(final JPanel parent, String appName) {
    	mBar = new JMenuBar();

    	// File menu
    	JMenu fileMenu = new JMenu("File");

    	JMenuItem exitBtn = new JMenuItem("Exit");
    	exitBtn.addActionListener(new ActionListener() {
    		public void actionPerformed(ActionEvent ev) {
    			System.exit(0);
    		}
    	});
    	fileMenu.add(exitBtn);

    	// Help menu
    	JMenu helpMenu = new JMenu("Help");
    	final String aboutString = "About " + appName;
    	JMenuItem aboutBtn = new JMenuItem(aboutString);
    	aboutBtn.addActionListener(new ActionListener() {
    		public void actionPerformed(ActionEvent ev) {
    			AboutBoxWindow about =
    				new AboutBoxWindow(parent,
    					               aboutString,
    					               authors,
    					               contributors);
    			about.setVisible(true);
    		}
    	});

    	aboutBtn.setEnabled(true);
    	helpMenu.add(aboutBtn);

    	mBar.add(fileMenu);
    	mBar.add(helpMenu);

    }

    public void addAuthor(String newAuth) {
    	String[] temp;

    	if (authors == null) {
    		temp = new String[] { newAuth };
    	}
    	else {
    		temp = new String[authors.length+1];
    		for(int i = 0; i < authors.length; i++) {
    			temp[i] = authors[i];
    		}
    		temp[authors.length] = newAuth;
    	}
    	
    	authors = temp;
    }

    public void addContributor(String newContrib) {
    	String[] temp;

    	if (contributors == null) {
    		temp = new String[] { newContrib };
    	}
    	else {
    		temp = new String[contributors.length+1];
    		for(int i = 0; i < contributors.length; i++) {
    			temp[i] = contributors[i];
    		}
    		temp[contributors.length] = newContrib;
    	}
    	
    	contributors = temp;
    }

    public JMenuBar graphic() {
    	return mBar;
    }

}
