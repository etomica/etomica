/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine;

import java.awt.Dimension;
import java.awt.ScrollPane;

import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;

import etomica.graph.engine.impl.EnvironmentImpl;
import etomica.graph.engine.impl.ParserImpl;

public class Launcher {

  public static void main(String[] args) {

    ConsoleWindow container = new ConsoleWindow(new ConsoleUI());
    new Interpreter(container, new EnvironmentImpl(), new ParserImpl());
    SwingUtilities.invokeLater(new ConsoleWindowRunnable(container));
  }

  public static class ConsoleWindowRunnable implements Runnable {

    private ConsoleWindow window;

    public ConsoleWindowRunnable(ConsoleWindow window) {

      this.window = window;
    }

    public void run() {

      window.setVisible(true);
    }
  }

  public static class ConsoleWindow extends JFrame implements ConsoleContainer {

    private static final long serialVersionUID = -3245527931727439689L;
    private ConsoleUI console;

    public ConsoleWindow(ConsoleUI console) {

      this.console = console;
      setPreferredSize(new Dimension(800, 600));
      setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
      setTitle("The Etomica Graph Console");

      // Add contents to the window.
      JScrollPane sp = new JScrollPane(console);
      sp.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
      sp.setPreferredSize(new Dimension(800, 600));
      add(sp);

      // Display the window.
      pack();
    }

    public Console getConsole() {

      return this.console;
    }

    public void quit() {

      setVisible(false);
    }
  }
}