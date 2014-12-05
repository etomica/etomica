/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine;

import java.awt.Color;
import java.awt.Font;
import java.awt.event.*;

import javax.swing.*;
import javax.swing.text.BadLocationException;

public class ConsoleUI extends JTextArea implements Console {

  private static final long serialVersionUID = 8775341949270012136L;

  private StringBuffer buffer = null;
  private int lineStart = 0;
  private ConsoleReader reader;

  public ConsoleUI() {

    setBackground(Color.black);
    setFont(new Font("monospaced", Font.PLAIN, 12));
    setForeground(Color.green);
    setCaretColor(Color.green);
    setLineWrap(true);
    setWrapStyleWord(true);
    addKeyListener(new KeyAdapter() {

      @Override
      public void keyPressed(KeyEvent e) {

        switch (e.getKeyCode()) {
          case KeyEvent.VK_LEFT:
          case KeyEvent.VK_BACK_SPACE:
            if (getCaretPosition() < lineStart || !canMoveLeft()) {
              e.consume();
            }
            break;
          case KeyEvent.VK_RIGHT:
          case KeyEvent.VK_DELETE:
            if (getCaretPosition() < lineStart || !canMoveRight()) {
              e.consume();
            }
            break;
          case KeyEvent.VK_ENTER:
            if (getCaretPosition() >= lineStart) {
              e.consume();
              processCommand();
            }
            break;
          case KeyEvent.VK_DOWN:
          case KeyEvent.VK_UP:
            if (getCaretPosition() >= lineStart) {
              e.consume();
            }
            break;
          case KeyEvent.VK_END:
          case KeyEvent.VK_HOME:
            break;
          default:
            e.consume();
        }
      }

      @Override
      public void keyTyped(KeyEvent e) {

        if (getCaretPosition() < lineStart) {
          e.consume();
        }
      }
    });
    printPrompt();
  }

  protected boolean canMoveLeft() {

    return getCaretPosition() > lineStart;
  }

  protected boolean canMoveRight() {

    return getCaretPosition() < getDocument().getLength();
  }

  private String extractCommand() {

    try {
      return getDocument().getText(lineStart, getDocument().getLength() - lineStart);
    }
    catch (BadLocationException e) {
      write(e);
    }
    return "";
  }

  protected String getLastLine() {

    return extractCommand();
  }

  public String getPrompt() {

    return ">";
  }

  private void printPrompt() {

    if (getDocument().getLength() != 0) {
      writeLn();
    }
    write(getPrompt());
    lineStart = getDocument().getLength();
    setCaretPosition(lineStart);
  }

  private void processCommand() {

    String command = extractCommand();
    try {
      if (reader != null) {
        reader.read(command);
      }
    }
    catch (Exception e) {
      write(e);
    }
    printPrompt();
  }

  public void setReader(ConsoleReader reader) {

    this.reader = reader;
  }

  public void updateBegin() {

    buffer = new StringBuffer("");
  }

  public void updateDone() {

    if (buffer != null) {
      append(buffer.toString());
      buffer = null;
      // updateUI, repaint, invalidate
      updateUI();
    }
  }

  private void localAppend(String value) {

    if (buffer != null) {
      buffer.append(value);
    }
    else {
      append(value);
    }
  }

  public void write(Exception e) {

    localAppend(e.getMessage());
  }

  public void write(final String value) {

    localAppend(value);
  }

  public void writeLn() {

    localAppend("\n");
  }

  public void clear() {

    setText(null);
  }
}