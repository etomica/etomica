package etomica.graph.engine.impl;

import etomica.graph.engine.Parser.Command;

public class CommandImpl implements Command {

  private String command;

  public CommandImpl(String command) {

    this.command = command;
  }

  public String getCommand() {

    return this.command;
  }
}