/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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