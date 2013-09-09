/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Munich
 *
 *  All rights reserved.
 *
 *  This file is part of EMPIRE.
 *
 *  EMPIRE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EMPIRE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with EMPIRE.  If not, see http://www.gnu.org/licenses/.
 */
#include "AbstractExtrapolator.h"
#include "ConnectionIO.h"
#include <assert.h>

namespace EMPIRE {

AbstractExtrapolator::AbstractExtrapolator(std::string _name) :
        name(_name), doExtrapolate(false) {
}

AbstractExtrapolator::~AbstractExtrapolator() {
}

  void AbstractExtrapolator::setInputAndOutput(const ConnectionIO *_input, ConnectionIO *_output) {
    assert(inputVec.size()==0);
    assert(outputVec.size()==0);
    assert(_input!=NULL);
    assert(_output!=NULL);
    inputVec.push_back(_input);
    outputVec.push_back(_output);
    assert(_input->type == _output->type);
  }

  void AbstractExtrapolator::addInput(const ConnectionIO *_input) {
    assert(_input!=NULL);
    inputVec.push_back(_input);

    if (inputVec.size() > 1) {
      assert(inputVec.back()->type == inputVec[0]->type);
    }
  }

  void AbstractExtrapolator::addOutput(ConnectionIO *_output) {
    assert(_output!=NULL);
    outputVec.push_back(_output);

    if (outputVec.size() > 1) {
      assert(outputVec.back()->type == outputVec[0]->type);
    }
  }

void AbstractExtrapolator::setDoExtrapolate() {
    doExtrapolate = true;
}

} /* namespace EMPIRE */
