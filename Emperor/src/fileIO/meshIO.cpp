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
#include "meshIO.h"

namespace EMPIRE {

  MeshIO::Node::Node(int _id, double *_coordinate) :
  id(_id), coordinate(_coordinate) {
  }

  MeshIO::Node::~Node() {
    delete[] coordinate;
  }

  MeshIO::Element::Element(int _id, int *_nodeIds, int _numberOfNodes) :
    id(_id), nodeIds(_nodeIds), numberOfNodes(_numberOfNodes) {
  }

  MeshIO::Element::~Element() {
    delete[] nodeIds;
  }

  bool MeshIO::compare_elem_by_num_node(Element *first, Element *second) {
    return (first->numberOfNodes < second->numberOfNodes) ? true : false;
  }

  bool MeshIO::compare_elem_by_id(Element *first, Element *second) {
    return (first->id < second->id) ? true : false;
  }

  bool MeshIO::compare_node_by_id(Node *first, Node *second) {
    return (first->id < second->id) ? true : false;
  }

} /* namespace EMPIRE */
