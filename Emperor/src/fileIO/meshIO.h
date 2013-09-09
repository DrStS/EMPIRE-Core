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
/***********************************************************************************************//**
 * \file meshIO.h
 * This file holds the class MeshIO. 
 * \date 01/30/2013
 **************************************************************************************************/
#ifndef MESHIO_H_
#define MESHIO_H_

namespace EMPIRE {
/********//**
 * \brief Class MeshIO contains a set of classes and functions for general mesh import/export.
 ***********/
class MeshIO {
public:

  /********//**
  * \brief Class Node stores node data
  ***********/
  class Node {
  public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _id node ID
     * \param[in] _coordinate node coordinates (x, y, z)
     * \author Tianyang Wang
     ***********/
    Node(int _id, double *_coordinate);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~Node();
    /// node ID
    const int id;
    /// node coordinate (x, y, z)
    const double * const coordinate;
  };

  /********//**
   * \brief Class Element stores element data
   ***********/
  class Element {
  public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _id element ID
     * \param[in] _nodeIds IDs of the element's nodes
     * \author Tianyang Wang, Michael Andre
     ***********/
    Element(int _id, int *_nodeIds, int _numberOfNodes);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~Element();
    /// element ID
    const int id;
    /// IDs of the element's nodes
    const int * nodeIds;
    /// number of nodes in this element
    const int numberOfNodes;
  };

  
  /***********************************************************************************************
   * \brief Compare two elements by their number of nodes
   * \param[in] first first element
   * \param[in] second second element
   * \return true if the first element has fewer nodes. Otherwise false.
   * \author Michael Andre
   ***********/
  static bool compare_elem_by_num_node(Element *first, Element *second);

  /***********************************************************************************************
   * \brief Compare two elements by their ids
   * \param[in] first first element
   * \param[in] second second element
   * \return true if the first element has a smaller id
   * \author Michael Andre
   ***********/
  static bool compare_elem_by_id(Element *first, Element *second);

  /***********************************************************************************************
   * \brief Compare two nodes by their ids
   * \param[in] first first node
   * \param[in] second second second
   * \return true if the first node has a smaller id
   * \author Michael Andre
   ***********/
  static bool compare_node_by_id(Node *first, Node *second);

};

} /* namespace EMPIRE */

#endif /* MESHIO_H_ */
