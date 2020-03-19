// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443211. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability see http://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

#ifndef MFEM_MESH_EXTRAS
#define MFEM_MESH_EXTRAS

#include "mfem.hpp"
#include <sstream>

namespace mfem
{

namespace common
{

class ElementMeshStream : public std::stringstream
{
public:
   ElementMeshStream(Element::Type e);
};

void
MergeMeshNodes(Mesh * mesh, int logging = 0);

void
IdentifyPeriodicMeshVertices(const Mesh & mesh,
                             const std::vector<Vector> & trans_vecs,
                             Array<int> & v2v,
                             int logging);
Mesh *
MakePeriodicMesh(Mesh * mesh, const Array<int> & v2v,
                 int logging = 0);

/// Convert a set of attribute numbers to a marker array
/** The marker array will be of size max_attr and it will contain only
    zeroes and ones.  Ones indicate which attribute numbers are
    present in the attrs array.  In the special case when attrs has a
    single entry equal to -1 the marker array will contain all ones.
*/
void AttrToMarker(int max_attr, const Array<int> &attrs, Array<int> &marker);

} // namespace common

} // namespace mfem

#endif
