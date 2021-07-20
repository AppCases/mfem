// Copyright (c) 2010-2020, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

// Implementation of class LinearForm

#include "fem.hpp"

namespace mfem
{

LinearForm::LinearForm(FiniteElementSpace *f, LinearForm *lf)
   : Vector(f->GetVSize())
{
   // Linear forms are stored on the device
   UseDevice(true);

   fes = f;
   extern_lfs = 1;

   // Copy the pointers to the integrators
   dlfi = lf->dlfi;

   dlfi_delta = lf->dlfi_delta;

   blfi = lf->blfi;

   flfi = lf->flfi;
   flfi_marker = lf->flfi_marker;
}

void LinearForm::AddDomainIntegrator(LinearFormIntegrator *lfi)
{
   DeltaLFIntegrator *maybe_delta =
      dynamic_cast<DeltaLFIntegrator *>(lfi);
   if (!maybe_delta || !maybe_delta->IsDelta())
   {
      dlfi.Append(lfi);
   }
   else
   {
      dlfi_delta.Append(maybe_delta);
   }
}

void LinearForm::AddBoundaryIntegrator (LinearFormIntegrator * lfi)
{
   blfi.Append (lfi);
   blfi_marker.Append(NULL); // NULL -> all attributes are active
}

void LinearForm::AddBoundaryIntegrator (LinearFormIntegrator * lfi,
                                        Array<int> &bdr_attr_marker)
{
   blfi.Append (lfi);
   blfi_marker.Append(&bdr_attr_marker);
}

void LinearForm::AddBdrFaceIntegrator (LinearFormIntegrator * lfi)
{
   flfi.Append(lfi);
   flfi_marker.Append(NULL); // NULL -> all attributes are active
}
  
void LinearForm::AddBdrFaceIntegrator(LinearFormIntegrator *lfi,
                                      Array<int> &bdr_attr_marker)
{
   flfi.Append(lfi);
   flfi_marker.Append(&bdr_attr_marker);
}

/* UW - GSJ */
void LinearForm::AddFreeSurfaceIntegrator(LinearFormIntegrator *lfi,
                                          Array<int> &bdr_attr_marker)
{
    fsi.Append(lfi);
    fsi_marker.Append(&bdr_attr_marker);
}

void LinearForm::AddFreeSurfaceTraceIntegrator(LinearFormIntegrator *lfi,
                                          Array<int> &bdr_attr_marker)
{
    fsti.Append(lfi);
    fsti_marker.Append(&bdr_attr_marker);
}

void LinearForm::AddTimeLevelQIntegrator(LinearFormIntegrator *lfi,
										Array<int> &bdr_attr_marker)
{
	tlqi.Append(lfi);
	tlqi_marker.Append(&bdr_attr_marker);
}
/* UW - GSJ ends */

/* UW */
void LinearForm::AddTimeLevelIntegrator (LinearFormIntegrator * lfi)
{
   tllfi.Append (lfi);
}

/* UW */
void LinearForm::AddTimeLevelRotationIntegrator (LinearFormIntegrator * lfi)
{
   tlrotlfi.Append (lfi);
}

/* UW */
void LinearForm::AddSktBoundaryNeumannIntegrator(LinearFormIntegrator * lfi)
{
   bdrsklneufi.Append (lfi);
   bdrsklneufi_marker.Append(NULL); // NULL -> all attributes are active
}

/* UW */
void LinearForm::AddSktBoundaryNeumannIntegrator(LinearFormIntegrator * lfi,
                                                 Array<int> &bdr_attr_marker)
{
   bdrsklneufi.Append (lfi);
   bdrsklneufi_marker.Append(&bdr_attr_marker); // NULL -> all attributes are active
}

/* UW */
void LinearForm::AddSktBoundaryNeumannIntegratorWithMesh (LinearFormIntegrator * lfi)
{
   bdrsklneu_mesh_fi.Append (lfi);
   bdrsklneu_mesh_fi_marker.Append(NULL); // NULL -> all attributes are active
}

/* UW */
void LinearForm::AddSktBoundaryNeumannIntegratorWithMesh (LinearFormIntegrator * lfi,
                                              Array<int> &bdr_attr_marker)
{
   bdrsklneu_mesh_fi.Append (lfi);
   bdrsklneu_mesh_fi_marker.Append(&bdr_attr_marker); // NULL -> all attributes are active
}
  
void LinearForm::Assemble()
{
   Array<int> vdofs;
   ElementTransformation *eltrans;
   Vector elemvect;

   int i;

   Vector::operator=(0.0);

   // The above operation is executed on device because of UseDevice().
   // The first use of AddElementVector() below will move it back to host
   // because both 'vdofs' and 'elemvect' are on host.

   if (dlfi.Size())
   {
      for (i = 0; i < fes -> GetNE(); i++)
      {
         fes -> GetElementVDofs (i, vdofs);
         eltrans = fes -> GetElementTransformation (i);
         for (int k=0; k < dlfi.Size(); k++)
         {
            dlfi[k]->AssembleRHSElementVect(*fes->GetFE(i), *eltrans, elemvect);
            AddElementVector (vdofs, elemvect);
         }
      }
   }
   AssembleDelta();

   if (blfi.Size())
   {
      Mesh *mesh = fes->GetMesh();

      // Which boundary attributes need to be processed?
      Array<int> bdr_attr_marker(mesh->bdr_attributes.Size() ?
                                 mesh->bdr_attributes.Max() : 0);
      bdr_attr_marker = 0;
      for (int k = 0; k < blfi.Size(); k++)
      {
         if (blfi_marker[k] == NULL)
         {
            bdr_attr_marker = 1;
            break;
         }
         Array<int> &bdr_marker = *blfi_marker[k];
         MFEM_ASSERT(bdr_marker.Size() == bdr_attr_marker.Size(),
                     "invalid boundary marker for boundary integrator #"
                     << k << ", counting from zero");
         for (int i = 0; i < bdr_attr_marker.Size(); i++)
         {
            bdr_attr_marker[i] |= bdr_marker[i];
         }
      }

      for (i = 0; i < fes -> GetNBE(); i++)
      {
         const int bdr_attr = mesh->GetBdrAttribute(i);
         if (bdr_attr_marker[bdr_attr-1] == 0) { continue; }
         fes -> GetBdrElementVDofs (i, vdofs);
         eltrans = fes -> GetBdrElementTransformation (i);
         for (int k=0; k < blfi.Size(); k++)
         {
            if (blfi_marker[k] &&
                (*blfi_marker[k])[bdr_attr-1] == 0) { continue; }

            blfi[k]->AssembleRHSElementVect(*fes->GetBE(i), *eltrans, elemvect);

            AddElementVector (vdofs, elemvect);
         }
      }
   }
   if (flfi.Size())
   {
      FaceElementTransformations *tr;
      Mesh *mesh = fes->GetMesh();

      // Which boundary attributes need to be processed?
      Array<int> bdr_attr_marker(mesh->bdr_attributes.Size() ?
                                 mesh->bdr_attributes.Max() : 0);
      bdr_attr_marker = 0;
      for (int k = 0; k < flfi.Size(); k++)
      {
         if (flfi_marker[k] == NULL)
         {
            bdr_attr_marker = 1;
            break;
         }
         Array<int> &bdr_marker = *flfi_marker[k];
         MFEM_ASSERT(bdr_marker.Size() == bdr_attr_marker.Size(),
                     "invalid boundary marker for boundary face integrator #"
                     << k << ", counting from zero");
         for (int i = 0; i < bdr_attr_marker.Size(); i++)
         {
            bdr_attr_marker[i] |= bdr_marker[i];
         }
      }

      for (i = 0; i < mesh->GetNBE(); i++)
      {
         const int bdr_attr = mesh->GetBdrAttribute(i);
         if (bdr_attr_marker[bdr_attr-1] == 0) { continue; }

         tr = mesh->GetBdrFaceTransformations(i);
         if (tr != NULL)
         {
            fes -> GetElementVDofs (tr -> Elem1No, vdofs);
            for (int k = 0; k < flfi.Size(); k++)
            {
               if (flfi_marker[k] &&
                   (*flfi_marker[k])[bdr_attr-1] == 0) { continue; }

               flfi[k] -> AssembleRHSElementVect (*fes->GetFE(tr -> Elem1No),
                                                  *tr, elemvect);
               AddElementVector (vdofs, elemvect);
            }
         }
      }
   }

   /* UW - GSJ */
   if (fsi.Size())
   {      
      Mesh *mesh = fes->GetMesh(); // This is meshM
      FaceElementTransformations *tr;
      
      // Which boundary attributes need to be processed?
      Array<int> bdr_attr_marker(mesh->bdr_attributes.Size() ?
                                 mesh->bdr_attributes.Max() : 0);
      bdr_attr_marker = 0;
      for (int k = 0; k < fsi.Size(); k++)
      {
         if (fsi_marker[k] == NULL)
         {
            bdr_attr_marker = 1;
            break;
         }
         Array<int> &bdr_marker = *fsi_marker[k];
         MFEM_ASSERT(bdr_marker.Size() == bdr_attr_marker.Size(),
                     "invalid boundary marker for boundary face integrator #"
                     << k << ", counting from zero");
         for (int i = 0; i < bdr_attr_marker.Size(); i++)
         {
            bdr_attr_marker[i] |= bdr_marker[i];
         }
      }

      for (i = 0; i < mesh->GetNBE(); i++)
      {
         const int bdr_attr = mesh->GetBdrAttribute(i);

         if (bdr_attr_marker[bdr_attr-1] == 0) { continue; }

         tr = mesh -> GetBdrFaceTransformations (i);



		 fes -> GetElementVDofs (tr -> Elem1No, vdofs);

         for (int k = 0; k < fsi.Size(); k++)
         {
            if (fsi_marker[k] &&
                (*fsi_marker[k])[bdr_attr-1] == 0) { continue; }

            fsi[k] -> AssembleRHSElementVectWithMesh(*fes->GetFE(tr -> Elem1No), *mesh,
                                                                  *tr, i, elemvect);

            AddElementVector (vdofs, elemvect);
         }
      }
   }
   
   if (fsti.Size())
   {
	   Mesh *mesh = fes->GetMesh(); // This is meshW
	   FaceElementTransformations *tr;

	   // Which boundary attributes need to be processed?
	   Array<int> bdr_attr_marker(mesh->bdr_attributes.Size() ?
								mesh->bdr_attributes.Max() : 0);
	   bdr_attr_marker = 0;
	   for (int k = 0; k < fsti.Size(); k++)
	   {
		   if (fsti_marker[k] == NULL)
		   {
			   bdr_attr_marker = 1;
			   break;
		   }
		   Array<int> &bdr_marker = *fsti_marker[k];
		   MFEM_ASSERT(bdr_marker.Size() == bdr_attr_marker.Size(),
				   "invalid boundary marker for boundary face integrator #"
						<< k << ", counting from zero");
		   for (int i = 0; i < bdr_attr_marker.Size(); i++)
		   {
			   bdr_attr_marker[i] |= bdr_marker[i];
		   }
	   }

	   for (i = 0; i < mesh->GetNBE(); i++)
	   {
		   const int bdr_attr = mesh->GetBdrAttribute(i);

		   if (bdr_attr_marker[bdr_attr-1] == 0) { continue; }

		   tr = mesh -> GetBdrFaceTransformations (i);

		   int face;
		   mesh->GetBdrFaceToEdge(i,&face);
		   fes -> GetFaceVDofs (face, vdofs);

		   for (int k = 0; k < fsti.Size(); k++)
		   {
			   if (fsti_marker[k] &&
				  (*fsti_marker[k])[bdr_attr-1] == 0) { continue; }

			   fsti[k] -> AssembleRHSElementVectWithMesh(*fes->GetFaceElement(face), *mesh,
					   	   	   	   	   	   	   	   	   	 *tr, i, elemvect);

			   AddElementVector (vdofs, elemvect);
		   }
	   }
   }

   if (tlqi.Size())
   {
   	   Mesh *mesh = fes->GetMesh(); // This is meshW
   	   FaceElementTransformations *tr;

   	   // Which boundary attributes need to be processed?
   	   Array<int> bdr_attr_marker(mesh->bdr_attributes.Size() ?
   								mesh->bdr_attributes.Max() : 0);
   	   bdr_attr_marker = 0;
   	   for (int k = 0; k < tlqi.Size(); k++)
   	   {
   		   if (tlqi_marker[k] == NULL)
   		   {
   			   bdr_attr_marker = 1;
   			   break;
   		   }
   		   Array<int> &bdr_marker = *tlqi_marker[k];
   		   MFEM_ASSERT(bdr_marker.Size() == bdr_attr_marker.Size(),
   				   "invalid boundary marker for boundary face integrator #"
   						<< k << ", counting from zero");
   		   for (int i = 0; i < bdr_attr_marker.Size(); i++)
   		   {
   			   bdr_attr_marker[i] |= bdr_marker[i];
   		   }
   	   }

   	   for (i = 0; i < mesh->GetNBE(); i++)
   	   {
   		   const int bdr_attr = mesh->GetBdrAttribute(i);

   		   if (bdr_attr_marker[bdr_attr-1] == 0) { continue; }

   		   tr = mesh -> GetBdrFaceTransformations (i);
   		   if (tr != NULL)
		   {
   			   int face;
			   mesh->GetBdrFaceToEdge(i,&face);
			   fes -> GetFaceVDofs (face, vdofs);
			   fes -> GetElementVDofs (tr -> Elem1No, vdofs);
   			   for (int k = 0; k < tlqi.Size(); k++)
			   {
   				   if (tlqi_marker[k] &&
   						   (*tlqi_marker[k])[bdr_attr-1] == 0) { continue; }

				   tlqi[k] -> AssembleRHSElementVectWithMesh(*fes->GetFE(tr -> Elem1No), *mesh,
															 *tr, tr -> Elem1No, elemvect);

				   AddElementVector (vdofs, elemvect);
			   }
		   }


   	   }
   }
   /* UW - GSJ ends */

   /* UW */
   if (tllfi.Size())
   {
      FaceElementTransformations *tr;
      Mesh *mesh = fes -> GetMesh();
      for (i = 0; i < mesh -> GetNBE(); i++)
      {
         tr = mesh -> GetBdrFaceTransformations (i);
         if (tr != NULL)
         {
            fes -> GetElementVDofs (tr -> Elem1No, vdofs);
            for (int k = 0; k < tllfi.Size(); k++)
            {
               tllfi[k] -> AssembleRHSElementVectWithElementIndex (*fes->GetFE(tr -> Elem1No),
                                                                  *tr, i, elemvect);
               AddElementVector (vdofs, elemvect);
            }
         }
      }
   }
   
   /* UW */
   if (tlrotlfi.Size())
   {
      FaceElementTransformations *tr;
      Mesh *mesh = fes -> GetMesh();
      for (i = 0; i < mesh -> GetNBE(); i++)
      {
         tr = mesh -> GetBdrFaceTransformations (i);
         if (tr != NULL)
         {
            fes -> GetElementVDofs (tr -> Elem1No, vdofs);
            for (int k = 0; k < tlrotlfi.Size(); k++)
            {
               tlrotlfi[k] -> AssembleRHSElementVectWithMesh (*fes->GetFE(tr -> Elem1No), *mesh,
                                                                  *tr, i, elemvect);
               AddElementVector (vdofs, elemvect);
            }
         }
      }
   }
   
   /* UW */
   if (bdrsklneufi.Size())
   {
      FaceElementTransformations *ftr;
      const FiniteElement *face_fe;
      int nbdrfaces = fes->GetNBE();
      Mesh *mesh = fes -> GetMesh();
      
      for (i = 0; i < nbdrfaces; i++)
      {
         int face;
         mesh->GetBdrFaceToEdge(i, &face);
         ftr = mesh->GetFaceElementTransformations(face); // the transformation of the face
         fes->GetFaceVDofs(face, vdofs);   // the degrees of freedom related to the face
         face_fe = fes->GetFaceElement(face);   // point face_fe to the FiniteElement over the edge
            
         if (ftr != NULL)
         {
            for (int k = 0; k < bdrsklneufi.Size(); k++) // Loop over the related interals
            {
               int compute = 0;
               if (bdrsklneufi_marker[k] == NULL)
                  compute = 1;
               else
               {
                  Array<int> &bdr_marker = *bdrsklneufi_marker[k];
                  const int bdr_attr = mesh->GetBdrAttribute(i);
                  if (bdr_marker[bdr_attr-1] == 1)
                     compute = 1;
               }
               if (compute)
               {
                  bdrsklneufi[k] -> AssembleRHSElementVect (*face_fe, *ftr, elemvect);
                  AddElementVector (vdofs, elemvect);
               }
            }
         }
	   
      }
   }


   /* UW */
   if (bdrsklneu_mesh_fi.Size())
   {
      FaceElementTransformations *ftr;
      const FiniteElement *face_fe;
      int nbdrfaces = fes->GetNBE();
      Mesh *mesh = fes -> GetMesh();
      
      for (i = 0; i < nbdrfaces; i++)
      {
         int face;
         mesh->GetBdrFaceToEdge(i, &face);
         ftr = mesh->GetBdrFaceTransformations(i); // the transformation of the face
         fes->GetFaceVDofs(face, vdofs);   // the degrees of freedom related to the face
         face_fe = fes->GetFaceElement(face);   // point face_fe to the FiniteElement over the edge
            
         if (ftr != NULL)
         {
            for (int k = 0; k < bdrsklneu_mesh_fi.Size(); k++) // Loop over the related interals
            {
               int compute = 0;
               if (bdrsklneu_mesh_fi_marker[k] == NULL)
                  compute = 1;
               else
               {
                  Array<int> &bdr_marker = *bdrsklneu_mesh_fi_marker[k];
                  const int bdr_attr = mesh->GetBdrAttribute(i);
                  if (bdr_marker[bdr_attr-1] == 1)
                     compute = 1;
               }
               if (compute)
               {
                  bdrsklneu_mesh_fi[k] -> AssembleRHSElementVectWithMesh (*face_fe, *mesh, 
                                                                          *ftr, i, elemvect);
                  AddElementVector (vdofs, elemvect);
               }
            }
         }
	   
      }
   }   
}

void LinearForm::Update(FiniteElementSpace *f, Vector &v, int v_offset)
{
   fes = f;
   NewDataAndSize((double *)v + v_offset, fes->GetVSize());
   ResetDeltaLocations();
}

void LinearForm::AssembleDelta()
{
   if (dlfi_delta.Size() == 0) { return; }

   if (!HaveDeltaLocations())
   {
      int sdim = fes->GetMesh()->SpaceDimension();
      Vector center;
      DenseMatrix centers(sdim, dlfi_delta.Size());
      for (int i = 0; i < centers.Width(); i++)
      {
         centers.GetColumnReference(i, center);
         dlfi_delta[i]->GetDeltaCenter(center);
         MFEM_VERIFY(center.Size() == sdim,
                     "Point dim " << center.Size() <<
                     " does not match space dim " << sdim);
      }
      fes->GetMesh()->FindPoints(centers, dlfi_delta_elem_id, dlfi_delta_ip);
   }

   Array<int> vdofs;
   Vector elemvect;
   for (int i = 0; i < dlfi_delta.Size(); i++)
   {
      int elem_id = dlfi_delta_elem_id[i];
      // The delta center may be outside of this sub-domain, or
      // (Par)Mesh::FindPoints() failed to find this point:
      if (elem_id < 0) { continue; }

      const IntegrationPoint &ip = dlfi_delta_ip[i];
      ElementTransformation &Trans = *fes->GetElementTransformation(elem_id);
      Trans.SetIntPoint(&ip);

      fes->GetElementVDofs(elem_id, vdofs);
      dlfi_delta[i]->AssembleDeltaElementVect(*fes->GetFE(elem_id), Trans,
                                              elemvect);
      AddElementVector(vdofs, elemvect);
   }
}

LinearForm & LinearForm::operator=(double value)
{
   Vector::operator=(value);
   return *this;
}

LinearForm & LinearForm::operator=(const Vector &v)
{
   MFEM_ASSERT(fes && v.Size() == fes->GetVSize(), "");
   Vector::operator=(v);
   return *this;
}

LinearForm::~LinearForm()
{
   if (!extern_lfs)
   {
      int k;
      for (k=0; k < dlfi_delta.Size(); k++) { delete dlfi_delta[k]; }
      for (k=0; k < dlfi.Size(); k++) { delete dlfi[k]; }
      for (k=0; k < blfi.Size(); k++) { delete blfi[k]; }
      for (k=0; k < flfi.Size(); k++) { delete flfi[k]; }
      /*UW*/
      for (k=0; k < tllfi.Size(); k++) { delete tllfi[k]; }
      for (k=0; k < tlrotlfi.Size(); k++) { delete tlrotlfi[k]; }
      for (k=0; k < bdrsklneufi.Size(); k++) { delete bdrsklneufi[k]; }
      for (k=0; k < bdrsklneu_mesh_fi.Size(); k++) { delete bdrsklneu_mesh_fi[k]; }
      /*UW - GSJ*/
      for (k=0; k < fsi.Size(); k++) { delete fsi[k]; }
      for (k=0; k < fsti.Size(); k++) { delete fsti[k]; }
      for (k=0; k < tlqi.Size(); k++) { delete tlqi[k]; }
   }
}

}
