// Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "mfem.hpp"
#include "unit_tests.hpp"
#include "linalg/dtensor.hpp"

using namespace mfem;

TEST_CASE("DenseMatrix init-list construction", "[DenseMatrix]")
{
   double ContigData[6] = {6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
   DenseMatrix Contiguous(ContigData, 2, 3);

   DenseMatrix Nested(
   {
      {6.0, 4.0, 2.0},
      {5.0, 3.0, 1.0}
   });

   for (int i = 0; i < Contiguous.Height(); i++)
   {
      for (int j = 0; j < Contiguous.Width(); j++)
      {
         REQUIRE(Nested(i,j) == Contiguous(i,j));
      }
   }
}

TEST_CASE("DenseMatrix LinearSolve methods",
          "[DenseMatrix]")
{
   SECTION("singular_system")
   {
      constexpr int N = 3;

      DenseMatrix A(N);
      A.SetRow(0, 0.0);
      A.SetRow(1, 0.0);
      A.SetRow(2, 0.0);

      double X[3];

      REQUIRE_FALSE(LinearSolve(A,X));
   }

   SECTION("1x1_system")
   {
      constexpr int N = 1;
      DenseMatrix A(N);
      A(0,0) = 2;

      double X[1] = { 12 };

      REQUIRE(LinearSolve(A,X));
      REQUIRE(X[0] == MFEM_Approx(6));
   }

   SECTION("2x2_system")
   {
      constexpr int N = 2;

      DenseMatrix A(N);
      A(0,0) = 2.0; A(0,1) = 1.0;
      A(1,0) = 3.0; A(1,1) = 4.0;

      double X[2] = { 1, 14 };

      REQUIRE(LinearSolve(A,X));
      REQUIRE(X[0] == MFEM_Approx(-2));
      REQUIRE(X[1] == MFEM_Approx(5));
   }

   SECTION("3x3_system")
   {
      constexpr int N = 3;

      DenseMatrix A(N);
      A(0,0) = 4; A(0,1) =  5; A(0,2) = -2;
      A(1,0) = 7; A(1,1) = -1; A(1,2) =  2;
      A(2,0) = 3; A(2,1) =  1; A(2,2) =  4;

      double X[3] = { -14, 42, 28 };

      REQUIRE(LinearSolve(A,X));
      REQUIRE(X[0] == MFEM_Approx(4));
      REQUIRE(X[1] == MFEM_Approx(-4));
      REQUIRE(X[2] == MFEM_Approx(5));
   }

}

TEST_CASE("DenseMatrix A*B^T methods",
          "[DenseMatrix]")
{
   double tol = 1e-12;

   double AtData[6] = {6.0, 5.0,
                       4.0, 3.0,
                       2.0, 1.0
                      };
   double BtData[12] = {1.0, 3.0, 5.0, 7.0,
                        2.0, 4.0, 6.0, 8.0,
                        1.0, 2.0, 3.0, 5.0
                       };

   DenseMatrix A(AtData, 2, 3);
   DenseMatrix B(BtData, 4, 3);
   DenseMatrix C(2,4);

   SECTION("MultABt")
   {
      double BData[12] = {1.0, 2.0, 1.0,
                          3.0, 4.0, 2.0,
                          5.0, 6.0, 3.0,
                          7.0, 8.0, 5.0
                         };
      DenseMatrix Bt(BData, 3, 4);

      double CtData[8] = {16.0, 12.0,
                          38.0, 29.0,
                          60.0, 46.0,
                          84.0, 64.0
                         };
      DenseMatrix Cexact(CtData, 2, 4);

      MultABt(A, B, C);
      C.Add(-1.0, Cexact);

      REQUIRE(C.MaxMaxNorm() < tol);

      Mult(A, Bt, Cexact);
      MultABt(A, B, C);
      C.Add(-1.0, Cexact);

      REQUIRE(C.MaxMaxNorm() < tol);
   }
   SECTION("MultADBt")
   {
      double DData[3] = {11.0, 7.0, 5.0};
      Vector D(DData, 3);

      double CtData[8] = {132.0, 102.0,
                          330.0, 259.0,
                          528.0, 416.0,
                          736.0, 578.0
                         };
      DenseMatrix Cexact(CtData, 2, 4);

      MultADBt(A, D, B, C);
      C.Add(-1.0, Cexact);

      REQUIRE(C.MaxMaxNorm() < tol);
   }
   SECTION("AddMultABt")
   {
      double CtData[8] = {17.0, 17.0,
                          40.0, 35.0,
                          63.0, 53.0,
                          88.0, 72.0
                         };
      DenseMatrix Cexact(CtData, 2, 4);

      C(0, 0) = 1.0; C(0, 1) = 2.0; C(0, 2) = 3.0; C(0, 3) = 4.0;
      C(1, 0) = 5.0; C(1, 1) = 6.0; C(1, 2) = 7.0; C(1, 3) = 8.0;

      AddMultABt(A, B, C);
      C.Add(-1.0, Cexact);

      REQUIRE(C.MaxMaxNorm() < tol);

      MultABt(A, B, C);
      C *= -1.0;
      AddMultABt(A, B, C);
      REQUIRE(C.MaxMaxNorm() < tol);
   }
   SECTION("AddMultADBt")
   {
      double DData[3] = {11.0, 7.0, 5.0};
      Vector D(DData, 3);

      double CtData[8] = {133.0, 107.0,
                          332.0, 265.0,
                          531.0, 423.0,
                          740.0, 586.0
                         };
      DenseMatrix Cexact(CtData, 2, 4);

      C(0, 0) = 1.0; C(0, 1) = 2.0; C(0, 2) = 3.0; C(0, 3) = 4.0;
      C(1, 0) = 5.0; C(1, 1) = 6.0; C(1, 2) = 7.0; C(1, 3) = 8.0;

      AddMultADBt(A, D, B, C);
      C.Add(-1.0, Cexact);

      REQUIRE(C.MaxMaxNorm() < tol);

      MultADBt(A, D, B, C);
      C *= -1.0;
      AddMultADBt(A, D, B, C);
      REQUIRE(C.MaxMaxNorm() < tol);

      DData[0] = 1.0; DData[1] = 1.0; DData[2] = 1.0;
      MultABt(A, B, C);
      C *= -1.0;
      AddMultADBt(A, D, B, C);
      REQUIRE(C.MaxMaxNorm() < tol);
   }
   SECTION("AddMult_a_ABt")
   {
      double a = 3.0;

      double CtData[8] = { 49.0,  41.0,
                           116.0,  93.0,
                           183.0, 145.0,
                           256.0, 200.0
                         };
      DenseMatrix Cexact(CtData, 2, 4);

      C(0, 0) = 1.0; C(0, 1) = 2.0; C(0, 2) = 3.0; C(0, 3) = 4.0;
      C(1, 0) = 5.0; C(1, 1) = 6.0; C(1, 2) = 7.0; C(1, 3) = 8.0;

      AddMult_a_ABt(a, A, B, C);
      C.Add(-1.0, Cexact);

      REQUIRE(C.MaxMaxNorm() < tol);

      MultABt(A, B, C);
      AddMult_a_ABt(-1.0, A, B, C);

      REQUIRE(C.MaxMaxNorm() < tol);
   }
}


TEST_CASE("LUFactors RightSolve", "[DenseMatrix]")
{
   double tol = 1e-12;

   // Zero on diagonal forces non-trivial pivot
   double AData[9] = { 0.0, 0.0, 3.0, 2.0, 2.0, 2.0, 2.0, 0.0, 4.0 };
   double BData[6] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
   int ipiv[3];

   DenseMatrix A(AData, 3, 3);
   DenseMatrix B(BData, 2, 3);

   DenseMatrixInverse Af1(A);
   DenseMatrix Ainv;
   Af1.GetInverseMatrix(Ainv);

   LUFactors Af2(AData, ipiv);
   Af2.Factor(3);

   DenseMatrix C(2,3);
   Mult(B, Ainv, C);
   Af2.RightSolve(3, 2, B.GetData());
   C -= B;

   REQUIRE(C.MaxMaxNorm() < tol);
}

TEST_CASE("DenseTensor LinearSolve methods",
          "[DenseMatrix]")
{

   int N = 3;
   DenseMatrix A(N);
   A(0,0) = 4; A(0,1) =  5; A(0,2) = -2;
   A(1,0) = 7; A(1,1) = -1; A(1,2) =  2;
   A(2,0) = 3; A(2,1) =  1; A(2,2) =  4;

   double X[3] = { -14, 42, 28 };

   int NE = 10;
   Vector X_batch(N*NE);
   DenseTensor A_batch(N,N,NE);

   auto a_batch = mfem::Reshape(A_batch.HostWrite(),N,N,NE);
   auto x_batch = mfem::Reshape(X_batch.HostWrite(),N,NE);
   // Column major
   for (int e=0; e<NE; ++e)
   {

      for (int r=0; r<N; ++r)
      {
         for (int c=0; c<N; ++c)
         {
            a_batch(c, r, e) = A.GetData()[c+r*N];
         }
         x_batch(r,e) = X[r];
      }
   }

   Array<int> P;
   BatchLUFactor(A_batch, P);
   BatchLUSolve(A_batch, P, X_batch);

   auto xans_batch = mfem::Reshape(X_batch.HostRead(),N,NE);
   REQUIRE(LinearSolve(A,X));
   for (int e=0; e<NE; ++e)
   {
      for (int r=0; r<N; ++r)
      {
         REQUIRE(xans_batch(r,e) == MFEM_Approx(X[r]));
      }
   }
}

TEST_CASE("DenseTensor copy", "[DenseMatrix][DenseTensor]")
{
   DenseTensor t1(2,3,4);
   for (int i=0; i<t1.TotalSize(); ++i)
   {
      t1.Data()[i] = i;
   }
   DenseTensor t2(t1);
   DenseTensor t3;
   t3 = t1;
   REQUIRE(t2.SizeI() == t1.SizeI());
   REQUIRE(t2.SizeJ() == t1.SizeJ());
   REQUIRE(t2.SizeK() == t1.SizeK());

   REQUIRE(t3.SizeI() == t1.SizeI());
   REQUIRE(t3.SizeJ() == t1.SizeJ());
   REQUIRE(t3.SizeK() == t1.SizeK());

   REQUIRE(t2.Data() != t1.Data());
   REQUIRE(t3.Data() != t1.Data());

   for (int i=0; i<t1.TotalSize(); ++i)
   {
      REQUIRE(t2.Data()[i] == t1.Data()[i]);
      REQUIRE(t3.Data()[i] == t1.Data()[i]);
   }
}

#ifdef MFEM_USE_LAPACK

enum class TestCase { GenEigSPD, GenEigGE, SVD};
std::string TestCaseName(TestCase testcase)
{
   switch (testcase)
   {
      case TestCase::GenEigSPD:
         return "Generalized Eigenvalue problem for an SPD matrix";
      case TestCase::GenEigGE:
         return "Generalized Eigenvalue problem for a general matrix";
      case TestCase::SVD:
         return "Singular Value Decomposition for a general matrix";
   }
   return "";
}

TEST_CASE("Eigensystem Problems",
          "[DenseMatrix]")
{
   auto testcase = GENERATE(TestCase::GenEigSPD, TestCase::GenEigGE,
                            TestCase::SVD);

   CAPTURE(TestCaseName(testcase));

   DenseMatrix M({{0.279841, 0.844288, 0.498302, 0.323955},
      {0.884680, 0.243511, 0.397405, 0.265708},
      {0.649685, 0.700754, 0.586396, 0.023724},
      {0.081588, 0.728236, 0.083123, 0.488041}
   });

   switch (testcase)
   {
      case TestCase::GenEigSPD:
      {
         DenseMatrix A({{0.56806, 0.29211, 0.48315, 0.70024},
            {0.29211, 0.85147, 0.68123, 0.70689},
            {0.48315, 0.68123, 1.07229, 1.02681},
            {0.70024, 0.70689, 1.02681, 1.15468}});

         DenseMatrixGeneralizedEigensystem geig(A,M,true,true);
         geig.Eval();
         Vector & Lambda = geig.EigenvaluesRealPart();
         DenseMatrix & V = geig.RightEigenvectors();
         DenseMatrix & W = geig.LeftEigenvectors();

         // check A * V - M * V * L
         DenseMatrix AV(4); Mult(A,V,AV);
         DenseMatrix MV(4); Mult(M,V,MV);
         MV.RightScaling(Lambda);  AV-=MV;

         REQUIRE(AV.MaxMaxNorm() == MFEM_Approx(0.));

         // check W^t * A - L * W^t * M
         DenseMatrix WtA(4); MultAtB(W,A,WtA);
         DenseMatrix WtM(4); MultAtB(W,M,WtM);
         WtM.LeftScaling(Lambda); WtA-=WtM;

         REQUIRE(WtA.MaxMaxNorm() == MFEM_Approx(0.));
      }
      break;
      case TestCase::GenEigGE:
      {
         DenseMatrix A({{0.486278, 0.041135, 0.480727, 0.616026},
            {0.523599, 0.119827, 0.087808, 0.415241},
            {0.214454, 0.661631, 0.909626, 0.744259},
            {0.107007, 0.630604, 0.077862, 0.221006}});

         DenseMatrixGeneralizedEigensystem geig(A,M,true,true);
         geig.Eval();
         Vector & Lambda_r = geig.EigenvaluesRealPart();
         Vector & Lambda_i = geig.EigenvaluesImagPart();
         DenseMatrix & V = geig.RightEigenvectors();

         DenseMatrix Vr(4), Vi(4);
         Vr.SetCol(0,V.GetColumn(0));
         Vr.SetCol(1,V.GetColumn(0));
         Vr.SetCol(2,V.GetColumn(2));
         Vr.SetCol(3,V.GetColumn(3));

         // Imag part of eigenvectors
         Vector vi(4); V.GetColumn(1,vi);
         Vi.SetCol(0,vi); vi *= -1.;
         Vi.SetCol(1,vi);
         Vi.SetCol(2,0.);
         Vi.SetCol(3,0.);

         // check A * V -  M * V * L
         // or  A * Vr = M*Vr * Lambda_r -  M*Vi * Lambda_i
         // and A * Vi = M*Vr * Lambda_i +  M*Vi * Lambda_r
         DenseMatrix AVr(4); Mult(A,Vr, AVr);
         DenseMatrix AVi(4); Mult(A,Vi, AVi);
         DenseMatrix MVr(4); Mult(M,Vr,MVr);
         DenseMatrix MVi(4); Mult(M,Vi,MVi);

         DenseMatrix MVrlr = MVr; MVrlr.RightScaling(Lambda_r);
         DenseMatrix MVrli = MVr; MVrli.RightScaling(Lambda_i);
         DenseMatrix MVilr = MVi; MVilr.RightScaling(Lambda_r);
         DenseMatrix MVili = MVi; MVili.RightScaling(Lambda_i);

         AVr -= MVrlr; AVr+= MVili;
         AVi -= MVrli; AVi-= MVilr;

         REQUIRE(AVr.MaxMaxNorm() == MFEM_Approx(0.));
         REQUIRE(AVi.MaxMaxNorm() == MFEM_Approx(0.));

         DenseMatrix & W = geig.LeftEigenvectors();
         DenseMatrix Wr(4), Wi(4);
         Wr.SetCol(0,W.GetColumn(0));
         Wr.SetCol(1,W.GetColumn(0));
         Wr.SetCol(2,W.GetColumn(2));
         Wr.SetCol(3,W.GetColumn(3));

         // Imag part of eigenvectors
         Vector wi(4); W.GetColumn(1,wi);
         Wi.SetCol(0,wi); wi *= -1.;
         Wi.SetCol(1,wi);
         Wi.SetCol(2,0.);
         Wi.SetCol(3,0.);

         // check W' * A - L * W' * M
         // or  Wr^t * A = Lambda_r * Wr^t * M + Lambda_i * Wi^t * M
         // and Wi^t * A = Lambda_r * Wi^t * M - Lambda_i * Wr^t * M
         DenseMatrix WrtA(4); MultAtB(Wr,A, WrtA);
         DenseMatrix WitA(4); MultAtB(Wi,A, WitA);
         DenseMatrix WrtM(4); MultAtB(Wr,M,WrtM);
         DenseMatrix WitM(4); MultAtB(Wi,M,WitM);

         DenseMatrix lrWrtM = WrtM; lrWrtM.LeftScaling(Lambda_r);
         DenseMatrix liWrtM = WrtM; liWrtM.LeftScaling(Lambda_i);
         DenseMatrix lrWitM = WitM; lrWitM.LeftScaling(Lambda_r);
         DenseMatrix liWitM = WitM; liWitM.LeftScaling(Lambda_i);

         WrtA -= lrWrtM; WrtA-= liWitM;
         WitA -= lrWitM; WitA+= liWrtM;

         REQUIRE(WrtA.MaxMaxNorm() == MFEM_Approx(0.));
         REQUIRE(WitA.MaxMaxNorm() == MFEM_Approx(0.));
      }
      break;
      case TestCase::SVD:
      {
         DenseMatrixSVD svd(M,true,true);
         svd.Eval(M);
         Vector &sigma = svd.Singularvalues();
         DenseMatrix &U = svd.LeftSingularvectors();
         DenseMatrix &V = svd.RightSingularvectors();

         DenseMatrix Vt(V); Vt.Transpose();
         DenseMatrix USVt(4); MultADBt(U,sigma,Vt,USVt);

         USVt -= M;

         REQUIRE(USVt.MaxMaxNorm() == MFEM_Approx(0.));
      }
      break;
   }
}

#endif // if MFEM_USE_LAPACK
