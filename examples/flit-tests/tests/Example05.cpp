#include "mfem.hpp"
#include "mfem-helper.h"

#include <flit.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
const std::string SEPARATOR = "#### MY SEPARATOR ####\n";

using namespace mfem;

void uFun_ex(const Vector & x, Vector & u)
{
   double xi(x(0));
   double yi(x(1));
   double zi(0.0);
   if (x.Size() == 3)
   {
      zi = x(2);
   }

   u(0) = - exp(xi)*sin(yi)*cos(zi);
   u(1) = - exp(xi)*cos(yi)*cos(zi);

   if (x.Size() == 3)
   {
      u(2) = exp(xi)*sin(yi)*sin(zi);
   }
}

// Change if needed
double pFun_ex(const Vector & x)
{
   double xi(x(0));
   double yi(x(1));
   double zi(0.0);

   if (x.Size() == 3)
   {
      zi = x(2);
   }

   return exp(xi)*sin(yi)*cos(zi);
}

void fFun(const Vector & x, Vector & f)
{
   FLIT_UNUSED(x);
   f = 0.0;
}

double gFun(const Vector & x)
{
   if (x.Size() == 3)
   {
      return -pFun_ex(x);
   }
   else
   {
      return 0;
   }
}

double f_natural(const Vector & x)
{
   return (-pFun_ex(x));
}
} // end of unnamed namespace

template <typename T>
class Example05 : public flit::TestBase<T> {
public:
  Example05(std::string id) : flit::TestBase<T>(std::move(id)) {}
  virtual size_t getInputsPerRun() override { return 0; }
  virtual std::vector<T> getDefaultInput() override { return {}; }

  virtual long double compare(const std::string &ground_truth,
                              const std::string &test_results) const override {
    FLIT_UNUSED(ground_truth);
    FLIT_UNUSED(test_results);

    auto gt = load_sim(ground_truth, SEPARATOR);
    auto res = load_sim(test_results, SEPARATOR);

    // TODO: use sol2 as well.  gt.sol is the u solution, gt.sol2 is the p solution
    auto& diff = gt.sol;
    diff -= res.sol;

    return diff.Norml2();
  }

protected:
  // Default implementation does nothing
  virtual flit::Variant run_impl(const std::vector<T>& ti) override {
    FLIT_UNUSED(ti);
    return flit::Variant();
  }

protected:
  using flit::TestBase<T>::id;
};

template<>
flit::Variant Example05<double>::run_impl(const std::vector<double>& ti) {
   FLIT_UNUSED(ti);
   using namespace mfem;

   StopWatch chrono;

   // 1. Use defulat command-line options
   const char *mesh_file = "../../data/star.mesh";
   int order = 1;
   bool visualization = false;

   // 2. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
   //    the same code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();

   // 3. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
   //    largest number that gives a final mesh with no more than 10,000
   //    elements.
   {
      int ref_levels =
         (int)floor(log(10000./mesh->GetNE())/log(2.)/dim);
      for (int l = 0; l < ref_levels; l++)
      {
         mesh->UniformRefinement();
      }
   }

   // 4. Define a finite element space on the mesh. Here we use the
   //    Raviart-Thomas finite elements of the specified order.
   FiniteElementCollection *hdiv_coll(new RT_FECollection(order, dim));
   FiniteElementCollection *l2_coll(new L2_FECollection(order, dim));

   FiniteElementSpace *R_space = new FiniteElementSpace(mesh, hdiv_coll);
   FiniteElementSpace *W_space = new FiniteElementSpace(mesh, l2_coll);

   // 5. Define the BlockStructure of the problem, i.e. define the array of
   //    offsets for each variable. The last component of the Array is the sum
   //    of the dimensions of each block.
   Array<int> block_offsets(3); // number of variables + 1
   block_offsets[0] = 0;
   block_offsets[1] = R_space->GetVSize();
   block_offsets[2] = W_space->GetVSize();
   block_offsets.PartialSum();

   std::cout << "***********************************************************\n";
   std::cout << "dim(R) = " << block_offsets[1] - block_offsets[0] << "\n";
   std::cout << "dim(W) = " << block_offsets[2] - block_offsets[1] << "\n";
   std::cout << "dim(R+W) = " << block_offsets.Last() << "\n";
   std::cout << "***********************************************************\n";

   // 6. Define the coefficients, analytical solution, and rhs of the PDE.
   ConstantCoefficient k(1.0);

   VectorFunctionCoefficient fcoeff(dim, fFun);
   FunctionCoefficient fnatcoeff(f_natural);
   FunctionCoefficient gcoeff(gFun);

   VectorFunctionCoefficient ucoeff(dim, uFun_ex);
   FunctionCoefficient pcoeff(pFun_ex);

   // 7. Allocate memory (x, rhs) for the analytical solution and the right hand
   //    side.  Define the GridFunction u,p for the finite element solution and
   //    linear forms fform and gform for the right hand side.  The data
   //    allocated by x and rhs are passed as a reference to the grid functions
   //    (u,p) and the linear forms (fform, gform).
   BlockVector x(block_offsets), rhs(block_offsets);

   LinearForm *fform(new LinearForm);
   fform->Update(R_space, rhs.GetBlock(0), 0);
   fform->AddDomainIntegrator(new VectorFEDomainLFIntegrator(fcoeff));
   fform->AddBoundaryIntegrator(new VectorFEBoundaryFluxLFIntegrator(fnatcoeff));
   fform->Assemble();

   LinearForm *gform(new LinearForm);
   gform->Update(W_space, rhs.GetBlock(1), 0);
   gform->AddDomainIntegrator(new DomainLFIntegrator(gcoeff));
   gform->Assemble();

   // 8. Assemble the finite element matrices for the Darcy operator
   //
   //                            D = [ M  B^T ]
   //                                [ B   0  ]
   //     where:
   //
   //     M = \int_\Omega k u_h \cdot v_h d\Omega   u_h, v_h \in R_h
   //     B   = -\int_\Omega \div u_h q_h d\Omega   u_h \in R_h, q_h \in W_h
   BilinearForm *mVarf(new BilinearForm(R_space));
   MixedBilinearForm *bVarf(new MixedBilinearForm(R_space, W_space));

   mVarf->AddDomainIntegrator(new VectorFEMassIntegrator(k));
   mVarf->Assemble();
   mVarf->Finalize();
   SparseMatrix &M(mVarf->SpMat());

   bVarf->AddDomainIntegrator(new VectorFEDivergenceIntegrator);
   bVarf->Assemble();
   bVarf->Finalize();
   SparseMatrix & B(bVarf->SpMat());
   B *= -1.;
   SparseMatrix *BT = Transpose(B);

   BlockMatrix darcyMatrix(block_offsets);
   darcyMatrix.SetBlock(0,0, &M);
   darcyMatrix.SetBlock(0,1, BT);
   darcyMatrix.SetBlock(1,0, &B);

   // 9. Construct the operators for preconditioner
   //
   //                 P = [ diag(M)         0         ]
   //                     [  0       B diag(M)^-1 B^T ]
   //
   //     Here we use Symmetric Gauss-Seidel to approximate the inverse of the
   //     pressure Schur Complement
   SparseMatrix *MinvBt = Transpose(B);
   Vector Md(M.Height());
   M.GetDiag(Md);
   for (int i = 0; i < Md.Size(); i++)
   {
      MinvBt->ScaleRow(i, 1./Md(i));
   }
   SparseMatrix *S = Mult(B, *MinvBt);

   Solver *invM, *invS;
   invM = new DSmoother(M);
#ifndef MFEM_USE_SUITESPARSE
   invS = new GSSmoother(*S);
#else
   invS = new UMFPackSolver(*S);
#endif

   invM->iterative_mode = false;
   invS->iterative_mode = false;

   BlockDiagonalPreconditioner darcyPrec(block_offsets);
   darcyPrec.SetDiagonalBlock(0, invM);
   darcyPrec.SetDiagonalBlock(1, invS);

   // 10. Solve the linear system with MINRES.
   //     Check the norm of the unpreconditioned residual.
   int maxIter(1000);
   double rtol(1.e-6);
   double atol(1.e-10);

   chrono.Clear();
   chrono.Start();
   MINRESSolver solver;
   solver.SetAbsTol(atol);
   solver.SetRelTol(rtol);
   solver.SetMaxIter(maxIter);
   solver.SetOperator(darcyMatrix);
   solver.SetPreconditioner(darcyPrec);
   solver.SetPrintLevel(1);
   x = 0.0;
   solver.Mult(rhs, x);
   chrono.Stop();

   if (solver.GetConverged())
      std::cout << "MINRES converged in " << solver.GetNumIterations()
                << " iterations with a residual norm of " << solver.GetFinalNorm() << ".\n";
   else
      std::cout << "MINRES did not converge in " << solver.GetNumIterations()
                << " iterations. Residual norm is " << solver.GetFinalNorm() << ".\n";
   std::cout << "MINRES solver took " << chrono.RealTime() << "s. \n";

   // 11. Create the grid functions u and p. Compute the L2 error norms.
   GridFunction u, p;
   u.MakeRef(R_space, x.GetBlock(0), 0);
   p.MakeRef(W_space, x.GetBlock(1), 0);

   int order_quad = std::max(2, 2*order+1);
   const IntegrationRule *irs[Geometry::NumGeom];
   for (int i=0; i < Geometry::NumGeom; ++i)
   {
      irs[i] = &(IntRules.Get(i, order_quad));
   }

   double err_u  = u.ComputeL2Error(ucoeff, irs);
   double norm_u = ComputeLpNorm(2., ucoeff, *mesh, irs);
   double err_p  = p.ComputeL2Error(pcoeff, irs);
   double norm_p = ComputeLpNorm(2., pcoeff, *mesh, irs);

   std::cout << "|| u_h - u_ex || / || u_ex || = " << err_u / norm_u << "\n";
   std::cout << "|| p_h - p_ex || / || p_ex || = " << err_p / norm_p << "\n";

   // 12. Save the mesh and the solution. This output can be viewed later using
   //     GLVis: "glvis -m ex5.mesh -g sol_u.gf" or "glvis -m ex5.mesh -g
   //     sol_p.gf".
   std::ostringstream out;
   {
      out.precision(17);
      mesh->Print(out);
      out << SEPARATOR;
      u.Save(out);
      out << SEPARATOR;
      p.Save(out);
   }

   // 14. Send the solution by socket to a GLVis server.
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream u_sock(vishost, visport);
      u_sock.precision(8);
      u_sock << "solution\n" << *mesh << u << "window_title 'Velocity'\n";
      socketstream p_sock(vishost, visport);
      p_sock.precision(8);
      p_sock << "solution\n" << *mesh << p << "window_title 'Pressure'\n";
   }


   // 15. Free the used memory.
   delete fform;
   delete gform;
   delete invM;
   delete invS;
   delete S;
   delete MinvBt;
   delete BT;
   delete mVarf;
   delete bVarf;
   delete W_space;
   delete R_space;
   delete l2_coll;
   delete hdiv_coll;
   delete mesh;

   return out.str();
}

REGISTER_TYPE(Example05)
