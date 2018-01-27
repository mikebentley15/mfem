#include "mfem.hpp"
#include "mfem-helper.h"

#include <flit.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
  const std::string SEPARATOR = "#### MY SEPARATOR ####\n";
} // end of unnamed namespace

template <typename T>
class Example14 : public flit::TestBase<T> {
public:
  Example14(std::string id) : flit::TestBase<T>(std::move(id)) {}
  virtual size_t getInputsPerRun() override { return 0; }
  virtual std::vector<T> getDefaultInput() override { return {}; }

  virtual long double compare(const std::string &ground_truth,
                              const std::string &test_results) const override {
    FLIT_UNUSED(ground_truth);
    FLIT_UNUSED(test_results);

    auto gt = load_sim(ground_truth, SEPARATOR);
    auto res = load_sim(test_results, SEPARATOR);

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
flit::Variant Example14<double>::run_impl(const std::vector<double>& ti) {
   FLIT_UNUSED(ti);
   using namespace std;
   using namespace mfem;

   // 1. Use defulat command-line options
   const char *mesh_file = "../../data/star.mesh";
   int ref_levels = -1;
   int order = 1;
   double sigma = -1.0;
   double kappa = -1.0;
   bool visualization = false;
   if (kappa < 0)
   {
      kappa = (order+1) * (order+1);
   }

   // 2. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral and hexahedral meshes with the same code.
   //    NURBS meshes are projected to second order meshes.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();

   // 3. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement. By default, or if ref_levels < 0,
   //    we choose it to be the largest number that gives a final mesh with no
   //    more than 50,000 elements.
   {
      if (ref_levels < 0)
      {
         ref_levels = (int)floor(log(50000./mesh->GetNE())/log(2.)/dim);
      }
      for (int l = 0; l < ref_levels; l++)
      {
         mesh->UniformRefinement();
      }
   }
   if (mesh->NURBSext)
   {
      mesh->SetCurvature(max(order, 1));
   }

   // 4. Define a finite element space on the mesh. Here we use discontinuous
   //    finite elements of the specified order >= 0.
   FiniteElementCollection *fec = new DG_FECollection(order, dim);
   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);
   cout << "Number of unknowns: " << fespace->GetVSize() << endl;

   // 5. Set up the linear form b(.) which corresponds to the right-hand side of
   //    the FEM linear system.
   LinearForm *b = new LinearForm(fespace);
   ConstantCoefficient one(1.0);
   ConstantCoefficient zero(0.0);
   b->AddDomainIntegrator(new DomainLFIntegrator(one));
   b->AddBdrFaceIntegrator(
      new DGDirichletLFIntegrator(zero, one, sigma, kappa));
   b->Assemble();

   // 6. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero.
   GridFunction x(fespace);
   x = 0.0;

   // 7. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //    domain integrator and the interior and boundary DG face integrators.
   //    Note that boundary conditions are imposed weakly in the form, so there
   //    is no need for dof elimination. After assembly and finalizing we
   //    extract the corresponding sparse matrix A.
   BilinearForm *a = new BilinearForm(fespace);
   a->AddDomainIntegrator(new DiffusionIntegrator(one));
   a->AddInteriorFaceIntegrator(new DGDiffusionIntegrator(one, sigma, kappa));
   a->AddBdrFaceIntegrator(new DGDiffusionIntegrator(one, sigma, kappa));
   a->Assemble();
   a->Finalize();
   const SparseMatrix &A = a->SpMat();

#ifndef MFEM_USE_SUITESPARSE
   // 8. Define a simple symmetric Gauss-Seidel preconditioner and use it to
   //    solve the system Ax=b with PCG in the symmetric case, and GMRES in the
   //    non-symmetric one.
   GSSmoother M(A);
   if (sigma == -1.0)
   {
      PCG(A, M, *b, x, 1, 500, 1e-12, 0.0);
   }
   else
   {
      GMRES(A, M, *b, x, 1, 500, 10, 1e-12, 0.0);
   }
#else
   // 8. If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
   UMFPackSolver umf_solver;
   umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
   umf_solver.SetOperator(A);
   umf_solver.Mult(*b, x);
#endif

   // 9. Save the refined mesh and the solution. This output can be viewed later
   //    using GLVis: "glvis -m refined.mesh -g sol.gf".
   std::ostringstream out;
   {
      out.precision(17);
      mesh->Print(out);
      out << SEPARATOR;
      x.Save(out);
   }

   // 10. Send the solution by socket to a GLVis server.
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << *mesh << x << flush;
   }

   // 11. Free the used memory.
   delete a;
   delete b;
   delete fespace;
   delete fec;
   delete mesh;

   return out.str();
}

REGISTER_TYPE(Example14)
