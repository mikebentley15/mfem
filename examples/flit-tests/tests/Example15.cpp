#include "mfem.hpp"
#include "mfem-helper.h"

#include <flit.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
  const std::string SEPARATOR = "#### MY SEPARATOR ####\n";

  using namespace std;
  using namespace mfem;

  // Choices for the problem setup. Affect bdr_func and rhs_func.
  int problem;
  int nfeatures;

  void UpdateProblem(Mesh &mesh, FiniteElementSpace &fespace,
                     GridFunction &x, BilinearForm &a, LinearForm &b)
  {
     FLIT_UNUSED(mesh);

     // Update the space: recalculate the number of DOFs and construct a matrix
     // that will adjust any GridFunctions to the new mesh state.
     fespace.Update();

     // Interpolate the solution on the new mesh by applying the transformation
     // matrix computed in the finite element space. Multiple GridFunctions could
     // be updated here.
     x.Update();

     // Free any transformation matrices to save memory.
     fespace.UpdatesFinished();

     // Inform the linear and bilinear forms that the space has changed.
     a.Update();
     b.Update();
  }

  const double alpha = 0.02;

  // Spherical front with a Gaussian cross section and radius t
  double front(double x, double y, double z, double t, int)
  {
     double r = sqrt(x*x + y*y + z*z);
     return exp(-0.5*pow((r - t)/alpha, 2));
  }

  double front_laplace(double x, double y, double z, double t, int dim)
  {
     double x2 = x*x, y2 = y*y, z2 = z*z, t2 = t*t;
     double r = sqrt(x2 + y2 + z2);
     double a2 = alpha*alpha, a4 = a2*a2;
     return -exp(-0.5*pow((r - t)/alpha, 2)) / a4 *
            (-2*t*(x2 + y2 + z2 - (dim-1)*a2/2)/r + x2 + y2 + z2 + t2 - dim*a2);
  }

  // Smooth spherical step function with radius t
  double ball(double x, double y, double z, double t, int)
  {
     double r = sqrt(x*x + y*y + z*z);
     return -atan(2*(r - t)/alpha);
  }

  double ball_laplace(double x, double y, double z, double t, int dim)
  {
     double x2 = x*x, y2 = y*y, z2 = z*z, t2 = 4*t*t;
     double r = sqrt(x2 + y2 + z2);
     double a2 = alpha*alpha;
     double den = pow(-a2 - 4*(x2 + y2 + z2 - 2*r*t) - t2, 2.0);
     return (dim == 2) ? 2*alpha*(a2 + t2 - 4*x2 - 4*y2)/r/den
            /*      */ : 4*alpha*(a2 + t2 - 4*r*t)/r/den;
  }

  // Composes several features into one function
  template<typename F0, typename F1>
  double composite_func(const Vector &pt, double t, F0 f0, F1 f1)
  {
     int dim = pt.Size();
     double x = pt(0), y = pt(1), z = 0.0;
     if (dim == 3) { z = pt(2); }

     if (problem == 0)
     {
        if (nfeatures <= 1)
        {
           return f0(x, y, z, t, dim);
        }
        else
        {
           double sum = 0.0;
           for (int i = 0; i < nfeatures; i++)
           {
              double x0 = 0.5*cos(2*M_PI * i / nfeatures);
              double y0 = 0.5*sin(2*M_PI * i / nfeatures);
              sum += f0(x - x0, y - y0, z, t, dim);
           }
           return sum;
        }
     }
     else
     {
        double sum = 0.0;
        for (int i = 0; i < nfeatures; i++)
        {
           double x0 = 0.5*cos(2*M_PI * i / nfeatures + M_PI*t);
           double y0 = 0.5*sin(2*M_PI * i / nfeatures + M_PI*t);
           sum += f1(x - x0, y - y0, z, 0.25, dim);
        }
        return sum;
     }
  }

  // Exact solution, used for the Dirichlet BC.
  double bdr_func(const Vector &pt, double t)
  {
     return composite_func(pt, t, front, ball);
  }

  // Laplace of the exact solution, used for the right hand side.
  double rhs_func(const Vector &pt, double t)
  {
     return composite_func(pt, t, front_laplace, ball_laplace);
  }
} // end of unnamed namespace

template <typename T>
class Example15 : public flit::TestBase<T> {
public:
  Example15(std::string id) : flit::TestBase<T>(std::move(id)) {}
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
flit::Variant Example15<double>::run_impl(const std::vector<double>& ti) {
   FLIT_UNUSED(ti);
   using namespace std;
   using namespace mfem;

   // 1. Use default command-line options.
   problem = 0;
   nfeatures = 1;
   const char *mesh_file = "../../data/star-hilbert.mesh";
   int order = 2;
   double t_final = 1.0;
   double max_elem_error = 5.0e-3;
   double hysteresis = 0.15; // derefinement safety coefficient
   int ref_levels = 0;
   int nc_limit = 3;         // maximum level of hanging nodes
   bool visualization = true;
   bool visit = false;

   // 2. Read the mesh from the given mesh file on all processors. We can handle
   //    triangular, quadrilateral, tetrahedral, hexahedral, surface and volume
   //    meshes with the same code.
   Mesh mesh(mesh_file, 1, 1);
   int dim = mesh.Dimension();
   int sdim = mesh.SpaceDimension();

   // 3. Project a NURBS mesh to a piecewise-quadratic curved mesh. Make sure
   //    that the mesh is non-conforming if it has quads or hexes and refine it.
   if (mesh.NURBSext)
   {
      mesh.UniformRefinement();
      if (ref_levels > 0) { ref_levels--; }
      mesh.SetCurvature(2);
   }
   mesh.EnsureNCMesh();
   for (int l = 0; l < ref_levels; l++)
   {
      mesh.UniformRefinement();
   }

   // 4. All boundary attributes will be used for essential (Dirichlet) BC.
   MFEM_VERIFY(mesh.bdr_attributes.Size() > 0,
               "Boundary attributes required in the mesh.");
   Array<int> ess_bdr(mesh.bdr_attributes.Max());
   ess_bdr = 1;

   // 5. Define a finite element space on the mesh. The polynomial order is one
   //    (linear) by default, but this can be changed on the command line.
   H1_FECollection fec(order, dim);
   FiniteElementSpace fespace(&mesh, &fec);

   // 6. As in Example 1p, we set up bilinear and linear forms corresponding to
   //    the Laplace problem -\Delta u = 1. We don't assemble the discrete
   //    problem yet, this will be done in the inner loop.
   BilinearForm a(&fespace);
   LinearForm b(&fespace);

   ConstantCoefficient one(1.0);
   FunctionCoefficient bdr(bdr_func);
   FunctionCoefficient rhs(rhs_func);

   BilinearFormIntegrator *integ = new DiffusionIntegrator(one);
   a.AddDomainIntegrator(integ);
   b.AddDomainIntegrator(new DomainLFIntegrator(rhs));

   // 7. The solution vector x and the associated finite element grid function
   //    will be maintained over the AMR iterations.
   GridFunction x(&fespace);

   // 8. Connect to GLVis. Prepare for VisIt output.
   char vishost[] = "localhost";
   int  visport   = 19916;

   socketstream sout;
   if (visualization)
   {
      sout.open(vishost, visport);
      if (!sout)
      {
         cout << "Unable to connect to GLVis server at "
              << vishost << ':' << visport << endl;
         cout << "GLVis visualization disabled.\n";
         visualization = false;
      }
      sout.precision(8);
   }

   VisItDataCollection visit_dc("Example15", &mesh);
   visit_dc.RegisterField("solution", &x);
   int vis_cycle = 0;

   // 9. As in Example 6, we set up a Zienkiewicz-Zhu estimator that will be
   //    used to obtain element error indicators. The integrator needs to
   //    provide the method ComputeElementFlux. The smoothed flux space is a
   //    vector valued H1 space here.
   FiniteElementSpace flux_fespace(&mesh, &fec, sdim);
   ZienkiewiczZhuEstimator estimator(*integ, x, flux_fespace);

   // 10. As in Example 6, we also need a refiner. This time the refinement
   //     strategy is based on a fixed threshold that is applied locally to each
   //     element. The global threshold is turned off by setting the total error
   //     fraction to zero. We also enforce a maximum refinement ratio between
   //     adjacent elements.
   ThresholdRefiner refiner(estimator);
   refiner.SetTotalErrorFraction(0.0); // use purely local threshold
   refiner.SetLocalErrorGoal(max_elem_error);
   refiner.PreferConformingRefinement();
   refiner.SetNCLimit(nc_limit);

   // 11. A derefiner selects groups of elements that can be coarsened to form
   //     a larger element. A conservative enough threshold needs to be set to
   //     prevent derefining elements that would immediately be refined again.
   ThresholdDerefiner derefiner(estimator);
   derefiner.SetThreshold(hysteresis * max_elem_error);
   derefiner.SetNCLimit(nc_limit);

   // 12. The outer time loop. In each iteration we update the right hand side,
   //     solve the problem on the current mesh, visualize the solution and
   //     refine the mesh as many times as necessary. Then we derefine any
   //     elements which have very small errors.
   x = 0.0;
   for (double time = 0.0; time < t_final + 1e-10; time += 0.01)
   {
      cout << "\nTime " << time << "\n\nRefinement:" << endl;

      // Set the current time in the coefficients.
      bdr.SetTime(time);
      rhs.SetTime(time);

      // Make sure errors will be recomputed in the following.
      refiner.Reset();
      derefiner.Reset();

      // 13. The inner refinement loop. At the end we want to have the current
      //     time step resolved to the prescribed tolerance in each element.
      for (int ref_it = 1; ; ref_it++)
      {
         cout << "Iteration: " << ref_it << ", number of unknowns: "
              << fespace.GetVSize() << endl;

         // 14. Recompute the field on the current mesh: assemble the stiffness
         //     matrix and the right-hand side.
         a.Assemble();
         b.Assemble();

         // 15. Project the exact solution to the essential boundary DOFs.
         x.ProjectBdrCoefficient(bdr, ess_bdr);

         // 16. Create and solve the linear system.
         Array<int> ess_tdof_list;
         fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

         SparseMatrix A;
         Vector B, X;
         a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

#ifndef MFEM_USE_SUITESPARSE
         GSSmoother M(A);
         PCG(A, M, B, X, 0, 500, 1e-12, 0.0);
#else
         UMFPackSolver umf_solver;
         umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
         umf_solver.SetOperator(A);
         umf_solver.Mult(B, X);
#endif

         // 17. Extract the local solution on each processor.
         a.RecoverFEMSolution(X, b, x);

         // 18. Send the solution by socket to a GLVis server and optionally
         //     save it in VisIt format.
         if (visualization)
         {
            sout.precision(8);
            sout << "solution\n" << mesh << x << flush;
         }
         if (visit)
         {
            visit_dc.SetCycle(vis_cycle++);
            visit_dc.SetTime(time);
            visit_dc.Save();
         }

         // 19. Apply the refiner on the mesh. The refiner calls the error
         //     estimator to obtain element errors, then it selects elements to
         //     be refined and finally it modifies the mesh. The Stop() method
         //     determines if all elements satisfy the local threshold.
         refiner.Apply(mesh);
         if (refiner.Stop())
         {
            break;
         }

         // 20. Update the space and interpolate the solution.
         UpdateProblem(mesh, fespace, x, a, b);
      }

      // 21. Use error estimates from the last inner iteration to check for
      //     possible derefinements. The derefiner works similarly as the
      //     refiner. The errors are not recomputed because the mesh did not
      //     change (and also the estimator was not Reset() at this time).
      if (derefiner.Apply(mesh))
      {
         cout << "\nDerefined elements." << endl;

         // 22. Update the space and interpolate the solution.
         UpdateProblem(mesh, fespace, x, a, b);
      }

      a.Update();
      b.Update();
   }

   // 22. Save the refined mesh and the solution. This output can be viewed later
   //     using GLVis: "glvis -m refined.mesh -g sol.gf".
   std::ostringstream out;
   out.precision(17);
   mesh.Print(out);
   out << SEPARATOR;
   x.Save(out);

   return out.str();
}

REGISTER_TYPE(Example15)
