#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/numerics/matrix_creator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <cmath>
#include <iostream>

namespace Step35
{
  using namespace dealii;

  // Since our method has several parameters that can be fine-tuned we put them
  // into an external file, so that they can be determined at run-time.
  //
  // This includes, in particular, the formulation of the equation for the
  // auxiliary variable $\phi$, for which we declare an <code>enum</code>. Next,
  // we declare a class that is going to read and store all the parameters that
  // our program needs to run.
  namespace RunTimeParameters
  {
    enum class Method
    {
      standard,
      rotational
    };

    class Data_Storage
    {
    public:
      Data_Storage();

      void read_data(const std::string &filename);

      Method form;

      double dt;
      double initial_time;
      double final_time;

      double Reynolds;

      unsigned int n_global_refines;

      unsigned int pressure_degree;

      unsigned int vel_max_iterations;
      unsigned int vel_Krylov_size;
      unsigned int vel_off_diagonals;
      unsigned int vel_update_prec;
      double       vel_eps;
      double       vel_diag_strength;

      bool         verbose;
      unsigned int output_interval;

    protected:
      ParameterHandler prm;
    };

    // In the constructor of this class we declare all the parameters.
    Data_Storage::Data_Storage()
      : form(Method::rotational)
      , dt(5e-4)
      , initial_time(0.)
      , final_time(1.)
      , Reynolds(1.)
      , n_global_refines(0)
      , pressure_degree(1)
      , vel_max_iterations(1000)
      , vel_Krylov_size(30)
      , vel_off_diagonals(60)
      , vel_update_prec(15)
      , vel_eps(1e-12)
      , vel_diag_strength(0.01)
      , verbose(true)
      , output_interval(15)
    {
      prm.declare_entry("Method_Form",
                        "rotational",
                        Patterns::Selection("rotational|standard"),
                        " Used to select the type of method that we are going "
                        "to use. ");
      prm.enter_subsection("Physical data");
      {
        prm.declare_entry("initial_time",
                          "0.",
                          Patterns::Double(0.),
                          " The initial time of the simulation. ");
        prm.declare_entry("final_time",
                          "1.",
                          Patterns::Double(0.),
                          " The final time of the simulation. ");
        prm.declare_entry("Reynolds",
                          "1.",
                          Patterns::Double(0.),
                          " The Reynolds number. ");
      }
      prm.leave_subsection();

      prm.enter_subsection("Time step data");
      {
        prm.declare_entry("dt",
                          "5e-4",
                          Patterns::Double(0.),
                          " The time step size. ");
      }
      prm.leave_subsection();

      prm.enter_subsection("Space discretization");
      {
        prm.declare_entry("n_of_refines",
                          "0",
                          Patterns::Integer(0, 15),
                          " The number of global refines we do on the mesh. ");
        prm.declare_entry("pressure_fe_degree",
                          "1",
                          Patterns::Integer(1, 5),
                          " The polynomial degree for the pressure space. ");
      }
      prm.leave_subsection();

      prm.enter_subsection("Data solve velocity");
      {
        prm.declare_entry(
          "max_iterations",
          "1000",
          Patterns::Integer(1, 1000),
          " The maximal number of iterations GMRES must make. ");
        prm.declare_entry("eps",
                          "1e-12",
                          Patterns::Double(0.),
                          " The stopping criterion. ");
        prm.declare_entry("Krylov_size",
                          "30",
                          Patterns::Integer(1),
                          " The size of the Krylov subspace to be used. ");
        prm.declare_entry("off_diagonals",
                          "60",
                          Patterns::Integer(0),
                          " The number of off-diagonal elements ILU must "
                          "compute. ");
        prm.declare_entry("diag_strength",
                          "0.01",
                          Patterns::Double(0.),
                          " Diagonal strengthening coefficient. ");
        prm.declare_entry("update_prec",
                          "15",
                          Patterns::Integer(1),
                          " This number indicates how often we need to "
                          "update the preconditioner");
      }
      prm.leave_subsection();

      prm.declare_entry("verbose",
                        "true",
                        Patterns::Bool(),
                        " This indicates whether the output of the solution "
                        "process should be verbose. ");

      prm.declare_entry("output_interval",
                        "1",
                        Patterns::Integer(1),
                        " This indicates between how many time steps we print "
                        "the solution. ");
    }

    void Data_Storage::read_data(const std::string &filename)
    {
      std::ifstream file(filename);
      AssertThrow(file, ExcFileNotOpen(filename));

      prm.parse_input(file);

      if (prm.get("Method_Form") == std::string("rotational"))
        form = Method::rotational;
      else
        form = Method::standard;

      prm.enter_subsection("Physical data");
      {
        initial_time = prm.get_double("initial_time");
        final_time   = prm.get_double("final_time");
        Reynolds     = prm.get_double("Reynolds");
      }
      prm.leave_subsection();

      prm.enter_subsection("Time step data");
      {
        dt = prm.get_double("dt");
      }
      prm.leave_subsection();

      prm.enter_subsection("Space discretization");
      {
        n_global_refines = prm.get_integer("n_of_refines");
        pressure_degree  = prm.get_integer("pressure_fe_degree");
      }
      prm.leave_subsection();

      prm.enter_subsection("Data solve velocity");
      {
        vel_max_iterations = prm.get_integer("max_iterations");
        vel_eps            = prm.get_double("eps");
        vel_Krylov_size    = prm.get_integer("Krylov_size");
        vel_off_diagonals  = prm.get_integer("off_diagonals");
        vel_diag_strength  = prm.get_double("diag_strength");
        vel_update_prec    = prm.get_integer("update_prec");
      }
      prm.leave_subsection();

      verbose = prm.get_bool("verbose");

      output_interval = prm.get_integer("output_interval");
    }
  } // namespace RunTimeParameters

  // In the next namespace, we declare the initial and boundary conditions:
  namespace EquationData
  {
    // As we have chosen a completely decoupled formulation, we will not take
    // advantage of deal.II's capabilities to handle vector valued problems. We
    // do, however, want to use an interface for the equation data that is
    // somehow dimension independent. To be able to do that, our functions
    // should be able to know on which spatial component we are currently
    // working, and we should be able to have a common interface to do that. The
    // following class is an attempt in that direction.
    template <int dim>
    class MultiComponentFunction : public Function<dim>
    {
    public:
      MultiComponentFunction(const double initial_time = 0.);
      void set_component(const unsigned int d);

    protected:
      unsigned int comp;
    };

    template <int dim>
    MultiComponentFunction<dim>::MultiComponentFunction(
      const double initial_time)
      : Function<dim>(1, initial_time)
      , comp(0)
    {}

    template <int dim>
    void MultiComponentFunction<dim>::set_component(const unsigned int d)
    {
      Assert(d < dim, ExcIndexRange(d, 0, dim));
      comp = d;
    }

    // With this class defined, we declare classes that describe the boundary
    // conditions for velocity and pressure:
    template <int dim>
    class Velocity : public MultiComponentFunction<dim>
    {
    public:
      Velocity(const double initial_time = 0.0);

      virtual double value(const Point<dim>  &p,
                           const unsigned int component = 0) const override;

      virtual void value_list(const std::vector<Point<dim>> &points,
                              std::vector<double>           &values,
                              const unsigned int component = 0) const override;
    };

    template <int dim>
    Velocity<dim>::Velocity(const double initial_time)
      : MultiComponentFunction<dim>(initial_time)
    {}

    template <int dim>
    void Velocity<dim>::value_list(const std::vector<Point<dim>> &points,
                                   std::vector<double>           &values,
                                   const unsigned int) const
    {
      const unsigned int n_points = points.size();
      AssertDimension(values.size(), n_points);
      for (unsigned int i = 0; i < n_points; ++i)
        values[i] = Velocity<dim>::value(points[i]);
    }

    template <int dim>
    double Velocity<dim>::value(const Point<dim> &p, const unsigned int) const
    {
      if (this->comp == 0)
        {
          const double Um = 1.5;
          const double H  = 4.1;
          return 4. * Um * p[1] * (H - p[1]) / (H * H);
        }
      else
        return 0.;
    }

    template <int dim>
    class Pressure : public Function<dim>
    {
    public:
      Pressure(const double initial_time = 0.0);

      virtual double value(const Point<dim>  &p,
                           const unsigned int component = 0) const override;

      virtual void value_list(const std::vector<Point<dim>> &points,
                              std::vector<double>           &values,
                              const unsigned int component = 0) const override;
    };

    template <int dim>
    Pressure<dim>::Pressure(const double initial_time)
      : Function<dim>(1, initial_time)
    {}

    template <int dim>
    double Pressure<dim>::value(const Point<dim>  &p,
                                const unsigned int component) const
    {
      (void)component;
      AssertIndexRange(component, 1);
      return 25. - p[0];
    }

    template <int dim>
    void Pressure<dim>::value_list(const std::vector<Point<dim>> &points,
                                   std::vector<double>           &values,
                                   const unsigned int component) const
    {
      (void)component;
      AssertIndexRange(component, 1);
      const unsigned int n_points = points.size();
      AssertDimension(values.size(), n_points);
      for (unsigned int i = 0; i < n_points; ++i)
        values[i] = Pressure<dim>::value(points[i]);
    }
  } // namespace EquationData

  // Now for the main class of the program. It implements the various versions
  // of the projection method for Navier-Stokes equations.
  template <int dim>
  class NavierStokesProjection
  {
  public:
    NavierStokesProjection(const RunTimeParameters::Data_Storage &data);

    void run(const bool verbose = false, const unsigned int n_plots = 10);

  protected:
    RunTimeParameters::Method type;

    const unsigned int deg;
    double             dt;
    const double       t_0;
    const double       T;
    const double       Re;

    EquationData::Velocity<dim>               vel_exact;
    std::map<types::global_dof_index, double> boundary_values;
    std::vector<types::boundary_id>           boundary_ids;

    Triangulation<dim> triangulation;

    const FE_Q<dim> fe_velocity;
    const FE_Q<dim> fe_pressure;

    DoFHandler<dim> dof_handler_velocity;
    DoFHandler<dim> dof_handler_pressure;

    const QGauss<dim> quadrature_pressure;
    const QGauss<dim> quadrature_velocity;

    SparsityPattern sparsity_pattern_velocity;
    SparsityPattern sparsity_pattern_pressure;
    SparsityPattern sparsity_pattern_pres_vel;

    SparseMatrix<double> vel_Laplace_plus_Mass;
    SparseMatrix<double> vel_it_matrix[dim];
    SparseMatrix<double> vel_Mass;
    SparseMatrix<double> vel_Laplace;
    SparseMatrix<double> vel_Advection;
    SparseMatrix<double> pres_Laplace;
    SparseMatrix<double> pres_Mass;
    SparseMatrix<double> pres_Diff[dim];
    SparseMatrix<double> pres_iterative;

    Vector<double> pres_n;
    Vector<double> pres_n_minus_1;
    Vector<double> phi_n;
    Vector<double> phi_n_minus_1;
    Vector<double> u_n[dim];
    Vector<double> u_n_minus_1[dim];
    Vector<double> u_star[dim];
    Vector<double> force[dim];
    Vector<double> v_tmp;
    Vector<double> pres_tmp;
    Vector<double> rot_u;

    SparseILU<double>   prec_velocity[dim];
    SparseILU<double>   prec_pres_Laplace;
    SparseDirectUMFPACK prec_mass;
    SparseDirectUMFPACK prec_vel_mass;

    DeclException2(ExcInvalidTimeStep,
                   double,
                   double,
                   << " The time step " << arg1 << " is out of range."
                   << std::endl
                   << " The permitted range is (0," << arg2 << ']');

    void create_triangulation_and_dofs(const unsigned int n_refines);

    void initialize();

    void interpolate_velocity();

    void diffusion_step(const bool reinit_prec);

    void projection_step(const bool reinit_prec);

    void update_pressure(const bool reinit_prec);

    // Function to estimate the time step
    double estimate_time_step(const Vector<double> u_n[dim], const Vector<double> u_n_minus_1[dim], const double &dt);

    // Function to save the time step values in a csv file
    void save_time_steps_to_csv(const std::string &filename) const;

  private:
    unsigned int vel_max_its;
    unsigned int vel_Krylov_size;
    unsigned int vel_off_diagonals;
    unsigned int vel_update_prec;
    double       vel_eps;
    double       vel_diag_strength;

    std::vector<double> dt_values; // This is used to store the time step values

    void initialize_velocity_matrices();

    void initialize_pressure_matrices();

    using IteratorTuple =
      std::tuple<typename DoFHandler<dim>::active_cell_iterator,
                 typename DoFHandler<dim>::active_cell_iterator>;

    using IteratorPair = SynchronousIterators<IteratorTuple>;

    void initialize_gradient_operator();

    struct InitGradPerTaskData
    {
      unsigned int                         d;
      unsigned int                         vel_dpc;
      unsigned int                         pres_dpc;
      FullMatrix<double>                   local_grad;
      std::vector<types::global_dof_index> vel_local_dof_indices;
      std::vector<types::global_dof_index> pres_local_dof_indices;

      InitGradPerTaskData(const unsigned int dd,
                          const unsigned int vdpc,
                          const unsigned int pdpc)
        : d(dd)
        , vel_dpc(vdpc)
        , pres_dpc(pdpc)
        , local_grad(vdpc, pdpc)
        , vel_local_dof_indices(vdpc)
        , pres_local_dof_indices(pdpc)
      {}
    };

    struct InitGradScratchData
    {
      unsigned int  nqp;
      FEValues<dim> fe_val_vel;
      FEValues<dim> fe_val_pres;
      InitGradScratchData(const FE_Q<dim>   &fe_v,
                          const FE_Q<dim>   &fe_p,
                          const QGauss<dim> &quad,
                          const UpdateFlags  flags_v,
                          const UpdateFlags  flags_p)
        : nqp(quad.size())
        , fe_val_vel(fe_v, quad, flags_v)
        , fe_val_pres(fe_p, quad, flags_p)
      {}
      InitGradScratchData(const InitGradScratchData &data)
        : nqp(data.nqp)
        , fe_val_vel(data.fe_val_vel.get_fe(),
                     data.fe_val_vel.get_quadrature(),
                     data.fe_val_vel.get_update_flags())
        , fe_val_pres(data.fe_val_pres.get_fe(),
                      data.fe_val_pres.get_quadrature(),
                      data.fe_val_pres.get_update_flags())
      {}
    };

    void assemble_one_cell_of_gradient(const IteratorPair  &SI,
                                       InitGradScratchData &scratch,
                                       InitGradPerTaskData &data);

    void copy_gradient_local_to_global(const InitGradPerTaskData &data);

    void assemble_advection_term();

    struct AdvectionPerTaskData
    {
      FullMatrix<double>                   local_advection;
      std::vector<types::global_dof_index> local_dof_indices;
      AdvectionPerTaskData(const unsigned int dpc)
        : local_advection(dpc, dpc)
        , local_dof_indices(dpc)
      {}
    };

    struct AdvectionScratchData
    {
      unsigned int                nqp;
      unsigned int                dpc;
      std::vector<Point<dim>>     u_star_local;
      std::vector<Tensor<1, dim>> grad_u_star;
      std::vector<double>         u_star_tmp;
      FEValues<dim>               fe_val;
      AdvectionScratchData(const FE_Q<dim>   &fe,
                           const QGauss<dim> &quad,
                           const UpdateFlags  flags)
        : nqp(quad.size())
        , dpc(fe.n_dofs_per_cell())
        , u_star_local(nqp)
        , grad_u_star(nqp)
        , u_star_tmp(nqp)
        , fe_val(fe, quad, flags)
      {}

      AdvectionScratchData(const AdvectionScratchData &data)
        : nqp(data.nqp)
        , dpc(data.dpc)
        , u_star_local(nqp)
        , grad_u_star(nqp)
        , u_star_tmp(nqp)
        , fe_val(data.fe_val.get_fe(),
                 data.fe_val.get_quadrature(),
                 data.fe_val.get_update_flags())
      {}
    };

    void assemble_one_cell_of_advection(
      const typename DoFHandler<dim>::active_cell_iterator &cell,
      AdvectionScratchData                                 &scratch,
      AdvectionPerTaskData                                 &data);

    void copy_advection_local_to_global(const AdvectionPerTaskData &data);

    // The final few functions implement the diffusion solve as well as
    // postprocessing the output, including computing the curl of the velocity:
    void diffusion_component_solve(const unsigned int d);

    void output_results(const unsigned int step);

    void assemble_vorticity(const bool reinit_prec);
  };

  // In the constructor, we just read all the data from the
  // <code>Data_Storage</code> object that is passed as an argument, verify that
  // the data we read is reasonable and, finally, create the triangulation and
  // load the initial data.
  template <int dim>
  NavierStokesProjection<dim>::NavierStokesProjection(
    const RunTimeParameters::Data_Storage &data)
    : type(data.form)
    , deg(data.pressure_degree)
    , dt(data.dt)
    , t_0(data.initial_time)
    , T(data.final_time)
    , Re(data.Reynolds)
    , vel_exact(data.initial_time)
    , fe_velocity(deg + 1)
    , fe_pressure(deg)
    , dof_handler_velocity(triangulation)
    , dof_handler_pressure(triangulation)
    , quadrature_pressure(deg + 1)
    , quadrature_velocity(deg + 2)
    , vel_max_its(data.vel_max_iterations)
    , vel_Krylov_size(data.vel_Krylov_size)
    , vel_off_diagonals(data.vel_off_diagonals)
    , vel_update_prec(data.vel_update_prec)
    , vel_eps(data.vel_eps)
    , vel_diag_strength(data.vel_diag_strength)
  {
    if (deg < 1)
      std::cout
        << " WARNING: The chosen pair of finite element spaces is not stable."
        << std::endl
        << " The obtained results will be nonsense" << std::endl;

    AssertThrow(!((dt <= 0.) || (dt > .5 * T)), ExcInvalidTimeStep(dt, .5 * T));

    create_triangulation_and_dofs(data.n_global_refines);
    initialize();
  }

  // The method that creates the triangulation and refines it the needed number
  // of times. After creating the triangulation, it creates the mesh dependent
  // data, i.e. it distributes degrees of freedom and renumbers them, and
  // initializes the matrices and vectors that we will use.
  template <int dim>
  void NavierStokesProjection<dim>::create_triangulation_and_dofs(
    const unsigned int n_refines)
  {

    // For flow past a cylinder
    // GridGenerator::channel_with_cylinder(triangulation, 0.03, 2, 2.0, true);

    GridIn<dim> grid_in;
    grid_in.attach_triangulation(triangulation);

    {
      std::string   filename = "nsbench2.inp";
      std::ifstream file(filename);
      Assert(file, ExcFileNotOpen(filename));
      grid_in.read_ucd(file);
    }

    std::cout << "Number of refines = " << n_refines << std::endl;
    triangulation.refine_global(n_refines);
    std::cout << "Number of active cells: " << triangulation.n_active_cells()
              << std::endl;

    boundary_ids = triangulation.get_boundary_ids();

    dof_handler_velocity.distribute_dofs(fe_velocity);
    DoFRenumbering::boost::Cuthill_McKee(dof_handler_velocity);
    dof_handler_pressure.distribute_dofs(fe_pressure);
    DoFRenumbering::boost::Cuthill_McKee(dof_handler_pressure);

    initialize_velocity_matrices();
    initialize_pressure_matrices();
    initialize_gradient_operator();

    pres_n.reinit(dof_handler_pressure.n_dofs());
    pres_n_minus_1.reinit(dof_handler_pressure.n_dofs());
    phi_n.reinit(dof_handler_pressure.n_dofs());
    phi_n_minus_1.reinit(dof_handler_pressure.n_dofs());
    pres_tmp.reinit(dof_handler_pressure.n_dofs());
    for (unsigned int d = 0; d < dim; ++d)
      {
        u_n[d].reinit(dof_handler_velocity.n_dofs());
        u_n_minus_1[d].reinit(dof_handler_velocity.n_dofs());
        u_star[d].reinit(dof_handler_velocity.n_dofs());
        force[d].reinit(dof_handler_velocity.n_dofs());
      }
    v_tmp.reinit(dof_handler_velocity.n_dofs());
    rot_u.reinit(dof_handler_velocity.n_dofs());

    std::cout << "dim (X_h) = " << (dof_handler_velocity.n_dofs() * dim) //
              << std::endl                                               //
              << "dim (M_h) = " << dof_handler_pressure.n_dofs()         //
              << std::endl                                               //
              << "Re        = " << Re << std::endl                       //
              << std::endl;
  }

  // This method creates the constant matrices and loads the initial data
  template <int dim>
  void NavierStokesProjection<dim>::initialize()
  {
    vel_Laplace_plus_Mass = 0.;
    vel_Laplace_plus_Mass.add(1. / Re, vel_Laplace);
    vel_Laplace_plus_Mass.add(1.5 / dt, vel_Mass);

    EquationData::Pressure<dim> pres(t_0);
    VectorTools::interpolate(dof_handler_pressure, pres, pres_n_minus_1);
    pres.advance_time(dt);
    VectorTools::interpolate(dof_handler_pressure, pres, pres_n);
    phi_n         = 0.;
    phi_n_minus_1 = 0.;
    for (unsigned int d = 0; d < dim; ++d)
      {
        vel_exact.set_time(t_0);
        vel_exact.set_component(d);
        VectorTools::interpolate(dof_handler_velocity,
                                 vel_exact,
                                 u_n_minus_1[d]);
        vel_exact.advance_time(dt);
        VectorTools::interpolate(dof_handler_velocity, vel_exact, u_n[d]);
      }
  }

  template <int dim>
  void NavierStokesProjection<dim>::initialize_velocity_matrices()
  {
    {
      DynamicSparsityPattern dsp(dof_handler_velocity.n_dofs(),
                                 dof_handler_velocity.n_dofs());
      DoFTools::make_sparsity_pattern(dof_handler_velocity, dsp);
      sparsity_pattern_velocity.copy_from(dsp);
    }
    vel_Laplace_plus_Mass.reinit(sparsity_pattern_velocity);
    for (unsigned int d = 0; d < dim; ++d)
      vel_it_matrix[d].reinit(sparsity_pattern_velocity);
    vel_Mass.reinit(sparsity_pattern_velocity);
    vel_Laplace.reinit(sparsity_pattern_velocity);
    vel_Advection.reinit(sparsity_pattern_velocity);

    MatrixCreator::create_mass_matrix(dof_handler_velocity,
                                      quadrature_velocity,
                                      vel_Mass);
    MatrixCreator::create_laplace_matrix(dof_handler_velocity,
                                         quadrature_velocity,
                                         vel_Laplace);
  }

  template <int dim>
  void NavierStokesProjection<dim>::initialize_pressure_matrices()
  {
    {
      DynamicSparsityPattern dsp(dof_handler_pressure.n_dofs(),
                                 dof_handler_pressure.n_dofs());
      DoFTools::make_sparsity_pattern(dof_handler_pressure, dsp);
      sparsity_pattern_pressure.copy_from(dsp);
    }

    pres_Laplace.reinit(sparsity_pattern_pressure);
    pres_iterative.reinit(sparsity_pattern_pressure);
    pres_Mass.reinit(sparsity_pattern_pressure);

    MatrixCreator::create_laplace_matrix(dof_handler_pressure,
                                         quadrature_pressure,
                                         pres_Laplace);
    MatrixCreator::create_mass_matrix(dof_handler_pressure,
                                      quadrature_pressure,
                                      pres_Mass);
  }

  // For the gradient operator, we start by initializing the sparsity pattern
  // and compressing it. It is important to notice here that the gradient
  // operator acts from the pressure space into the velocity space, so we have
  // to deal with two different finite element spaces. To keep the loops
  // synchronized, we use the alias that we have defined before, namely
  // <code>PairedIterators</code> and <code>IteratorPair</code>.
  template <int dim>
  void NavierStokesProjection<dim>::initialize_gradient_operator()
  {
    {
      DynamicSparsityPattern dsp(dof_handler_velocity.n_dofs(),
                                 dof_handler_pressure.n_dofs());
      DoFTools::make_sparsity_pattern(dof_handler_velocity,
                                      dof_handler_pressure,
                                      dsp);
      sparsity_pattern_pres_vel.copy_from(dsp);
    }

    InitGradPerTaskData per_task_data(0,
                                      fe_velocity.n_dofs_per_cell(),
                                      fe_pressure.n_dofs_per_cell());
    InitGradScratchData scratch_data(fe_velocity,
                                     fe_pressure,
                                     quadrature_velocity,
                                     update_gradients | update_JxW_values,
                                     update_values);

    for (unsigned int d = 0; d < dim; ++d)
      {
        pres_Diff[d].reinit(sparsity_pattern_pres_vel);
        per_task_data.d = d;
        WorkStream::run(
          IteratorPair(IteratorTuple(dof_handler_velocity.begin_active(),
                                     dof_handler_pressure.begin_active())),
          IteratorPair(IteratorTuple(dof_handler_velocity.end(),
                                     dof_handler_pressure.end())),
          *this,
          &NavierStokesProjection<dim>::assemble_one_cell_of_gradient,
          &NavierStokesProjection<dim>::copy_gradient_local_to_global,
          scratch_data,
          per_task_data);
      }
  }

  template <int dim>
  void NavierStokesProjection<dim>::assemble_one_cell_of_gradient(
    const IteratorPair  &SI,
    InitGradScratchData &scratch,
    InitGradPerTaskData &data)
  {
    scratch.fe_val_vel.reinit(std::get<0>(*SI));
    scratch.fe_val_pres.reinit(std::get<1>(*SI));

    std::get<0>(*SI)->get_dof_indices(data.vel_local_dof_indices);
    std::get<1>(*SI)->get_dof_indices(data.pres_local_dof_indices);

    data.local_grad = 0.;
    for (unsigned int q = 0; q < scratch.nqp; ++q)
      {
        for (unsigned int i = 0; i < data.vel_dpc; ++i)
          for (unsigned int j = 0; j < data.pres_dpc; ++j)
            data.local_grad(i, j) +=
              -scratch.fe_val_vel.JxW(q) *
              scratch.fe_val_vel.shape_grad(i, q)[data.d] *
              scratch.fe_val_pres.shape_value(j, q);
      }
  }

  template <int dim>
  void NavierStokesProjection<dim>::copy_gradient_local_to_global(
    const InitGradPerTaskData &data)
  {
    for (unsigned int i = 0; i < data.vel_dpc; ++i)
      for (unsigned int j = 0; j < data.pres_dpc; ++j)
        pres_Diff[data.d].add(data.vel_local_dof_indices[i],
                              data.pres_local_dof_indices[j],
                              data.local_grad(i, j));
  }

  // Function to estimate the time step
  template <int dim>
  double NavierStokesProjection<dim>::estimate_time_step(const Vector<double> u_n[dim], const Vector<double> u_n_minus_1[dim], const double &dt)
  {
    double TOL = 0.1; // Tolerance value
    double alpha = 0.5; // Parameter in the range (0,1]

    // Compute the time error indicator and the normalization factor
    double time_error_indicator = 0.;
    double normalization = 0.;

    FEValues<dim> fe_values(fe_velocity, quadrature_velocity, update_gradients | update_JxW_values);
    
    const unsigned int n_q_points = quadrature_velocity.size();
    std::vector<Tensor<1, dim>> gradients_u_diff(n_q_points);
    std::vector<Tensor<1, dim>> gradients_u_n(n_q_points);

    for (unsigned int d = 0; d < dim; ++d)
    {
      for (const auto &cell : dof_handler_velocity.active_cell_iterators())
      {
        fe_values.reinit(cell);
        
        Vector<double> u_diff = u_n[d];
        u_diff -= u_n_minus_1[d];

        fe_values.get_function_gradients(u_diff, gradients_u_diff);
        fe_values.get_function_gradients(u_n[d], gradients_u_n);

        for (unsigned int q = 0; q < n_q_points; ++q)
        {
          time_error_indicator += gradients_u_diff[q] * gradients_u_diff[q] * fe_values.JxW(q);
          normalization += gradients_u_n[q] * gradients_u_n[q] * fe_values.JxW(q);
        }
      }
    }

    // Multiply by dt to account for temporal integration
    time_error_indicator *= dt;
    normalization *= dt;

    // Debugging information
    //std::cout << "Time error indicator = " << time_error_indicator << std::endl;
    //std::cout << "Normalization = " << normalization << std::endl;
      
    // Compute the next time step value based on the error indicator
    if (time_error_indicator > (normalization*((1.+alpha)*(1.+alpha))*TOL*TOL))
      {
        double den1 = std::sqrt(normalization)*(TOL*TOL)*(1.+alpha/2)*(1.+alpha/2);
        return std::max(dt / std::min(2., time_error_indicator/den1), 0.0001);
      }
    else if (time_error_indicator < (normalization*((1.-alpha)*(1.-alpha))*TOL*TOL))
      {
        double den2 = std::sqrt(normalization)*(TOL*TOL)*(1.-alpha/2)*(1.-alpha/2);
        return std::min(dt / std::max(0.5, time_error_indicator/den2), 0.1);
      }
    else
      return dt;
  }

  // Function to save the time step values in a csv file
  template <int dim>
  void NavierStokesProjection<dim>::save_time_steps_to_csv(const std::string &filename) const
  {
    std::ofstream file(filename);
    if (file.is_open())
      {
        file << "Step" << "," << "Time step" << "," << "Time" << std::endl;
        unsigned int step = 1;
        double time = t_0;
        for (const auto &dt : dt_values)
        {
          time += dt;
          file << step << "," << dt << "," << time << std::endl;
          ++step;
        }
        file.close();
      }
    else
    {
      std::cerr << "Unable to open file" << std::endl;
    }
  }

  // This is the time marching function
  template <int dim>
  void NavierStokesProjection<dim>::run(const bool         verbose,
                                        const unsigned int output_interval)
  {
    ConditionalOStream verbose_cout(std::cout, verbose);

    double time = t_0; // Initialize the time
    unsigned int step = 1;
    vel_exact.set_time(2. * dt);
    output_results(step);
    dt_values.push_back(dt); // Store the initial time step value
    time += dt; // Update the time
    while (time <= T)
      {
        ++step; // Update the step number
        time += dt; // Update the time

        if (step % output_interval == 0)
          {
            verbose_cout << "Plotting Solution" << std::endl;
            output_results(step);
          }
        std::cout << "Step = " << step << " Time = " << time << std::endl;
        verbose_cout << "  Interpolating the velocity " << std::endl;
        interpolate_velocity();

        verbose_cout << "  Diffusion Step" << std::endl;
        if (step % vel_update_prec == 0)
          verbose_cout << "    With reinitialization of the preconditioner"
                       << std::endl;
        diffusion_step((step % vel_update_prec == 0) || (step == 2));

        verbose_cout << "  Projection Step" << std::endl;
        projection_step((step == 2));

        verbose_cout << "  Updating the Pressure" << std::endl;
        update_pressure((step == 2));
        
        vel_exact.advance_time(dt);

        // Estimate the new time step
        double new_dt = estimate_time_step(u_n, u_n_minus_1, dt);
        dt = new_dt;
        dt_values.push_back(dt); // Store the time step value

        // Print the new time step value
        std::cout << "New time step = " << dt << std::endl;
        std::cout << "------------------------------" << std::endl;

        // Update matrix vel_Laplace_plus_Mass (since dt has changed)
        vel_Laplace_plus_Mass = 0.;
        vel_Laplace_plus_Mass.add(1. / Re, vel_Laplace);
        vel_Laplace_plus_Mass.add(1.5 / dt, vel_Mass);
      }
    output_results(step);

    // Save the time step values in a csv file
    save_time_steps_to_csv("time_steps.csv");
  }

  template <int dim>
  void NavierStokesProjection<dim>::interpolate_velocity()
  {
    for (unsigned int d = 0; d < dim; ++d)
      {
        u_star[d].equ(2., u_n[d]);
        u_star[d] -= u_n_minus_1[d];
      }
  }

  // The implementation of a diffusion step.
  template <int dim>
  void NavierStokesProjection<dim>::diffusion_step(const bool reinit_prec)
  {
    pres_tmp.equ(-1., pres_n);
    pres_tmp.add(-4. / 3., phi_n, 1. / 3., phi_n_minus_1);

    assemble_advection_term();

    for (unsigned int d = 0; d < dim; ++d)
      {
        force[d] = 0.;
        v_tmp.equ(2. / dt, u_n[d]);
        v_tmp.add(-.5 / dt, u_n_minus_1[d]);
        vel_Mass.vmult_add(force[d], v_tmp);

        pres_Diff[d].vmult_add(force[d], pres_tmp);
        u_n_minus_1[d] = u_n[d];

        vel_it_matrix[d].copy_from(vel_Laplace_plus_Mass);
        vel_it_matrix[d].add(1., vel_Advection);

        vel_exact.set_component(d);
        boundary_values.clear();
        for (const auto &boundary_id : boundary_ids)
          {
            switch (boundary_id)
              {
                case 1: // 3
                  VectorTools::interpolate_boundary_values(
                    dof_handler_velocity,
                    boundary_id,
                    Functions::ZeroFunction<dim>(),
                    boundary_values);
                  break;
                case 2: // 0
                  VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                                           boundary_id,
                                                           vel_exact,
                                                           boundary_values);
                  break;
                case 3: // 1
                  if (d != 0)
                    VectorTools::interpolate_boundary_values(
                      dof_handler_velocity,
                      boundary_id,
                      Functions::ZeroFunction<dim>(),
                      boundary_values);
                  break;
                case 4: // 2
                  VectorTools::interpolate_boundary_values(
                    dof_handler_velocity,
                    boundary_id,
                    Functions::ZeroFunction<dim>(),
                    boundary_values);
                  break;
                default:
                  DEAL_II_NOT_IMPLEMENTED();
              }
          }
        MatrixTools::apply_boundary_values(boundary_values,
                                           vel_it_matrix[d],
                                           u_n[d],
                                           force[d]);
      }

    Threads::TaskGroup<void> tasks;
    for (unsigned int d = 0; d < dim; ++d)
      {
        if (reinit_prec)
          prec_velocity[d].initialize(vel_it_matrix[d],
                                      SparseILU<double>::AdditionalData(
                                        vel_diag_strength, vel_off_diagonals));
        tasks += Threads::new_task(
          &NavierStokesProjection<dim>::diffusion_component_solve, *this, d);
      }
    tasks.join_all();
  }

  template <int dim>
  void
  NavierStokesProjection<dim>::diffusion_component_solve(const unsigned int d)
  {
    SolverControl solver_control(vel_max_its, vel_eps * force[d].l2_norm());
    SolverGMRES<Vector<double>> gmres(
      solver_control,
      SolverGMRES<Vector<double>>::AdditionalData(vel_Krylov_size));
    gmres.solve(vel_it_matrix[d], u_n[d], force[d], prec_velocity[d]);
  }

  // The following few functions deal with assembling the advection terms.
  template <int dim>
  void NavierStokesProjection<dim>::assemble_advection_term()
  {
    vel_Advection = 0.;
    AdvectionPerTaskData data(fe_velocity.n_dofs_per_cell());
    AdvectionScratchData scratch(fe_velocity,
                                 quadrature_velocity,
                                 update_values | update_JxW_values |
                                   update_gradients);
    WorkStream::run(
      dof_handler_velocity.begin_active(),
      dof_handler_velocity.end(),
      *this,
      &NavierStokesProjection<dim>::assemble_one_cell_of_advection,
      &NavierStokesProjection<dim>::copy_advection_local_to_global,
      scratch,
      data);
  }

  template <int dim>
  void NavierStokesProjection<dim>::assemble_one_cell_of_advection(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    AdvectionScratchData                                 &scratch,
    AdvectionPerTaskData                                 &data)
  {
    scratch.fe_val.reinit(cell);
    cell->get_dof_indices(data.local_dof_indices);
    for (unsigned int d = 0; d < dim; ++d)
      {
        scratch.fe_val.get_function_values(u_star[d], scratch.u_star_tmp);
        for (unsigned int q = 0; q < scratch.nqp; ++q)
          scratch.u_star_local[q][d] = scratch.u_star_tmp[q];
      }

    for (unsigned int d = 0; d < dim; ++d)
      {
        scratch.fe_val.get_function_gradients(u_star[d], scratch.grad_u_star);
        for (unsigned int q = 0; q < scratch.nqp; ++q)
          {
            if (d == 0)
              scratch.u_star_tmp[q] = 0.;
            scratch.u_star_tmp[q] += scratch.grad_u_star[q][d];
          }
      }

    data.local_advection = 0.;
    for (unsigned int q = 0; q < scratch.nqp; ++q)
      for (unsigned int i = 0; i < scratch.dpc; ++i)
        for (unsigned int j = 0; j < scratch.dpc; ++j)
          data.local_advection(i, j) += (scratch.u_star_local[q] *            //
                                           scratch.fe_val.shape_grad(j, q) *  //
                                           scratch.fe_val.shape_value(i, q)   //
                                         +                                    //
                                         0.5 *                                //
                                           scratch.u_star_tmp[q] *            //
                                           scratch.fe_val.shape_value(i, q) * //
                                           scratch.fe_val.shape_value(j, q))  //
                                        * scratch.fe_val.JxW(q);
  }

  template <int dim>
  void NavierStokesProjection<dim>::copy_advection_local_to_global(
    const AdvectionPerTaskData &data)
  {
    for (unsigned int i = 0; i < fe_velocity.n_dofs_per_cell(); ++i)
      for (unsigned int j = 0; j < fe_velocity.n_dofs_per_cell(); ++j)
        vel_Advection.add(data.local_dof_indices[i],
                          data.local_dof_indices[j],
                          data.local_advection(i, j));
  }

  // This implements the projection step:
  template <int dim>
  void NavierStokesProjection<dim>::projection_step(const bool reinit_prec)
  {
    pres_iterative.copy_from(pres_Laplace);

    pres_tmp = 0.;
    for (unsigned d = 0; d < dim; ++d)
      pres_Diff[d].Tvmult_add(pres_tmp, u_n[d]);

    phi_n_minus_1 = phi_n;

    static std::map<types::global_dof_index, double> bval;
    if (reinit_prec)
      VectorTools::interpolate_boundary_values(dof_handler_pressure,
                                               3,
                                               Functions::ZeroFunction<dim>(),
                                               bval);

    MatrixTools::apply_boundary_values(bval, pres_iterative, phi_n, pres_tmp);

    if (reinit_prec)
      prec_pres_Laplace.initialize(pres_iterative,
                                   SparseILU<double>::AdditionalData(
                                     vel_diag_strength, vel_off_diagonals));

    SolverControl solvercontrol(vel_max_its, vel_eps * pres_tmp.l2_norm());
    SolverCG<Vector<double>> cg(solvercontrol);
    cg.solve(pres_iterative, phi_n, pres_tmp, prec_pres_Laplace);

    phi_n *= 1.5 / dt;
  }

  // This is the pressure update step of the projection method. It implements
  // the standard formulation of the method, that is @f[ p^{n+1} = p^n +
  // \phi^{n+1}, @f] or the rotational form, which is @f[ p^{n+1} = p^n +
  // \phi^{n+1} - \frac{1}{Re} \nabla\cdot u^{n+1}. @f]
  template <int dim>
  void NavierStokesProjection<dim>::update_pressure(const bool reinit_prec)
  {
    pres_n_minus_1 = pres_n;
    switch (type)
      {
        case RunTimeParameters::Method::standard:
          pres_n += phi_n;
          break;
        case RunTimeParameters::Method::rotational:
          if (reinit_prec)
            prec_mass.initialize(pres_Mass);
          pres_n = pres_tmp;
          prec_mass.solve(pres_n);
          pres_n.sadd(1. / Re, 1., pres_n_minus_1);
          pres_n += phi_n;
          break;
        default:
          DEAL_II_NOT_IMPLEMENTED();
      };
  }

  // This method plots the current solution. The main difficulty is that we want
  // to create a single output file that contains the data for all velocity
  // components, the pressure, and also the vorticity of the flow. On the other
  // hand, velocities and the pressure live on separate DoFHandler objects, and
  // so can't be written to the same file using a single DataOut object.
  template <int dim>
  void NavierStokesProjection<dim>::output_results(const unsigned int step)
  {
    assemble_vorticity((step == 1));
    const FESystem<dim> joint_fe(fe_velocity ^ dim, fe_pressure, fe_velocity);
    DoFHandler<dim>     joint_dof_handler(triangulation);
    joint_dof_handler.distribute_dofs(joint_fe);
    Assert(joint_dof_handler.n_dofs() ==
             ((dim + 1) * dof_handler_velocity.n_dofs() +
              dof_handler_pressure.n_dofs()),
           ExcInternalError());
    Vector<double> joint_solution(joint_dof_handler.n_dofs());
    std::vector<types::global_dof_index> loc_joint_dof_indices(
      joint_fe.n_dofs_per_cell()),
      loc_vel_dof_indices(fe_velocity.n_dofs_per_cell()),
      loc_pres_dof_indices(fe_pressure.n_dofs_per_cell());
    typename DoFHandler<dim>::active_cell_iterator
      joint_cell = joint_dof_handler.begin_active(),
      joint_endc = joint_dof_handler.end(),
      vel_cell   = dof_handler_velocity.begin_active(),
      pres_cell  = dof_handler_pressure.begin_active();
    for (; joint_cell != joint_endc; ++joint_cell, ++vel_cell, ++pres_cell)
      {
        joint_cell->get_dof_indices(loc_joint_dof_indices);
        vel_cell->get_dof_indices(loc_vel_dof_indices);
        pres_cell->get_dof_indices(loc_pres_dof_indices);
        for (unsigned int i = 0; i < joint_fe.n_dofs_per_cell(); ++i)
          switch (joint_fe.system_to_base_index(i).first.first)
            {
              case 0:
                Assert(joint_fe.system_to_base_index(i).first.second < dim,
                       ExcInternalError());
                joint_solution(loc_joint_dof_indices[i]) =
                  u_n[joint_fe.system_to_base_index(i).first.second](
                    loc_vel_dof_indices[joint_fe.system_to_base_index(i)
                                          .second]);
                break;
              case 1:
                Assert(joint_fe.system_to_base_index(i).first.second == 0,
                       ExcInternalError());
                joint_solution(loc_joint_dof_indices[i]) =
                  pres_n(loc_pres_dof_indices[joint_fe.system_to_base_index(i)
                                                .second]);
                break;
              case 2:
                Assert(joint_fe.system_to_base_index(i).first.second == 0,
                       ExcInternalError());
                joint_solution(loc_joint_dof_indices[i]) = rot_u(
                  loc_vel_dof_indices[joint_fe.system_to_base_index(i).second]);
                break;
              default:
                DEAL_II_ASSERT_UNREACHABLE();
            }
      }
    std::vector<std::string> joint_solution_names(dim, "v");
    joint_solution_names.emplace_back("p");
    joint_solution_names.emplace_back("rot_u");
    DataOut<dim> data_out;
    data_out.attach_dof_handler(joint_dof_handler);
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      component_interpretation(
        dim + 2, DataComponentInterpretation::component_is_part_of_vector);
    component_interpretation[dim] =
      DataComponentInterpretation::component_is_scalar;
    component_interpretation[dim + 1] =
      DataComponentInterpretation::component_is_scalar;
    data_out.add_data_vector(joint_solution,
                             joint_solution_names,
                             DataOut<dim>::type_dof_data,
                             component_interpretation);
    data_out.build_patches(deg + 1);
    std::ofstream output("solution-" + Utilities::int_to_string(step, 5) +
                         ".vtk");
    data_out.write_vtk(output);
  }



  // Following is the helper function that computes the vorticity by projecting
  // the term $\text{curl} u$ onto the finite element space used for the
  // components of the velocity.
  template <int dim>
  void NavierStokesProjection<dim>::assemble_vorticity(const bool reinit_prec)
  {
    Assert(dim == 2, ExcNotImplemented());
    if (reinit_prec)
      prec_vel_mass.initialize(vel_Mass);

    FEValues<dim>      fe_val_vel(fe_velocity,
                             quadrature_velocity,
                             update_gradients | update_JxW_values |
                               update_values);
    const unsigned int dpc = fe_velocity.n_dofs_per_cell(),
                       nqp = quadrature_velocity.size();
    std::vector<types::global_dof_index> ldi(dpc);
    Vector<double>                       loc_rot(dpc);

    std::vector<Tensor<1, dim>> grad_u1(nqp), grad_u2(nqp);
    rot_u = 0.;

    for (const auto &cell : dof_handler_velocity.active_cell_iterators())
      {
        fe_val_vel.reinit(cell);
        cell->get_dof_indices(ldi);
        fe_val_vel.get_function_gradients(u_n[0], grad_u1);
        fe_val_vel.get_function_gradients(u_n[1], grad_u2);
        loc_rot = 0.;
        for (unsigned int q = 0; q < nqp; ++q)
          for (unsigned int i = 0; i < dpc; ++i)
            loc_rot(i) += (grad_u2[q][0] - grad_u1[q][1]) * //
                          fe_val_vel.shape_value(i, q) *    //
                          fe_val_vel.JxW(q);

        for (unsigned int i = 0; i < dpc; ++i)
          rot_u(ldi[i]) += loc_rot(i);
      }

    prec_vel_mass.solve(rot_u);
  }
} // namespace Step35

int main()
{
  try
    {
      using namespace Step35;

      RunTimeParameters::Data_Storage data;
      data.read_data("parameter-file.prm");

      deallog.depth_console(data.verbose ? 2 : 0);

      std::cout << "Adaptive time-step version of the code" << std::endl;
      std::cout << "--------------------------------------" << std::endl;

      NavierStokesProjection<2> test(data);
      test.run(data.verbose, data.output_interval);
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  std::cout << "----------------------------------------------------"
            << std::endl
            << "Apparently everything went fine!" << std::endl
            << "Now have a look at the results :-)" << std::endl
            << std::endl;
  return 0;
}
